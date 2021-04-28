import os, re, csv, subprocess, pickle, copyreg, ssl, json
import lxml.etree as ET
from glob import glob
from Bio import SearchIO
from urllib.request import urlopen
#from Bio.Blast.Applications import NcbiblastpCommandline ## Biopython for BLASTP
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(funcName)s:%(message)s')

file_handler = logging.FileHandler('logging_ComplexPortal.log', mode ='w')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

CP_ftp = r'/nfs/ftp/pub/databases/intact/complex/current/complextab/'
#CP_ftp = "/Users/neli/EBI//annotations/data/cpx/"
#CP_ftp = "/Users/amudha/project/cpx_data/complextab/"
MIN_SCORE = 0.5

def overlap(set1, set2):
    L = len(set1.union(set2))
    match = len(set1.intersection(set2))
    return float(match)/float(L)

class CPX:
    """
    Complex Portal entry obtained from their FTP area
    """

    def __init__(self, row):
        self.cpx_id = row[0]
        self.name = row[1]
        self.taxonomy = row[3]
        self.uniprot = set()
        self.identifiers = re.sub(r'\(\d+\)', '', row[4]).split('|')
        self.confidence = row[5]
        self.GO = re.sub(r'\(.+?\)', '', row[7]).split('|')
        self.cross_ref = re.findall(r':(.*?)\(', row[8], re.S)

        for idt in self.identifiers:
            if 'CHEBI:' in idt:
                continue
            if '-PRO_' in idt:
                self.uniprot.add(idt.split('-')[0])
                continue
            if  '_' in idt:
                continue
            self.uniprot.add(idt)


class EMDB_complex:
    """
    EMDB complex sample obtained from the header files in the Uniprot mapping
    """

    def __init__(self, emdb_id, sample_id, complex_sample_id):
        self.emdb_id = emdb_id
        self.sample_id = sample_id
        self.complex_sample_id = complex_sample_id
        self.cpx_list = []
        self.proteins = set()
        self.method = ""
        self.score = 0.0

    def add_protein(self, uniprot_id):
        self.proteins.add(uniprot_id)


class CPX_database:
    """
    Class containing all the CPX entries from their FTP area
    """

    def __init__(self):
        self.entries = {}
        self.id_mapped = {}
        self.uniprot_map = {}

    def add_row(self, row):
        cpx = CPX(row)
        self.entries[cpx.cpx_id] = cpx
        for identifier in cpx.identifiers:
            if identifier in self.id_mapped:
                self.id_mapped[identifier].add(cpx.cpx_id)
            else:
                self.id_mapped[identifier] = set([cpx.cpx_id])
        for unp_id in cpx.uniprot:
            if unp_id in self.uniprot_map:
                self.uniprot_map[unp_id].add(cpx.cpx_id)
            else:
                self.uniprot_map[unp_id] = set([cpx.cpx_id])

    def get_from_cpx(self, cpx_id):
        if cpx_id in self.entries:
            return self.entries[cpx_id]
        return None

    # Use this method to return CPX entry based on Uniprot, CHEBI or PubMed ids
    def get_from_identifier(self, ext_id):
        if ext_id in self.id_mapped:
            return self.id_mapped[ext_id]
        return None

    # Use this method to return CPX entry based on Uniprot ID
    def get_from_uniprot(self, unp_id):
        if unp_id in self.uniprot_map:
            return self.uniprot_map[unp_id]
        return None

class CPMapping:
    """
    Extracting the IDs from the EMDB entry and for map only entries running the BLASTP search to get UNIPROT IDs.
    Querying with the extracted IDs for mapping the EMDB entries to the complex portal if exists.
    """

    def __init__(self, workDir, headerDir, PDBeDir, proteins):
        self.cpx_db = CPX_database()
        self.emdb_complexes = {}
        self.annotations = []
        self.workDir = workDir
        self.headerDir = headerDir
        self.tax_ids = []

        # Parse Complex Portal tables
        for fn in glob(os.path.join(str(CP_ftp), '*.tsv')):
            self.tax_ids.append(os.path.basename(fn).split('.')[0])
            with open(fn, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader, None)  # skip the headers
                batch_data = list(reader)
                for line in batch_data:
                    self.cpx_db.add_row(line)

        # Parse EMDB proteins
        for protein in proteins:
            #protein is part of a complex and contains a Uniprot ID
            if len(protein.sample_complexes) > 0 and protein.uniprot_id:
                emdb_id = protein.emdb_id
                sample_id = protein.sample_id
                uniprot_id = protein.uniprot_id
                for complex_id in protein.sample_complexes:
                    emdb_complex_id = emdb_id + "_" + str(complex_id)
                    if emdb_complex_id in self.emdb_complexes:
                        self.emdb_complexes[emdb_complex_id].add_protein(uniprot_id)
                    else:
                        emdb_cpx = EMDB_complex(emdb_id, complex_id, emdb_complex_id)
                        emdb_cpx.add_protein(uniprot_id)
                        self.emdb_complexes[emdb_complex_id] = emdb_cpx

    def execute(self):
        for emdb_complex in self.emdb_complexes.values():
            emdb_protein_list = emdb_complex.proteins
            cpx_found = set()
            for uniprot_id in emdb_protein_list:
                cpx_complexes = self.cpx_db.get_from_uniprot(uniprot_id)
                if cpx_complexes:
                    for cpx_id in cpx_complexes:
                        cpx_found.add(cpx_id)

            max_score = 0.0
            best_hits = []
            for cpx_id in cpx_found:
                cpx = self.cpx_db.get_from_cpx(cpx_id)
                cpx_uniprot = cpx.uniprot
                overlap_score = overlap(emdb_protein_list,cpx_uniprot)
                if overlap_score == max_score:
                    best_hits.append(cpx)
                elif overlap_score > max_score:
                    best_hits = [cpx]
                    max_score = overlap_score

            if max_score >= MIN_SCORE:
                emdb_complex.cpx_list = best_hits
                emdb_complex.score = max_score
                emdb_complex.method = "UNIPROT"

                self.annotations.append(emdb_complex)

    def write_cpx_map(self):
        filepath = os.path.join(self.workDir, "emdb_cpx.tsv")
        with open(filepath, 'w') as f:
            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "EMDB_SAMPLE_ID", "CPX_ID", "CPX_TITLE", "METHOD", "SCORE"))
            for emcpx in self.annotations:
                for cpx in emcpx.cpx_list:
                    f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (emcpx.emdb_id, emcpx.sample_id, cpx.cpx_id, cpx.name, emcpx.method, emcpx.score))




