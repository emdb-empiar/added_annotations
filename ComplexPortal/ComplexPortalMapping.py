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

pdb_sifts_ftp = r'/Users/amudha/project/uniprot_pdb.csv'
#pdb_sifts_ftp = r'/nfs/ftp/pub/databases/msd/sifts/cvs/uniprot_pdb.csv'
#CP_ftp = r'/nfs/ftp/pub/databases/intact/complex/current/complextab/'
CP_ftp = "/Users/amudha/project/cpx_data/complextab/"
pdb_baseurl = r'https://www.ebi.ac.uk/pdbe/graph-api/mappings/uniprot/'
BLAST_DB = "uniprot_sprot"  # uniprotkb_swissprot

class PDBeCPX:
    """
    PDBeCPX Mapping extracted from complex_portal_output_complete_complexes.csv
    """

    def __init__(self, row):
        self.pdbe_complex_id = row[0]
        self.cpx_id = row[1]
        self.components = row[2]
        self.pdb_list = row[3].split(',')
        self.emdb_list = row[5].split(',')
        self.taxonomy = row[6]

        self.pdb_list = [x for x in self.pdb_list if x != '']
        self.emdb_list = [x for x in self.emdb_list if x != '']

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s" % (
        self.pdbe_complex_id, self.cpx_id, self.components, str(self.pdb_list), str(self.emdb_list), self.taxonomy)

    def get_annotations(self):
        annotations = []
        for emdb_id in self.emdb_list:
            annotation = Annotation(emdb_id, self.pdbe_complex_id, self.cpx_id, "PDB_ID")
            annotations.append(annotation)
            logger.debug(annotation)
        return annotations

    def get_pdb_relations(self):
        relations = {}
        for pdb_id in self.pdb_list:
            relations[pdb_id.lower()] = self.cpx_id
        return relations


class MolCPX:
    """
    PDBeCPX Mapping extracted from pdb_complex_protein_details_complete_complexes.csv
    """

    def __init__(self, row):
        self.database = row[0]
        self.accession = row[1]
        self.mol_name = row[2]
        self.stoichiometry = row[3]
        self.polymer_type = row[4]
        self.taxonomy = row[5]
        self.cpx_id = row[6]

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s" % (
        self.database, self.accession, self.mol_name, self.polymer_type, self.cpx_id, self.taxonomy)

    def is_uniprot(self):
        return self.database == 'UNP'


class CPX:
    """
    Complex Portal entry obtained from their FTP area
    """

    def __init__(self, row):
        self.cpx_id = row[0]
        self.name = row[1]
        self.taxonomy = row[3]
        self.identifiers = re.sub(r'\(\d+\)', '', row[4]).split('|')
        self.confidence = row[5]
        self.GO = re.sub(r'\(.+?\)', '', row[7]).split('|')
        self.cross_ref = re.findall(r':(.*?)\(', row[8], re.S)

class CPX_database:
    """
    Class containing all the CPX entries from their FTP area
    """

    def __init__(self):
        self.entries = {}
        self.id_mapped = {}
        self.emdb_id_mapped = {}

    def add_row(self, row):
        cpx = CPX(row)
        self.entries[cpx.cpx_id] = cpx
        for identifier in cpx.identifiers:
            self.id_mapped[identifier] = cpx.cpx_id
        for cross_ref in cpx.cross_ref:
            self.emdb_id_mapped[cross_ref] = cpx.cpx_id

    def get_from_cpx(self, cpx_id):
        if cpx_id in self.entries:
            return self.entries[cpx_id]
        return None

    # Use this method to return CPX entry based on Uniprot, CHEBI or PubMed ids
    def get_from_identifier(self, ext_id):
        if ext_id in self.id_mapped:
            return self.id_mapped[ext_id]
        return None

    # Use this method to return CPX entry based on EMDB_ID
    def get_from_cross_ref(self, emd_id):
        if emd_id in self.emdb_id_mapped:
            return self.emdb_id_mapped[emd_id]
        return None


class Annotation:
    """
    EMDB-CPX annotation
    """

    def __init__(self, emdb_id, group_by, cpx_id, method, cpx_title=""):
        self.emdb_id = emdb_id
        self.group_by = group_by
        self.cpx_id = cpx_id
        self.cpx_title = cpx_title  # Maybe remove or add in the end
        self.method = method

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s" % (self.emdb_id, self.group_by, self.cpx_id, self.cpx_title, self.method)

class CPMapping:
    """
    Extracting the IDs from the EMDB entry and for map only entries running the BLASTP search to get UNIPROT IDs.
    Querying with the extracted IDs for mapping the EMDB entries to the complex portal if exists.
    """

    def __init__(self, workDir, headerDir, PDBeDir):
        self.cpx_db = CPX_database()
        self.annotations = []
        self.workDir = workDir
        self.headerDir = headerDir
        self.pdb_relations = {}  # Pdb_id -> cpx_id
        self.cpx_relations = {}  # pdb_cpx_id -> cpx_id
        self.uniprot_relations = {}  # Unp_id -> pdb_cpx_id
        self.uniprot_map = set()
        self.tax_ids = []

        self.pdbefile = os.path.join(PDBeDir, "complex_portal_output_complete_complexes.csv")
        self.pdbeSciFile = os.path.join(PDBeDir, "pdb_complex_protein_details_complete_complexes.csv")

        # Parse Complex Portal tables
        for fn in glob(os.path.join(str(CP_ftp), '*.tsv')):
            self.tax_ids.append(os.path.basename(fn).split('.')[0])
            with open(fn, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader, None)  # skip the headers
                batch_data = list(reader)
                for line in batch_data:
                    self.cpx_db.add_row(line)

        # Parse complex_portal_output_complete_complexes.csv
        with open(self.pdbefile, 'r') as f:
            reader = csv.reader(f)
            next(reader, None)  # skip the headers
            batch_data = list(reader)
            for line in batch_data:
                row = PDBeCPX(line)
                self.annotations += row.get_annotations()
                # Merge two dictionaries, requires python >= 3.5
                self.pdb_relations = {**self.pdb_relations, **row.get_pdb_relations()}
                self.cpx_relations[row.pdbe_complex_id] = row.cpx_id

        # Parse pdb_complex_protein_details_complete_complexes.csv
        #
        # I had to remove this block for two reasons:
        # 1: There was a bug mapping PDB-COMPLEX_ID instead of COMPLEX-ID
        # 2: This was mapping Uniprot->Complex Portal, therefore it will add complexes
        # unrelated to the structure.
        #
        # with open(self.pdbeSciFile, 'r') as f:
        #     reader = csv.reader(f)
        #     next(reader, None)  # skip the headers
        #     batch_data = list(reader)
        #     for line in batch_data:
        #         row = MolCPX(line)
        #         if row.is_uniprot():
        #             unp_id = row.accession
        #             if unp_id in self.uniprot_relations:
        #                 self.uniprot_relations[unp_id].add(row.cpx_id)
        #             else:
        #                 self.uniprot_relations[unp_id] = set([row.cpx_id])

    def execute(self):
        ###### Fetch header files for query ########
        for fn in glob(os.path.join(str(self.headerDir), '*')):
            print(fn)
            id_num = fn.split('-')[1]
            xml_filename = "emd-" + id_num + "-v30.xml"
            xml_dirpath = os.path.join(str(self.headerDir), fn, "header")
            xml_filepath = os.path.join(xml_dirpath, xml_filename)

            ####### Extract the ids and store author provided annotations #######
            uniprot_ids, pdb_ids = self.extracting_IDs(xml_filepath)

            ####### Fetching Uniprot from PDB ids #######
            for pdb_id in pdb_ids:
                ###self.extracting_UniprotFromPDBe("EMD-" + id_num, pdb_id) ## PDBe API call
                self.extracting_Uniprot_PDBeSifts("EMD-" + id_num, pdb_id)

            ####### CREATING FASTA SEQUENCE AND RUNNING BLAST ########
            self.uniprot_map = self.uniprot_map.union(self.create_fasta_run_blastp(id_num, xml_filepath))

            ####### CONVERTING UNIPROT DATA TO CPX ########
            self.map_uniprot_to_cpx("EMD-" + id_num, uniprot_ids)

            ####### CONVERTING EMDB DATA TO CPX ########
            self.map_emdb_to_cpx("EMD-" + id_num)

        self.parse_annotations()

    def extracting_IDs(self, xml_filepath):
        """
        Extract the IDs (EMDB, PDB from both header file and PDBe files and UNIPROT)
        """
        uniprot_ids = set()
        pdb_ids = set()

        with open(xml_filepath, 'r') as filexml:
            tree = ET.parse(filexml)
            root = tree.getroot()
            a = root.attrib
            emd_id = a.get('emdb_id')
            for x in list(root.iter('pdb_reference')):
                qs = x.find('pdb_id').text.lower()
                pdb_ids.add(qs)
                if qs in self.pdb_relations:
                    cpx_id = self.pdb_relations[qs]
                    self.annotations.append(Annotation(emd_id, qs, cpx_id, "PDBe"))
            if list(root.iter('protein_or_peptide')):
                for x in list(root.iter('protein_or_peptide')):
                    qs = x.find('sequence')
                    if qs.find('external_references') is not None:
                        if qs.find('external_references').attrib['type'] == 'UNIPROTKB':
                            q = qs.find('external_references').text
                            uniprot_ids.add(q)
                            self.uniprot_map.add((emd_id, q, "AUTHOR"))
        return uniprot_ids, pdb_ids

    def create_fasta_run_blastp(self, emd_id, xml_filepath):
        """
        Creating the fasta format file to run BLASTP
        """
        filexml = open(xml_filepath, 'r')
        tree = ET.parse(filexml)
        root = tree.getroot()
        uniprot_map = set()
        if list(root.iter('protein_or_peptide')):
            for x in list(root.iter('protein_or_peptide')):
                qs = x.find('sequence')
                if qs.find('string') is not None:
                    fasta_file = os.path.join(self.workDir, emd_id + ".fasta")
                    with open(fasta_file, "w") as f:
                        seq = qs.find('string').text
                        seq = re.sub(r'\(\s*UNK\s*\)', 'X', seq)
                        seq = seq.replace("\n", "")
                        f.write(">seq\n%s" % seq)

                    db_path = os.path.join(self.workDir, BLAST_DB)  #
                    qout = os.path.join(self.workDir, emd_id + ".out")
                    #### Using biopython for Blastp ##
                    #blastp_command = NcbiblastpCommandline(query=fasta_file, db=db_path, out=qout, evalue='1e-40')
                    #blastp_command()
                    blastp_command = ["blastp", "-query", fasta_file, "-db", db_path, "-out", qout, "-evalue", "1e-40"]
                    subprocess.call(blastp_command)
                    if os.path.isfile(qout):
                        uniprot_id = self.extract_uniprot_from_blast(qout)
                        if uniprot_id:
                            uniprot_map.add(("EMD-" + emd_id, uniprot_id, "BLAST"))

        return uniprot_map

    def extract_uniprot_from_blast(self, out_file):
        """
        Extracting the UNIPROT ID from BLASTP output
        """
        qresult = SearchIO.read(out_file, 'blast-text')
        seq_len = qresult.seq_len
        for result in qresult:
            desc = result._description
            uniprot_id = result.id.replace("SP:", "")
            tax_id = re.findall('OX=(.*)=', desc, re.S)[0].split("GN")[0].strip()
            if tax_id in self.tax_ids:
                for hit in result:
                    ident_num = hit.ident_num
                    coverage = ident_num / seq_len
                    if coverage < 1.0:
                        continue
                    return uniprot_id
        return None

    def extracting_UniprotFromPDBe(self, emdb_id, pdb_id):
        """
        Extracting the UNIPROT ID from PDBe API if model exists for the entry
        """
        url = pdb_baseurl + pdb_id

        def save_sslcontext(obj):
            return obj.__class__, (obj.protocol,)

        copyreg.pickle(ssl.SSLContext, save_sslcontext)
        context = ssl.create_default_context()
        foo = pickle.dumps(context)
        gcontext = pickle.loads(foo)
        try:
            pdbjson = urlopen(url, context=gcontext).read()
            pdbjdata = json.loads(pdbjson.decode('utf-8'))
            jdata = pdbjdata[pdb_id]['UniProt']
            for key in jdata.keys():
                self.uniprot_map.add((emdb_id, key, "UNIPROT_from_PDBe"))
        except Exception as e:
            logger.WARNING(e)

    def extracting_Uniprot_PDBeSifts(self, emdb_id, pdb_id):
        """
        Extracting the UNIPROT ID from PDBe SIFTS files if model exists for the entry
        """
        with open(pdb_sifts_ftp, 'r') as sifts_file:
            reader = csv.reader(sifts_file, delimiter=',')
            next(reader, None)
            for row in reader:
                for column in row:
                    if pdb_id in column:
                        self.uniprot_map.add((emdb_id, row[0], "UNIPROT_PDBe_SIFTS"))

    def parse_annotations(self):
        for emdb_id, uniprot, method in self.uniprot_map:
            if uniprot in self.uniprot_relations:
                cpx_id = self.uniprot_relations[uniprot]
                self.annotations.append(Annotation(emdb_id, uniprot, cpx_id, method))

    def map_uniprot_to_cpx(self, emdb_id, uniprot_ids):
        for uniprot in uniprot_ids:
            cpx = self.cpx_db.get_from_identifier(uniprot)
            if cpx:
                self.annotations.append(Annotation(emdb_id, uniprot, cpx, "UNIPROT"))
                logger.debug(Annotation(emdb_id, uniprot, cpx, "UNIPROT"))

    def map_emdb_to_cpx(self, emdb_id):
        cpx = self.cpx_db.get_from_cross_ref(emdb_id)
        if cpx:
            self.annotations.append(Annotation(emdb_id, emdb_id, cpx, "EMDB_ID"))
            logger.debug(Annotation(emdb_id, emdb_id, cpx, "EMDB_ID"))

    def write_cpx_map(self):
        filepath = os.path.join(self.workDir, "emdb_cpx.tsv")
        with open(filepath, 'w') as f:
            f.write("%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "QUERY_ID", "CPX_ID", "CPX_TITLE", "METHOD"))
            for cpx in self.annotations:
                cpx_obj = self.cpx_db.get_from_cpx(cpx.cpx_id)
                if cpx_obj:
                    cpx_title = cpx_obj.name
                else:
                    cpx_title = ""
                f.write("%s\t%s\t%s\t%s\t%s\n" % (cpx.emdb_id, cpx.group_by, cpx.cpx_id, cpx_title, cpx.method))

    def write_uniprot_map(self):
        filepath = os.path.join(self.workDir, "emdb_unp.tsv")
        with open(filepath, 'w') as f:
            f.write("%s\t%s\t%s\n" % ("EMDB_ID", "UNIPROT_ID", "QUERY_METHOD"))
            for emdb_id, uniprot, method in self.uniprot_map:
                f.write("%s\t%s\t%s\n" % (emdb_id, uniprot, method))