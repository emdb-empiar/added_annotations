#!/usr/bin/env/python

import os, re, sys
import argparse
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError
import json
import ssl
from lxml import etree
import lxml.etree as ET
import subprocess
import operator
import glob
import pickle, copyreg
import csv

#CP_baseurl = r'https://www.ebi.ac.uk/intact/complex-ws/search/'
#pdb_baseurl = r'https://www.ebi.ac.uk/pdbe/graph-api/mappings/uniprot/'
#workDir = r'/homes/amudha/project/'
###headerDir = r'/nfs/nobackup/msd/emdep_test/EM_prepare2/structures/'
#headerDir = r'/nfs/nobackup/msd/pdb_em_test/EM_ftp/structures/'
#PDBeDir = r"/nfs/public/ro/pdbe/release-data/complex/"
###PDBeDir = r"/nfs/msd/work2/msdsd/complex_portal/"
###PDBeDir = r"/nfs/nobackup/msd/complexes/"
#
CP_baseurl = r'https://www.ebi.ac.uk/intact/complex-ws/search/'
pdb_baseurl = r'https://www.ebi.ac.uk/pdbe/graph-api/mappings/uniprot/'

class CPMapping:
    """
    Extracting the IDs from the EMDB entry and for map only entries running the BLASTP search to get UNIPROT IDs.
    Querying with the extracted IDs for mapping the EMDB entries to the complex portal if exists.
    """

    def __init__(self, workDir, headerDir, PDBeDir):
        ###### Path for the files ##########
        self.workDir = workDir
        self.mapFile = os.path.join(str(workDir) + "/ComplexPortal/" + "EMDB_CPX_map.tsv")
        self.sort_mapFile = os.path.join(str(workDir) + "/ComplexPortal/" + "EMDB_CPX_map_sorted.tsv")
        self.outFile_uniprot = os.path.join(str(workDir) + "/ComplexPortal/" + "Extracted_UNIPROT_PBDe_BLASTP.tsv")
        self.pdbefile = os.path.join(str(PDBeDir) + "/" + "complex_portal_output_complete_complexes.csv")
        self.pdbeSciFile = os.path.join(str(PDBeDir) + "/" + "pdb_complex_protein_details_complete_complexes.csv")
        self.taxidFile = os.path.join(str(workDir) + "/" + 'ncbi_taxid_scientific_names.txt')
        ##### Remove if files exists ##########
        if os.path.exists(self.outFile_uniprot):
            os.remove(self.outFile_uniprot)
        else:
            pass
        if os.path.exists(self.mapFile):
            os.remove(self.mapFile)
        if os.path.exists(self.sort_mapFile):
            os.remove(self.sort_mapFile)
        ###### Fetch header files for query ########
        for self.fn in glob.glob(os.path.join(str(headerDir), '*')):
                self.id_num = self.fn.split("-",1)[1]
                xml_filename = "emd-" + self.id_num + "-v30.xml"
                xml_dirpath = os.path.join(str(headerDir), self.fn, "header")
                self.xml_filepath = os.path.join(xml_dirpath, xml_filename)
                self.xml_cp_filename = "EMD-" + self.id_num
                xml_cp_dirpath = os.path.join(str(workDir), "out")
                self.xml_cp_filepath = os.path.join(xml_cp_dirpath, self.xml_cp_filename)
                if not os.path.exists(self.xml_cp_filepath):
                    os.mkdir(self.xml_cp_filepath)
                print(self.xml_filepath)
                self.out_cp_filepath = (self.xml_cp_filepath + "/")

                ######## Extract IDs (EMDB, PDB and UNIPROT) from header file ###########
                self.query = self.extracting_IDs(self.xml_filepath, self.pdbefile)

                ####### CREATING FASTA FILE FOR BLAST ########
                if 'PDB_ID' not in self.query[1:]:
                    #print(self.query[self.query.index('PDB_ID') - 1])
                    self.create_fasta_run_blastp(self.out_cp_filepath, self.xml_filepath)
                    out_fasta = os.path.join(self.out_cp_filepath, "*.out")
                    for self.out_fasta_file in glob.glob(out_fasta):
                    #os.chdir(self.out_cp_filepath)
                    #for self.out_fasta_file in glob.glob("*.out"):
                        print(self.out_fasta_file)
                        self.split_files(self.out_fasta_file)
                        self.purge(self.out_cp_filepath)
                    os.chdir(self.out_cp_filepath)
                    all_uni = []
                    ######## Extract UNIPROT #######
                    for filena in glob.glob("*.txt"):
                        self.filena = filena
                        self.fs = (self.out_cp_filepath + filena)
                        #print("CHE", self.outFile_uniprot, self.fs)
                        bl_uniprot = self.extract_uniprot_from_blast(self.outFile_uniprot, self.fs)
                        #print("UNI", filena, bl_uniprot)
                        all_uni.extend(bl_uniprot)
                    self.query.extend(all_uni)

                if 'PDB_ID' in self.query[1:]:
                    #print(self.query[self.query.index('PDB_ID') - 1])
                    self.uni_ID_name = self.extracting_UniprotFromPDBe(self.query, self.outFile_uniprot)
                self.query.extend(self.uni_ID_name)

                #print("QUERY+METHOD", self.query)
                ############## QUERING EXTRACTED EMDB, PDB and UNIPROT IDs ###########
                QI_CP_ID_name = self.quering_IDs(self.query, self.mapFile)
                self.CPID = QI_CP_ID_name[0]
                self.CPname = QI_CP_ID_name[1]
                self.writing_outFile(self.out_cp_filepath, self.fn, self.CPID, self.CPname)


                self.scientific_name_query = self.extracting_scientific_name(self.xml_filepath)
                #print("Number of SciName queries", len(self.scientific_name_query))
                if len(self.scientific_name_query) > 0:
                    QP_CP_ID_name = self.quering_pdbeSciFile(self.query, self.xml_filepath, self.mapFile,
                                                               self.scientific_name_query, self.pdbefile,
                                                               self.pdbeSciFile)
                    ############## IF COMPLEX PORTAL REST API NEEDS TO BE USED FOR QUERY ###############
                    #self.CP_ID_name = self.quering_sciName(self.query, self.xml_filepath,
                    #                                       self.mapFile,  self.scientific_name_query)
                    self.CPID = QP_CP_ID_name[0]
                    self.CPname = QP_CP_ID_name[1]
                    self.writing_outFile(self.out_cp_filepath, self.fn, self.CPID, self.CPname)

        self.sort_outFile(self.mapFile, self.sort_mapFile)

    def EMD_entry_ID(self, xml_filepath):
        """
        Get the EMDB entry ID number
        """
        with open(xml_filepath, 'r') as filexml:
            tree = ET.parse(filexml)
            root = tree.getroot()
            emdb_value = root.attrib
            emd_id = emdb_value.get('emdb_id')
        return emd_id

    def extracting_IDs(self, xml_filepath, pdbefile):
        """
        Extract the IDs (EMDB, PDB from both header file and PDBe files and UNIPROT)
        """
        query = []
        with open(xml_filepath, 'r') as filexml:
            tree = ET.parse(filexml)
            root = tree.getroot()
            a = root.attrib
            emd_id = a.get('emdb_id')
            query.extend([emd_id, 'EMDB_ID'])
            for x in list(root.iter('pdb_reference')):
                qs = x.find('pdb_id').text
                query.extend([qs, "PDB_ID"])
                if pdbefile:
                    with open(pdbefile, 'r') as pdbef:
                        for line in pdbef:
                            if qs not in line:
                                pass
                            else:
                                if qs in line:
                                    pdbe_com = line.split(',')[1]
                                    cid = "complex_id:" + pdbe_com
                                    query.extend([cid, "PDBe"])
            if list(root.iter('protein_or_peptide')):
               for x in list(root.iter('protein_or_peptide')):
                   qs = x.find('sequence')
                   if qs.find('external_references') is not None:
                        if qs.find('external_references').attrib['type'] == 'UNIPROTKB':
                            q = qs.find('external_references').text
                            query.extend([q, "UNIPROT"])
                   else:
                       pass
        print("Extracted_IDs + METHOD", query)
        return query

    def extracting_scientific_name(self, xml_filepath):
        """
        Extracting the scientific name from the header file
        """
        filexml = open(xml_filepath, 'r')
        tree = ET.parse(filexml)
        root = tree.getroot()
        query = []
        if list(root.iter('protein_or_peptide')):
            for x in list(root.iter('protein_or_peptide')):
                mol_id = x.attrib['macromolecule_id']
                qs = x.find('name').text
                for ns in (root.iter('natural_source')):
                   for child in root:
                       if 'organism' in child.attrib:
                           nat_sou = ns.find('organism').text
                           tid = ns.find('organism')
                           taxid = tid.attrib['ncbi']
                           query.extend([qs, nat_sou, taxid])
        print("SCINAME", query)
        return query


    def create_fasta_run_blastp(self, out_cp_filepath, xml_filepath):
        """
        Creating the fasta format file to run BLASTP
        """
        filexml = open(xml_filepath, 'r')
        tree = ET.parse(filexml)
        root = tree.getroot()
        a = root.attrib
        emd_id = a.get('emdb_id')
        fasta_file = (out_cp_filepath + emd_id + ".fasta")
        with open(fasta_file, "w") as f:
            if list(root.iter('protein_or_peptide')):
                for x in list(root.iter('protein_or_peptide')):
                    mol_id = x.attrib['macromolecule_id']
                    n = x.find('name').text
                    nat_sou = "UNKNOWN"
                    for ns in list(root.iter('natural_source')):
                        if ns.find('organism').text is not None:
                            nat_sou = ns.find('organism').text
                        else:
                            nat_sou = "unidentified"
                    qs = x.find('sequence')
                    if qs.find('string') is not None:
                        seq = qs.find('string').text
                        unk = re.findall('(UNK)+', seq)
                        if unk:
                            pass
                        else:
                            if not unk:
                                seq = seq.replace("(UNK)", "X")
                                seq = seq.replace("UNK", "X")
                                with open(fasta_file, "a") as file:
                                    file.write(">" + emd_id + "|" + mol_id + "|" + nat_sou + "|" + n + '\n' + seq + '\n')
                    else:
                        if os.path.exists(fasta_file):
                            os.remove(fasta_file)
            else:
                if fasta_file:
                    os.remove(fasta_file)
            if os.path.exists(fasta_file):
                db_path = (str(self.workDir) + "/" + "uniprot_sprot")
                qout = (out_cp_filepath + emd_id + ".out")
                command = ["blastp", "-query", "%s" % fasta_file, "-db", "%s" % db_path, "-out",
                           "%s" % qout, "-evalue",
                           "1e-40"]
                subprocess.call(command)
        return

    def files(self, out_cp_filepath, xml_cp_filename):
        """
        Create filenames with running index as prefix
        """
        n = 0
        while True:
            n += 1
            sp_fn = (out_cp_filepath + '%d_' % n + xml_cp_filename + ".txt")
            yield open(sp_fn, 'w')

    def split_files(self, out_fasta_file):
        """
        Split as separate files from fasta file for each sequence
        """
        pat = 'Query='
        fs = self.files(self.out_cp_filepath, self.xml_cp_filename)
        outfile = next(fs)
        with open(out_fasta_file) as wfile:
            for line in wfile:
                if pat not in line:
                    outfile.write(line)
                else:
                    items = line.split(pat)
                    outfile.write(items[0])
                    for item in items[1:]:
                        outfile = next(fs)
                        outfile.write(pat + item)

    def purge(self, out_cp_filepath):
        """
        Purge if file exists
        """
        pattern = "1_EMD"
        for f in os.listdir(out_cp_filepath):
            if re.search(pattern, f):
                os.remove(os.path.join(out_cp_filepath, f))

    def extract_uniprot_from_blast(self, outFile_uniprot, fs):
        """
        Extracting the UNIPROT ID from BLASTP output
        """
        ext_uniprot = []
        with open(fs) as f:
            print("File_blast", fs)
            emd_seq_tag = os.path.splitext(os.path.basename(fs))[0]
            spe_arr = ['Homo sapiens', 'Mus musculus', 'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)',
                       'Arabidopsis thaliana', 'Escherichia coli (strain K12)', 'Caenorhabditis elegans',
                       'Rattus norvegicus', 'Gallus gallus', 'Bos taurus', 'Drosophila melanogaster',
                       'Schizosaccharomyces pombe (strain 972 / ATCC 24843)', 'Canis lupus familiaris',
                       'Danio rerio', 'Xenopus laevis', 'Oryctolagus cuniculus', 'Sus scrofa', 'Lymnaea stagnalis',
                       'Pseudomonas aeruginosa (strain ATCC 15692)', 'Tetronarce californica', 'Torpedo marmorata']
            resl = [next(f) for x in range(6)]
            strings = re.findall(r'Length=\d*', str(resl), re.S)
            seq_len = strings[0].split("=")[-1]
            sep_n = re.findall('EMD-(.*)Length=', str(resl), re.S)
            spe_name = sep_n[0].split("|")[-2]
            for result in re.findall('>(.*?)Query', f.read(), re.S):
                ident = re.findall(r'Identities = \d*/\d*', result) #Identities = \d*/\d* \(\d+%\)
                ident_len = int(ident[0].split("/")[-1])
                iden = re.findall(r'Identities = \d*/\d* \(\d*', result, re.S)
                ident_per = int(iden[0].split("(")[-1])
                if int(seq_len) == ident_len:
                #if percentage_match <= ident_len:    #### FOR PERCENTAGE MATCHING OF INPUT SEQ
                    if ident_per >= 100:
                        unip = re.findall('sp(.*)OS=', result, re.S)
                        uniprot = unip[0].split("|")[-2]
                        uniprot_n = unip[0].split("|")[-1]
                        uniprot_name = uniprot_n.replace('\n', '')
                        ev = re.findall('Expect =(.*)Method', result, re.S)
                        evalue = ("Evalue=" + ev[0].split(",")[0])
                        sp_n = re.findall('OS=(.*)=', result, re.S)
                        bsp_name = sp_n[0].split("OX")[0]
                        bsp_name = bsp_name.replace('\n', '')
                        print("SPECIES_NAME", spe_name, "VS",  bsp_name)
                        if any(x in bsp_name for x in spe_arr):
                              #unip_out = (emd_seq_tag + "," + uniprot + "," + uniprot_name + "," +
                                          #bsp_name + "," + evalue + '\n')
                              #uni_file += uniprot
                              ext_uniprot.extend([uniprot, 'BLAST'])
                        with open(outFile_uniprot, "a") as fo:
                            EMDBID = emd_seq_tag.split('_')[1]
                            fo.write(EMDBID + "\t" + uniprot + "\n")
        print("Extracted_uniprot", ext_uniprot)
        return ext_uniprot

    def grep_identities(self, fs):
        """
        Grep the Identities from the header file
        """
        file_object = open(self.fs, 'r')
        strings = re.findall(r'Identities = \d*/\d* \(\d+%\)', file_object.read())
        print(strings)

    def grep_seq_len_name(self, fs):
        """
        Grep the sequence length
        """
        global whole, sep_name
        file_object = open(self.fs, 'r')
        for result in re.findall('Query=(.*?)Sequences', file_object.read(), re.S):
            strings = re.findall(r'Length=\d*', result, re.S)
            whole = strings[0].split("=")[-1]
            sep_n = re.findall('EMD-(.*)Length=', result, re.S)
            sep_name = sep_n[0].split("|")[-2]
        return whole, sep_name

    def quering_IDs(self, query, mapFile):
        """
        Querying the complex portal database for extracted IDs
        """
        if len(self.query) <= 0:
            pass
        else:
            CP_ID = []
            CP_name = []
            CP_n = ""
            queryID = self.query[::2]
            meth = self.query[1::2]
            for y in range(len(queryID)):
                EMDBID = queryID[0]
                queryString = queryID[y]
                method = meth[y]
                url = CP_baseurl + (queryString)
                def save_sslcontext(obj):
                    return obj.__class__, (obj.protocol,)
                copyreg.pickle(ssl.SSLContext, save_sslcontext)
                context = ssl.create_default_context()
                foo = pickle.dumps(context)
                gcontext = pickle.loads(foo)
                cpjson = urlopen(url, context=gcontext).read()
                cpjdata = json.loads(cpjson.decode('utf-8'))
                size = cpjdata['size']
                for ind in range(size):
                    CP_ind = cpjdata['elements'][ind]['complexAC']
                    #if CP_ind not in [str(y) for x in CP_ID for y in x.split()]:
                    CP_ID.append(CP_ind)
                    CP_n = cpjdata['elements'][ind]['complexName']
                    CP_name.append(CP_n)
                    with open(self.mapFile, 'a') as f:
                      f.write(EMDBID + "\t" + queryString + "\t" + CP_ind + "\t" + CP_n + "\t" + method + '\n')
            print("OUTPUT_TO_FILE", CP_ID, CP_name)
        return CP_ID, CP_name

    def extracting_UniprotFromPDBe(self, query, outFile_uniprot):
        """
        Extracting the UNIPROT ID from PDBe API if model exists for the entry
        """
        if len(self.query) <= 0:
            pass
        else:
            if len(self.query) > 0:
                keys = []
                uni_ID = []
                queryID = []
                EMDBID = (self.query[self.query.index('EMDB_ID') - 1])
                if 'PDB_ID' in self.query[1:]:
                    queryID.append(self.query[self.query.index('PDB_ID') - 1])
                    for y in range(len(queryID)):
                        queryString = queryID[y]
                        url = pdb_baseurl + (queryString)
                        def save_sslcontext(obj):
                            return obj.__class__, (obj.protocol,)
                        copyreg.pickle(ssl.SSLContext, save_sslcontext)
                        context = ssl.create_default_context()
                        foo = pickle.dumps(context)
                        gcontext = pickle.loads(foo)
                        try:
                            pdbjson = urlopen(url, context=gcontext).read()
                            pdbjdata = json.loads(pdbjson.decode('utf-8'))
                            jdata = pdbjdata[queryString]['UniProt']
                            for key in jdata.keys():
                                uni_ID.extend([key, "UNIPROT_from_PDBe"])
                                keys.append(key)
                            with open(outFile_uniprot, "a") as fo:
                                fo.write(EMDBID + '\t' + str(keys) + '\n')
                        except HTTPError:
                            pass
                    print("UNI_ID", uni_ID)
        return uni_ID

    def quering_sciName(self, query, xml_filepath, mapFile, scientific_name_query):
        """
        Querying the complex portal for the extracted scientific name
        """
        cpj_quering_ID = self.quering_IDs(self.query, self.mapFile)
        print("Query_sciName", cpj_quering_ID)
        CP_ID = cpj_quering_ID[0]
        CP_name = cpj_quering_ID[1]
        CP_n = ""
        print(CP_ID, CP_name)
        if len(self.scientific_name_query) <= 0:
            pass
        else:
            mol_id = len(self.scientific_name_query)
            ind = int(mol_id/3)
            for y in range(ind):
                sci_ind = int((3*y))
                spe_ind = int((3*y)+1)
                tax_id = int((3*y)+2)
                sci_qs = self.scientific_name_query[sci_ind]
                sci_queryString = sci_qs.replace(' ', '%20')
                sci_queryString = sci_queryString.replace("'", "%27")
                sci_queryString = sci_queryString.replace("(", "%28")
                sci_queryString = sci_queryString.replace(")", "%29")
                sci_queryString = sci_queryString.replace(",", "%2C")
                sci_queryString = sci_queryString.replace("-", "%2D")
                sci_queryString = sci_queryString.replace("/", "%2F")
                spe_qs = self.scientific_name_query[spe_ind]
                spe_queryString = spe_qs.replace(' ', '%20')
                spe_queryString = spe_queryString.replace('(', '%28')
                spe_queryString = spe_queryString.replace(')', '%29')
                spe_queryString = spe_queryString.replace("'", "%27")
                spe_queryString = spe_queryString.replace(',', '%2C')
                spe_queryString = spe_queryString.replace('-', '%2D')
                spe_queryString = spe_queryString.replace('/', '%2F')
                url = CP_baseurl + 'complex_alias:' + "%22" + (sci_queryString) + "%22" + \
                      "?first=0&number=30&filters=species_f:%28%22" + (spe_queryString) +"%22%29&facets=species_f"
                def save_sslcontext(obj):
                    return obj.__class__, (obj.protocol,)
                copyreg.pickle(ssl.SSLContext, save_sslcontext)
                context = ssl.create_default_context()
                foo = pickle.dumps(context)
                gcontext = pickle.loads(foo)
                cpjson = urlopen(url, context=gcontext).read()
                cpjdata = json.loads(cpjson.decode('utf-8'))
                print(cpjdata)
                size = cpjdata['size']
                for ind in range(size):
                    CP_ind = cpjdata['elements'][ind]['complexAC']
                    if CP_ind not in [str(y) for x in CP_ID for y in x.split()]:
                       CP_ID.append(CP_ind)
                       CP_n = cpjdata['elements'][ind]['complexName']
                       CP_name.append(CP_n)
                    EMDBID = self.EMD_entry_ID(self.xml_filepath)
                    with open(self.mapFile, 'a') as f:
                      f.write(EMDBID + "\t" + sci_queryString + "\t" + CP_ind + "\t" + CP_n + "\t"
                              + "Scientific name" + '\n')
        print("FINAL", CP_ID, CP_name)
        return CP_ID, CP_name

    def quering_pdbeSciFile(self, query,xml_filepath, mapFile, scientific_name_query, pdbefile, pdbeSciFile):
        """
        Querying the files created by the PDBe to check for scientific name of the EMDB entry
        """
        cpj_quering_ID = self.quering_IDs(self.query, self.mapFile)
        print("QUERY_PDBe", cpj_quering_ID)
        CP_ID = cpj_quering_ID[0]
        CP_name = cpj_quering_ID[1]
        CP_n = ""
        if len(self.scientific_name_query) <= 0:
            pass
        else:
            mol_id = len(self.scientific_name_query)
            ind = int(mol_id/3)
            for y in range(ind):
                sci_ind = int((3*y))
                spe_ind = int((3*y)+1)
                tax_id = int((3*y)+2)
                sci_queryString = self.scientific_name_query[sci_ind]
                spe_qs = self.scientific_name_query[spe_ind]
                taxid = self.scientific_name_query[tax_id]
                print("QUERY_SCINAME:" + sci_queryString + "," + spe_qs + "," + taxid)
                pdbesf = open(self.pdbeSciFile, 'r')
                for ln in pdbesf:
                    if sci_queryString not in ln:
                        pass
                    else:
                        t = ln.split(',')[-2]
                        pcpx = ln.split(',')[-1]
                        if len(t) <= 1:
                            pass
                        else:
                            if int(t) == int(taxid):
                                f = open(pdbefile, 'r')
                                for x in f:
                                    if pcpx not in x:
                                        pass
                                    else:
                                        pdbcpx = x.split(',')[1]
                                        cid = "complex_id:" + pdbcpx
                                        url = CP_baseurl + (cid)
                                        def save_sslcontext(obj):
                                            return obj.__class__, (obj.protocol,)
                                        copyreg.pickle(ssl.SSLContext, save_sslcontext)
                                        context = ssl.create_default_context()
                                        foo = pickle.dumps(context)
                                        gcontext = pickle.loads(foo)
                                        cpjson = urlopen(url, context=gcontext).read()
                                        cpjdata = json.loads(cpjson.decode('utf-8'))
                                        print(cpjdata)
                                        size = cpjdata['size']
                                        for ind in range(size):
                                            CP_ind = cpjdata['elements'][ind]['complexAC']
                                            if CP_ind not in [str(y) for x in CP_ID for y in x.split()]:
                                               CP_ID.append(CP_ind)
                                               CP_n = cpjdata['elements'][ind]['complexName']
                                               CP_name.append(CP_n)
                                            EMDBID = self.EMD_entry_ID(self.xml_filepath)
                                            with open(self.mapFile, 'a') as f:
                                              f.write(EMDBID + "\t" + sci_queryString + "\t" + CP_ind + "\t" +
                                                      CP_n + "\t" + "PDBe" + '\n')
        print("FINAL", CP_ID, CP_name)
        return CP_ID, CP_name

    def sort_outFile(self, mapFile, sort_mapFile):
        """
        Sort the EMDB mapping to complex portal output file with EMDB_ID numbers
        """
        data = csv.reader(open(self.mapFile), delimiter='\t')
        sortedlist = sorted(data, key=operator.itemgetter(0))
        with open(self.sort_mapFile, "a") as f:
            f.write("EMDB_ID" + "\t" + "queryString" + "\t" + "Complex_ID" + "\t" +
                    "Complex_Name" + "\t" + "MethodOfQuery" + '\n')
            fileWriter = csv.writer(f, delimiter='\t')
            for row in sortedlist:
                fileWriter.writerow(row)
        return

    def addNextCP_ID(self, xml_filepath, CP_ID, CP_name, size):
        """
        Adding complex portal ID and name after mapping to the EMDB header file
        """
        filexml = open(xml_filepath, 'r')
        parser = ET.XMLParser(remove_blank_text=True)
        tree = ET.parse(filexml, parser)
        root = tree.getroot()
        a = root.find('sample')
        if size == 0:
            pass
        else:
            b = ET.SubElement(a, "complex_reference")
            for ind in range(size):
                cp = ET.SubElement(b, "complex_ID")
                cp.text = CP_ID[ind]
                cpn = ET.SubElement(b, "complex_name")
                cpn.text = CP_name[ind]
        return tree

    def writing_outFile(self, out_cp_filepath, id_num, CPID, CPname):
        """
        Writing the EMDB mapping to complex portal entries
        """
        out_file = (self.out_cp_filepath + "CP_EMD-" + self.id_num)
        print(out_file)
        if os.path.exists(out_file):
            os.remove(out_file)
        if len(CPID) <= 0:
            cpString = etree.tostring(self.addNextCP_ID(self.xml_filepath, "None", "None", len(CPID)), pretty_print=True)
            with open(out_file, "ab") as file:
              file.write(cpString)
        else:
            cpString = etree.tostring(self.addNextCP_ID(self.xml_filepath, CPID, CPname, len(CPID)), pretty_print=True)
            with open(out_file, "ab") as file:
              file.write(cpString)

if __name__== "__main__":
    ######### Command : python /Users/amudha/project/ComplexPortal/ComplexPortalMapping.py -w /Users/amudha/project/
    # -f /Users/amudha/project/EMD_XML/ -p /Users/amudha/project/pdbeFiles/ ############################

    prog = "ComplexPortalMapping"
    usage = """
            Mapping EMDB entries to Complex portal.

            Example:
            python ComplexPortalMapping.py -w '[{"/path/to/working/folder"}]'
            -f '[{"/path/to/EMDB/header/files/folder"}]'
            -p '[{"/path/to/PDBe/files/folder"}]'
          """
    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path .")
    parser.add_argument('-f', '--headerDir', type=Path, help="Directory path to the EMDB version 3.0 header files.")
    parser.add_argument('-p', '--PDBeDir', type=Path, help="Directory path to the PDBe files.")
    args = parser.parse_args()
    CPMapping(args.workDir, args.headerDir, args.PDBeDir)