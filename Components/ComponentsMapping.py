import os, sys, csv
import lxml.etree as ET
from glob import glob
import logging
from gemmi import cif
from ComplexPortal.ComplexPortalMapping import CPMapping

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(funcName)s:%(message)s')
file_handler = logging.FileHandler('logging_components.log')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

#chEMBL_ftp = r'/nfs/ftp/pub/databases/chembl/ChEMBLdb/latest/'
chEMBL_ftp = r'/Users/amudha/project/chEMBLdb/latest/'
cif_filepath = r'/Users/amudha/project/'
#cif_filepath = r'/nfs/ftp/pub/databases/msd/pdbechem_v2/'

### TO DO LIST
#### Replace (logger.debug(HET, "NOT IN PDB_CCD") with corresponding resource API, as of now no entry has HET which is not in CCD #

class ComponentsMap:
    """
    Extracting PDB_IDs from header and get the HET code from PDBe components.cif file.
    Querying with the extracted HET code for mapping the EMDB entries to the Chemical Component Dictionary (CCD) for
    mapping to various database like ChEMBL, ChEBI and DrugBank.
    """

    def __init__(self, workDir, headerDir):
        self.annotations = []
        self.workDir = workDir
        self.headerDir = headerDir
        self.CCD_HET = set()
        self.HET_map = set()
        self.chembl_map = set()
        self.chebi_map = set()
        self.drugbank_map = set()
        self.componentsDir = os.path.join(self.workDir, "git_code/added_annotations/Components")

    def execute_annotations(self):
        ####### Extract only the HET_code, resource name and IDs from the PDBe componenets.cif file #####
        self.extract_resources_from_cif(cif_filepath, self.componentsDir)

        ###### Fetch header files for query ########
        for fn in glob(os.path.join(str(self.headerDir), '*')):
            print(fn)
            id_num = fn.split('-')[1]
            xml_filename = "emd-" + id_num + "-v30.xml"
            xml_dirpath = os.path.join(str(self.headerDir), fn, "header")
            xml_filepath = os.path.join(xml_dirpath, xml_filename)

            ####### List of all HET codes from PDB CCD ##########
            self.get_HET_codes(cif_filepath)

            ####### Extract the pdb_ids and store author provided annotations #######
            HETs, pdb_ids = self.extracting_IDs(xml_filepath)

            ######## If HET_code in CCD exists, then map the ChEMBL, ChEBI and DrugBank to EMDB and write to file #####
            for HET in HETs:
                if HET in self.CCD_HET:
                    self.external_mapping_from_cif(self.componentsDir, "EMD-" + id_num, HET)
                    self.write_chembl_map()
                    self.sort_emdb_chembl_map()
                    self.write_chebi_map()
                    self.sort_emdb_chebi_map()
                    self.write_drugbank_map()
                    self.sort_emdb_drugbank_map()
                if not HET in self.CCD_HET:
                    logger.debug(HET, "NOT IN PDB_CCD")  #### Replace with corresponding resource API

    def extracting_IDs(self, xml_filepath):
        """
        Extract the IDs (EMDB, PDB, HET_code from both EMDB header file. If model exists, HET_CODE from PDBE cif file)
        """
        pdb_ids = set()
        HET = set()

        with open(xml_filepath, 'r') as filexml:
            tree = ET.parse(filexml)
            root = tree.getroot()
            a = root.attrib
            emd_id = a.get('emdb_id')
            for x in list(root.iter('pdb_reference')):
                model = x.find('pdb_id').text.lower()
                pdb_ids.add(model)
                try:
                    cifFile = os.path.join(str(cif_filepath), "mmCIF", model + "_updated.cif")
                    doc = cif.read_file(cifFile)  # copy all the data from mmCIF file
                    block = doc.sole_block()  # mmCIF has exactly one block
                    for element in block.find_loop("_pdbx_entity_nonpoly.comp_id"):
                        if element not in HET:
                            logger.debug(element)
                            HET.add(element)
                            self.HET_map.add((emd_id, element, "CCD"))
                except Exception as e:
                    logger.warning(e)
            if len(list(root.iter('pdb_reference'))) == 0:
                if list(root.iter('ligand')):
                    for x in list(root.iter('ligand')):
                        if x is not None:
                            compound = x.find('formula').text
                            HET.add(compound)
                            self.HET_map.add((emd_id, compound, "AUTHOR"))
        return HET, pdb_ids

    def get_HET_codes(self, cif_filepath):
        """
        Extract only the HET_CODEs from the pdbe components.cif file
        """
        filecif = os.path.join(str(cif_filepath), "components.cif")
        doc = cif.read_file(filecif)
        for x in range(len(doc)):
            block = doc[x]
            HET = block.name
            self.CCD_HET.add(HET)

    def extract_resources_from_cif(self, cif_filepath, componentsDir):
        """
        Extract only the external mapping for the HET_CODE from the pdbe components.cif file and write to new file
        """
        with open(os.path.join(str(componentsDir), "components_resources_ID.tsv"), 'w') as fileID:
            fileID.write("%s\t%s\t%s\n" % ("HET_CODE", "SOURCE", "ID"))
            filecif = os.path.join(str(cif_filepath), "components.cif")
            doc = cif.read_file(filecif)
            for x in range(len(doc)):
                block = doc[x]
                HET_name = block.find_value('_chem_comp.name')
                for element in block.find('_pdbe_chem_comp_external_mappings.', ['comp_id', 'resource', 'resource_id']):
                    fileID.write("%s\t%s\t%s\t%s\n" % (element[0], HET_name, element[1], element[2]))

    def external_mapping_from_cif(self, componentsDir, emdb_id, HET):
        """
        Annotating the extracted HET_CODE to various database
        """
        with open(os.path.join(str(componentsDir), "components_resources_ID.tsv"), 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader, None)
            for row in reader:
                formula = row[0]
                if (HET == formula and row[2] == "ChEMBL"):
                    self.chembl_map.add((emdb_id, row[0], row[1], row[2], "CCD"))
                    logger.debug((emdb_id, row[0], row[1], row[2], "CCD"))
                if (HET == formula and row[2] == "ChEBI"):
                    self.chebi_map.add((emdb_id, row[0], row[1], row[2], "CCD"))
                    logger.debug((emdb_id, row[0], row[1], row[2], "CCD"))
                if (HET == formula and row[2] == "DrugBank"):
                    self.drugbank_map.add((emdb_id, row[0], row[1], row[2], "CCD"))
                    logger.debug((emdb_id, row[0], row[1], row[2], "CCD"))
                    print(emdb_id, row[0], row[1], row[2], "CCD")

    def write_chembl_map(self):
        """
        Write ChEMBL annotations to a file
        """
        filepath = os.path.join(self.componentsDir, "emdb_chembl.tsv")
        with open(filepath, 'w') as f:
            for emdb_id, HETs, HET_name, chembl_id, method in self.chembl_map:
                f.write("%s\t%s\t%s\t%s\t%s\n" % (emdb_id, HETs, HET_name, chembl_id, method))

    def sort_emdb_chembl_map(self):
        """
        Sort the ChEMBL annotations with respect to EMDB_ID to a file
        """
        with open(os.path.join(self.componentsDir, "emdb_chEMBL.tsv"), 'r') as lines:
            with open(os.path.join(self.componentsDir, "sorted_emdb_chEMBL.tsv"), 'w') as sort_file:
                sort_file.write("%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "HET_CODE", "COMP_NAME", "ChEMBL_ID", "QUERY_METHOD"))
                for line in sorted(lines, key=lambda line: line.split()[0]):
                    sort_file.write(line)

    def write_chebi_map(self):
        """
        Write ChEBI annotations to a file
        """
        filepath = os.path.join(self.componentsDir, "emdb_chebi.tsv")
        with open(filepath, 'w') as f:
            for emdb_id, HETs, HET_name, chebi_id, method in self.chebi_map:
                f.write("%s\t%s\t%s\t%s\t%s\n" % (emdb_id, HETs, HET_name, chebi_id, method))

    def sort_emdb_chebi_map(self):
        """
        Sort ChEBI annotations with respect to EMDB_ID to a file
        """
        with open(os.path.join(self.componentsDir, "emdb_chebi.tsv"), 'r') as lines:
            with open(os.path.join(self.componentsDir, "sorted_emdb_chebi.tsv"), 'w') as sort_file:
                sort_file.write("%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "HET_CODE", "COMP_NAME", "ChEBI_ID", "QUERY_METHOD"))
                for line in sorted(lines, key=lambda line: line.split()[0]):
                    sort_file.write(line)

    def write_drugbank_map(self):
        """
        Write DrugBank annotations to a file
        """
        filepath = os.path.join(self.componentsDir, "emdb_drugbank.tsv")
        with open(filepath, 'w') as f:
            for emdb_id, HETs, HET_name, drugbank_id, method in self.drugbank_map:
                f.write("%s\t%s\t%s\t%s\t%s\n" % (emdb_id, HETs, HET_name, drugbank_id, method))

    def sort_emdb_drugbank_map(self):
        """
        Sort DrugBank annotations with respect to EMDB_ID to a file
        """
        with open(os.path.join(self.componentsDir, "emdb_drugbank.tsv"), 'r') as lines:
            with open(os.path.join(self.componentsDir, "sorted_emdb_drugbank.tsv"), 'w') as sort_file:
                sort_file.write("%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "HET_CODE", "COMP_NAME", "DrugBank_ID", "QUERY_METHOD"))
                for line in sorted(lines, key=lambda line: line.split()[0]):
                    sort_file.write(line)