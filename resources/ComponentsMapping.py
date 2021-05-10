import os
import logging
from gemmi import cif

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
#### Replace (logger.debug(HET, "NOT IN PDB_CCD") with corresponding resource API,
# as of now no entry has HET which is not in CCD #

class ComponentsMap:
    """
    Extracting PDB_IDs from header and get the HET code from PDBe components.cif file.
    Querying with the extracted HET code for mapping the EMDB entries to the Chemical Component Dictionary (CCD) for
    mapping to various database like ChEMBL, ChEBI and DrugBank.
    """

    def __init__(self, workDir, ligands):
        self.annotations = []
        self.workDir = workDir
        self.CCD_HET = set()
        self.HET_map = set()
        self.HET_info = set()
        self.chembl_map = set()
        self.chebi_map = set()
        self.drugbank_map = set()
        self.ligands = ligands
        self.componentsDir = os.path.join(self.workDir, "git_code/added_annotations/Components")

    def execute_annotations(self):
        ####### Extract only the HET_code, resource name and IDs from the PDBe componenets.cif file #####
        self.extract_resources_from_cif(cif_filepath)

        ####### List of all HET codes from PDB CCD ##########
        self.get_HET_codes(cif_filepath)

        ###### Mapping HET_CODE TO CHEMBL, CHEBI and DRUGBANK ########
        for ligand in self.ligands:
            HET = ligand.HET
            if HET in self.CCD_HET:
                self.external_mapping_from_cif(ligand.emdb_id, ligand.sample_id, HET)
            if not HET in self.CCD_HET:
                logger.debug(HET, "NOT IN PDB_CCD")  #### Replace with corresponding resource API

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

    def extract_resources_from_cif(self, cif_filepath):
        """
        Extract only the external mapping for the HET_CODE from the pdbe components.cif file
        """
        filecif = os.path.join(str(cif_filepath), "components.cif")
        doc = cif.read_file(filecif)
        for x in range(len(doc)):
            block = doc[x]
            HET_name = block.find_value('_chem_comp.name')
            for element in block.find('_pdbe_chem_comp_external_mappings.', ['comp_id', 'resource', 'resource_id']):
                self.HET_info.add((element[0], HET_name, element[1], element[2]))

    def external_mapping_from_cif(self, emdb_id, lig_id, HET):
        """
        Annotating the extracted HET_CODE to various database
        """
        for row in self.HET_info:
            formula = row[0]
            if (HET == formula and row[2] == "ChEMBL"):
                self.chembl_map.add((emdb_id, lig_id, row[0], row[1], row[3], "CCD"))
                logger.debug((emdb_id, lig_id, row[0], row[1], row[3], "CCD"))
            if (HET == formula and row[2] == "ChEBI"):
                self.chebi_map.add((emdb_id, lig_id, row[0], row[1], row[3], "CCD"))
                logger.debug((emdb_id, lig_id, row[0], row[1], row[3], "CCD"))
            if (HET == formula and row[2] == "DrugBank"):
                self.drugbank_map.add((emdb_id, lig_id, row[0], row[1], row[3], "CCD"))
                logger.debug((emdb_id, lig_id, row[0], row[1], row[3], "CCD"))

    def write_chembl_map(self):
        """
        Write ChEMBL annotations to a file (sorted by EMDB_ID)
        """
        filepath = os.path.join(self.componentsDir, "emdb_chembl.tsv")
        with open(filepath, 'w') as f:
            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "SAMPLE_ID", "HET_CODE", "COMP_NAME", "ChEMBL_ID", "PROVENANCE"))
            for emdb_id, lig_id, HETs, HET_name, chembl_id, method in sorted(self.chembl_map, key=lambda x: x[0]):
                f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (emdb_id, lig_id, HETs, HET_name, chembl_id, method))

    def write_chebi_map(self):
        """
        Write ChEBI annotations to a file (sorted by EMDB_ID)
        """
        filepath = os.path.join(self.componentsDir, "emdb_chebi.tsv")
        with open(filepath, 'w') as f:
            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "SAMPLE_ID", "HET_CODE", "COMP_NAME", "ChEBI_ID", "PROVENANCE"))
            for emdb_id, lig_id, HETs, HET_name, chebi_id, method in sorted(self.chebi_map, key=lambda x: x[0]):
                f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (emdb_id, lig_id, HETs, HET_name, chebi_id, method))

    def write_drugbank_map(self):
        """
        Write DrugBank annotations to a file (sorted by EMDB_ID)
        """
        filepath = os.path.join(self.componentsDir, "emdb_drugbank.tsv")
        with open(filepath, 'w') as f:
            f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "SAMPLE_ID", "HET_CODE", "COMP_NAME", "DrugBank_ID", "PROVENANCE"))
            for emdb_id, lig_id, HETs, HET_name, drugbank_id, method in sorted(self.drugbank_map, key=lambda x: x[0]):
                f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (emdb_id, lig_id, HETs, HET_name, drugbank_id, method))
