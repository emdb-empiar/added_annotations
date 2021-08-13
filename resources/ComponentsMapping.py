import os
import logging
from gemmi import cif
from multiprocessing import Pool

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(funcName)s:%(message)s')
file_handler = logging.FileHandler('logging_components.log')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

### TO DO LIST
#### Replace (logger.debug(HET, "NOT IN PDB_CCD") with corresponding resource API,
# as of now no entry has HET which is not in CCD #

class ComponentsMapping:
    """
    Extracting PDB_IDs from header and get the HET code from PDBe components.cif file.
    Querying with the extracted HET code for mapping the EMDB entries to the Chemical Component Dictionary (CCD) for
    mapping to various database like ChEMBL, ChEBI and DrugBank.
    """

    def __init__(self, workDir, ligands, components_cif):
        self.workDir = workDir
        self.ligands = ligands
        self.chembl_map = {}
        self.chebi_map = {}
        self.drugbank_map = {}
        self.components_cif = components_cif

    def execute(self, threads):
        ####### Extract only the HET_code, resource name and IDs from the PDBe componenets.cif file #####
        self.chembl_map, self.chebi_map, self.drugbank_map = self.extract_resources_from_cif()

        ###### Mapping HET_CODE TO CHEMBL, CHEBI and DRUGBANK ########
        with Pool(processes=threads) as pool:
            self.ligands = pool.map(self.worker, self.ligands)
        return self.ligands

    def worker(self, ligand):
        HET = ligand.HET
        if ligand.provenance == "AUTHOR":
            return ligand
        else:
            if HET in self.chembl_map:
                ligand.chembl_id = self.chembl_map[HET]
                ligand.provenance = "CCD"
            if HET in self.chebi_map:
                ligand.chebi_id = self.chebi_map[HET]
                ligand.provenance = "CCD"
            if HET in self.drugbank_map:
                ligand.drugbank_id = self.drugbank_map[HET]
                ligand.provenance = "CCD"
            else:
                logger.debug("NOT IN CCD %s" % (HET))  #### Replace with corresponding resource API
        return ligand

    def extract_resources_from_cif(self):
        """
        Extract only the external mapping for the HET_CODE from the pdbe components.cif file
        """
        chembl_map = {}
        chebi_map = {}
        drugbank_map = {}

        doc = cif.read_file(self.components_cif)
        for block in doc:
            HET = block.name
            name = block.find_value('_chem_comp.name')
            for db,value in block.find('_pdbe_chem_comp_external_mappings.', ['resource', 'resource_id']):
                if db == "ChEMBL":
                    chembl_map[HET] = value
                elif db == "ChEBI":
                    chebi_map[HET] = value
                elif db == "DrugBank":
                    drugbank_map[HET] = value
        return chembl_map, chebi_map, drugbank_map

    def write_ligands(self):
        chembl_file = os.path.join(self.workDir, "emdb_chembl.tsv")
        chebi_file = os.path.join(self.workDir, "emdb_chebi.tsv")
        db_file = os.path.join(self.workDir, "emdb_drugbank.tsv")
        with open(chembl_file, 'w') as f1, open(chebi_file, 'w') as f2, open(db_file, 'w') as f3:
            f1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "SAMPLE_ID", "HET_CODE", "COMP_NAME", "COMP_COPIES", "ChEMBL_ID", "PROVENANCE"))
            f2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "SAMPLE_ID", "HET_CODE", "COMP_NAME", "COMP_COPIES", "ChEBI_ID", "PROVENANCE"))
            f3.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("EMDB_ID", "SAMPLE_ID", "HET_CODE", "COMP_NAME", "COMP_COPIES", "DrugBank_ID", "PROVENANCE"))

            for ligand in self.ligands:
                f1.write(ligand.get_chembl_tsv())
                f2.write(ligand.get_chebi_tsv())
                f3.write(ligand.get_drugbank_tsv())
