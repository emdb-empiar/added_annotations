import os
import logging
from gemmi import cif
from multiprocessing import Pool
from EICSS import EICSS

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(funcName)s:%(message)s')
file_handler = logging.FileHandler('logging_components.log')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

components_cif = r'/Users/amudha/project/ftp_data/pdbe/components.cif'
#components_cif = r'/nfs/ftp/pub/databases/msd/pdbechem_v2/components.cif'
#components_cif = "/Users/neli/Downloads/components.cif"

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
        self.workDir = workDir
        self.ligands = ligands
        self.chembl_map = {}
        self.chebi_map = {}
        self.drugbank_map = {}

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

        doc = cif.read_file(components_cif)
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

    # def writeXML_ligands(self):
    #     component_DB = set()
    #     headerXML = EICSS.eicss()
    #     DBs_list = EICSS.DBs_listType()
    #     list_macro_molecules = EICSS.list_macro_moleculesType()
    #     sample_annotation = EICSS.sample_annotationType()
    #
    #     em_id = self.ligands[0].emdb_id
    #     # print("FIR", em_id)
    #     for ligand in self.ligands:
    #         headerXML.set_EMDB_ID(ligand.emdb_id)
    #         # print("AC", ligand.emdb_id, ligand.HET)
    #         macro_molecule_annotation = EICSS.macro_molecule_annotationType()
    #         if ligand.emdb_id == em_id:
    #             # print("K", em_id)
    #             macro_molecule_annotation.set_macro_kind("%s" % "ligand")
    #             macro_molecule_annotation.set_macro_ID(int(ligand.sample_id))
    #             macro_molecule_annotation.set_macro_copies(int(ligand.lig_copies))
    #             macro_molecule_annotation.set_macro_name("%s" % ligand.lig_name)
    #             list_macro_molecules.add_macro_molecule_annotation(macro_molecule_annotation)
    #         if ligand.chembl_id:
    #             if "CHEMBL" not in component_DB:
    #                 DB = EICSS.DBType()
    #                 DB.set_DB_source("%s" % "CHEMBL")
    #                 DB.set_DB_version("%s" % "4.2.0")
    #                 DBs_list.add_DB(DB)
    #
    #             # crossRefDB = EICSS.crossRefDBType()
    #             # crossRefDB.set_DB_source("%s" % "ChEMBL")
    #             # crossRefDB.set_provenance("%s" % ligand.provenance)
    #             # crossRefDB.set_DB_accession_ID("%s" % ligand.chembl_id)
    #             # list_crossRefDBs = EICSS.list_crossRefDBsType()
    #             # list_crossRefDBs.add_crossRefDB(crossRefDB)
    #             # macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
    #
    #         component_DB.add("CHEMBL")
    #         if ligand.chebi_id:
    #             if "CHEBI" not in component_DB:
    #                 DB = EICSS.DBType()
    #                 DB.set_DB_source("%s" % "CHEBI")
    #                 DB.set_DB_version("%s" % "15.21")
    #                 DBs_list.add_DB(DB)
    #         component_DB.add("CHEBI")
    #         if ligand.drugbank_id:
    #             if "DRUGBANK" not in component_DB:
    #                 DB = EICSS.DBType()
    #                 DB.set_DB_source("%s" % "DRUGBANK")
    #                 DB.set_DB_version("%s" % "2021.03.30")
    #                 DBs_list.add_DB(DB)
    #         component_DB.add("DRUGBANK")
    #
    #         # print("BF", em_id, ligand.emdb_id)
    #         em_id = ligand.emdb_id
    #         # print("FIN", em_id, ligand.emdb_id)
    #         headerXML.set_DBs_list(DBs_list)
    #         headerXML.set_sample_annotation(list_macro_molecules)
    #
    #         xmlFile = os.path.join(self.workDir, ligand.emdb_id + "_eicss.xml")
    #         with open(xmlFile, 'w') as f:
    #             headerXML.export(f, 0, name_='eicss')
