from gemmi import cif

### TO DO LIST: ##
#### Replace (logger.debug(HET, "NOT IN PDB_CCD") with corresponding resource API,
##### as of now no entry has HET which is not in CCD ####

class ComponentsMapping:
    """
    Extracting PDB_IDs from header and get the HET code from PDBe components.cif file.
    Querying with the extracted HET code for mapping the EMDB entries to the Chemical Component Dictionary (CCD) for
    mapping to various database like ChEMBL, ChEBI and DrugBank.
    """

    def __init__(self, ligands, components_cif):
        self.ligands = ligands
        self.chembl_map = {}
        self.chebi_map = {}
        self.drugbank_map = {}
        self.components_cif = components_cif

    def execute(self):
        ####### Extract only the HET_code, resource name and IDs from the PDBe componenets.cif file #####
        self.chembl_map, self.chebi_map, self.drugbank_map = self.extract_resources_from_cif()

        ###### Mapping HET_CODE TO CHEMBL, CHEBI and DRUGBANK ########
        for ligand in self.ligands:
            ligand = self.worker(ligand)
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
                print("NOT IN CCD %s" % (HET))  #### Replace with corresponding resource API
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

    def export_tsv(self, chembl_logger, chebi_logger, drugbank_logger):
        for ligand in self.ligands:
            chembl_row = f"{ligand.emdb_id}\t{ligand.sample_id}\t{ligand.HET}\t{ligand.lig_name}\t{ligand.lig_copies}\t{ligand.chembl_id}\t{ligand.provenance}"
            chembl_logger.info(chembl_row)
            chebi_row = f"{ligand.emdb_id}\t{ligand.sample_id}\t{ligand.HET}\t{ligand.lig_name}\t{ligand.lig_copies}\t{ligand.chebi_id}\t{ligand.provenance}"
            chebi_logger.info(chebi_row)
            db_row = f"{ligand.emdb_id}\t{ligand.sample_id}\t{ligand.HET}\t{ligand.lig_name}\t{ligand.lig_copies}\t{ligand.drugbank_id}\t{ligand.provenance}"
            drugbank_logger.info(db_row)

