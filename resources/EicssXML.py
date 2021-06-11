import os
from EICSS import EICSS

class EicssXML:
    "Writing annotations to output xml file according to the EMDB_EICSS.xsd schema "

    def __init__(self, workDir, proteins, lig_map, models):
        self.workDir = workDir
        self.proteins = proteins
        self.lig_map = lig_map
        self.models = models

    def execute(self):
        self.lig_eicss = self.ligands_eicss()
        self.writeXML_ligands()

    def ligands_eicss(self):
        eicss_lig = {}
        for ligand in self.lig_map:
            if ligand.emdb_id not in eicss_lig:
                eicss_lig[ligand.emdb_id] = {}
            if ligand.emdb_id not in eicss_lig[ligand.emdb_id]:
                eicss_lig[ligand.emdb_id][ligand.sample_id] = ligand.__dict__
            else:
                eicss_lig[ligand.emdb_id][ligand.sample_id] += ligand.__dict__
        return eicss_lig


    def writeXML_ligands(self):
        component_DB = set()
        id_list = set()
        headerXML = EICSS.eicss()
        DBs_list = EICSS.DBs_listType()
        list_macro_molecules = EICSS.list_macro_moleculesType()
        list_crossRefDBs = EICSS.list_crossRefDBsType()

        for em_id, val in self.lig_eicss.items():
            headerXML.set_EMDB_ID(em_id)
            for samp_id in val.keys():
                lig_copies = val.get(samp_id, {}).get('lig_copies')
                lig_name = val.get(samp_id, {}).get('lig_name')
                chembl_id = val.get(samp_id, {}).get('chembl_id')
                chebi_id = val.get(samp_id, {}).get('chebi_id')
                drugbank_id = val.get(samp_id, {}).get('drugbank_id')
                provenance = val.get(samp_id, {}).get('provenance')


                macro_molecule_annotation = EICSS.macro_molecule_annotationType()
                macro_molecule_annotation.set_macro_kind("%s" % "ligand")
                macro_molecule_annotation.set_macro_ID(int(samp_id))
                macro_molecule_annotation.set_macro_copies(int(lig_copies))
                macro_molecule_annotation.set_macro_name("%s" % lig_name)
                list_macro_molecules.add_macro_molecule_annotation(macro_molecule_annotation)
                if chembl_id:
                    if "CHEMBL" not in component_DB:
                        DB = EICSS.DBType()
                        DB.set_DB_source("%s" % "CHEMBL")
                        DB.set_DB_version("%s" % "4.2.0")
                        DBs_list.add_DB(DB)
                    if chebi_id not in id_list:
                        crossRefDB = EICSS.crossRefDBType()
                        crossRefDB.set_DB_source("%s" % "ChEMBL")
                        crossRefDB.set_provenance("%s" % provenance)
                        crossRefDB.set_DB_accession_ID("%s" % chembl_id)
                        list_crossRefDBs.add_crossRefDB(crossRefDB)
                        macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                id_list.add(chembl_id)
                component_DB.add("CHEMBL")
                if chebi_id:
                    if "CHEBI" not in component_DB:
                        DB = EICSS.DBType()
                        DB.set_DB_source("%s" % "CHEBI")
                        DB.set_DB_version("%s" % "15.21")
                        DBs_list.add_DB(DB)
                    if chebi_id not in id_list:
                        crossRefDB = EICSS.crossRefDBType()
                        crossRefDB.set_DB_source("%s" % "ChEBI")
                        crossRefDB.set_provenance("%s" % provenance)
                        crossRefDB.set_DB_accession_ID("%s" % chebi_id)
                        list_crossRefDBs.add_crossRefDB(crossRefDB)
                        macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                id_list.add(chebi_id)
                component_DB.add("CHEBI")
                if drugbank_id:
                    if "DRUGBANK" not in component_DB:
                        DB = EICSS.DBType()
                        DB.set_DB_source("%s" % "DRUGBANK")
                        DB.set_DB_version("%s" % "2021.03.30")
                        DBs_list.add_DB(DB)
                    if drugbank_id not in id_list:
                        crossRefDB = EICSS.crossRefDBType()
                        crossRefDB.set_DB_source("%s" % "DrugBank")
                        crossRefDB.set_provenance("%s" % provenance)
                        crossRefDB.set_DB_accession_ID("%s" % drugbank_id)
                        list_crossRefDBs.add_crossRefDB(crossRefDB)
                        macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                id_list.add(drugbank_id)
                component_DB.add("DRUGBANK")

            headerXML.set_DBs_list(DBs_list)
            headerXML.set_sample_annotation(list_macro_molecules)

            xmlFile = os.path.join(self.workDir, em_id + "_eicss.xml")
            with open(xmlFile, 'w') as f:
                headerXML.export(f, 0, name_='eicss')