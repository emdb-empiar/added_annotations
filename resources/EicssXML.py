import os, re
from EICSS import EICSS

class EicssXML:
    "Writing annotations to output xml file according to the EMDB_EICSS.xsd schema "

    def __init__(self, workDir, proteins, lig_map, mw_map):
        self.workDir = workDir
        self.proteins = proteins
        self.lig_map = lig_map
        self.mw_map = mw_map

    def execute(self):
        self.eicss_annotation = self.dict_eicss()
        self.writeXML_ligands()

    def dict_eicss(self):
        eicss_dict = {}
        for ligand in self.lig_map:
            if ligand.emdb_id not in eicss_dict:
                eicss_dict[ligand.emdb_id] = {}
            if ligand.emdb_id not in eicss_dict[ligand.emdb_id]:
                eicss_dict[ligand.emdb_id][ligand.sample_id] = ligand.__dict__
            else:
                eicss_dict[ligand.emdb_id][ligand.sample_id] += ligand.__dict__
        for mw in self.mw_map:
            if mw.emdb_id not in eicss_dict:
                eicss_dict[mw.emdb_id] = {}
            if mw.emdb_id not in eicss_dict[mw.emdb_id]:
                eicss_dict[mw.emdb_id][mw.pdb_id] = mw.__dict__
            else:
                eicss_dict[mw.emdb_id][mw.pdb_id] += mw.__dict__
        return eicss_dict


    def writeXML_ligands(self):
        # print(self.eicss_annotation)
        for em_id, val in self.eicss_annotation.items():
            all_DB = set()
            headerXML = EICSS.eicss()
            DBs_list = EICSS.DBs_listType()
            models_list = EICSS.models_listType()
            list_macro_molecules = EICSS.list_macro_moleculesType()

            headerXML.set_EMDB_ID(em_id)

            for samp_id in val.keys():
                if samp_id.isalnum() and not samp_id.isalpha() and not samp_id.isnumeric():
                    pdb_id = val.get(samp_id, {}).get('pdb_id')
                    assembly = val.get(samp_id, {}).get('assembly')
                    mw = val.get(samp_id, {}).get('molecular_weight')
                    if pdb_id:
                        if "PDBe" not in all_DB:
                            DB = EICSS.DBType()
                            DB.set_DB_source("%s" % "PDBe")
                            DB.set_DB_version("%s" % "2.0")
                            DBs_list.add_DB(DB)
                    all_DB.add("PDBe")
                    model_annotation = EICSS.model_annotationType()
                    model_annotation.set_PDBID("%s" % pdb_id)
                    model_annotation.set_assemblies(int(assembly))
                    model_annotation.set_weight(float(mw))
                    model_annotation.set_units("%s" % "Da")
                    model_annotation.set_provenance("%s" % "PDBe")
                    models_list.add_model_annotation(model_annotation)
                if samp_id.isnumeric():
                    list_crossRefDBs = EICSS.list_crossRefDBsType()
                    lig_copies = val.get(samp_id, {}).get('lig_copies')
                    lig_name = val.get(samp_id, {}).get('lig_name')
                    HET = val.get(samp_id, {}).get('HET')
                    chembl_id = val.get(samp_id, {}).get('chembl_id')
                    chebi_id = val.get(samp_id, {}).get('chebi_id')
                    drugbank_id = val.get(samp_id, {}).get('drugbank_id')
                    provenance = val.get(samp_id, {}).get('provenance')

                    macro_molecule_annotation = EICSS.macro_molecule_annotationType()
                    macro_molecule_annotation.set_macro_kind("%s" % "ligand")
                    macro_molecule_annotation.set_macro_ID(int(samp_id))
                    macro_molecule_annotation.set_macro_CCD_ID("%s" % HET)
                    macro_molecule_annotation.set_macro_copies(int(lig_copies))
                    macro_molecule_annotation.set_macro_name("%s" % lig_name)
                    list_macro_molecules.add_macro_molecule_annotation(macro_molecule_annotation)

                    if chembl_id:
                        if "CHEMBL" not in all_DB:
                            DB = EICSS.DBType()
                            DB.set_DB_source("%s" % "CHEMBL")
                            DB.set_DB_version("%s" % "4.2.0")
                            DBs_list.add_DB(DB)
                        crossRefDB = EICSS.crossRefDBType()
                        crossRefDB.set_DB_source("%s" % "ChEMBL")
                        crossRefDB.set_provenance("%s" % provenance)
                        crossRefDB.set_DB_accession_ID("%s" % chembl_id)
                        list_crossRefDBs.add_crossRefDB(crossRefDB)
                        macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                    all_DB.add("CHEMBL")
                    if chebi_id:
                        if "CHEBI" not in all_DB:
                            DB = EICSS.DBType()
                            DB.set_DB_source("%s" % "CHEBI")
                            DB.set_DB_version("%s" % "15.21")
                            DBs_list.add_DB(DB)
                        crossRefDB = EICSS.crossRefDBType()
                        crossRefDB.set_DB_source("%s" % "ChEBI")
                        crossRefDB.set_provenance("%s" % provenance)
                        crossRefDB.set_DB_accession_ID("%s" % chebi_id)
                        list_crossRefDBs.add_crossRefDB(crossRefDB)
                        macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                    all_DB.add("CHEBI")
                    if drugbank_id:
                        if "DRUGBANK" not in all_DB:
                            DB = EICSS.DBType()
                            DB.set_DB_source("%s" % "DRUGBANK")
                            DB.set_DB_version("%s" % "2021.03.30")
                            DBs_list.add_DB(DB)
                        crossRefDB = EICSS.crossRefDBType()
                        crossRefDB.set_DB_source("%s" % "DrugBank")
                        crossRefDB.set_provenance("%s" % provenance)
                        crossRefDB.set_DB_accession_ID("%s" % drugbank_id)
                        list_crossRefDBs.add_crossRefDB(crossRefDB)
                        macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                    all_DB.add("DRUGBANK")

            headerXML.set_DBs_list(DBs_list)

            headerXML.set_sample_annotation(list_macro_molecules)
            headerXML.set_models_list(models_list)
            xmlFile = os.path.join(self.workDir, em_id + "_eicss.xml")
            with open(xmlFile, 'w') as f:
                headerXML.export(f, 0, name_='eicss')