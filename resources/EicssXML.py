import os
from EICSS import EICSS

class EicssXML:
    "Writing annotations to output xml file according to the EMDB_EICSS.xsd schema "

    def __init__(self, workDir, proteins, ligands, models):
        self.workDir = workDir
        self.proteins = proteins
        self.ligands = ligands
        self.models = models

    def execute(self):
        self.writeXML_ligands()

    def writeXML_ligands(self):
        headerXML = EICSS.eicss()
        DBs_list = EICSS.DBs_listType()
        #DB_source = EICSS.DB_source_type()
        list_macro_molecules = EICSS.list_macro_moleculesType()
        sample_annotation = EICSS.sample_annotationType()
        macro_molecule_annotation = EICSS.macro_molecule_annotationType()
        list_crossRefDBs = EICSS.list_crossRefDBsType()
        crossRefDB = EICSS.crossRefDBType()

        for ligand in self.ligands:
            print(type(ligand))
            print("ligand:", ligand)
            headerXML.set_EMDB_ID(ligand.emdb_id)

            DB = EICSS.DBType()
            if ligand.chembl_id:
                print(ligand.chembl_id)
                DB.set_DB_source("%s" % "CHEMBL")
                DB.set_DB_version("%s" % "4.2.0")
                print("chembl_id db:", DB.__dict__)
                DBs_list.add_DB(DB)
            if ligand.chebi_id:
                DB.set_DB_source("%s" % "CHEBI")
                DB.set_DB_version("%s" % "15.21")
                print("chebi_id db:", DB.__dict__)
                DBs_list.add_DB(DB)
            if ligand.drugbank_id:
                DB.set_DB_source("%s" % "DRUGBANK")
                DB.set_DB_version("%s" % "2021.03.30")
                print("drugbank_id db:", DB.__dict__)
                DBs_list.add_DB(DB)
            print("LIST", DBs_list.__dict__)
            headerXML.set_DBs_list(DBs_list)


                # macro_molecule_annotation.set_macro_kind("%s" % "ligand")
                # macro_molecule_annotation.set_macro_ID(int(ligand.sample_id))
                # macro_molecule_annotation.set_macro_copies(int(ligand.lig_copies))
                # macro_molecule_annotation.set_macro_name("%s" % ligand.lig_name)
                # crossRefDB.set_DB_source("%s" % "ChEMBL")
                # crossRefDB.set_provenance("%s" % ligand.provenance)
                # crossRefDB.set_DB_accession_ID("%s" % ligand.chembl_id)
                # list_crossRefDBs.add_crossRefDB(crossRefDB)
                # macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                # list_macro_molecules.add_macro_molecule_annotation(macro_molecule_annotation)
                # headerXML.set_sample_annotation(list_macro_molecules)
            #print(ligand.emdb_id, ligand
        xmlFile = os.path.join(self.workDir, ligand.emdb_id + "_eicss.xml")
        with open(xmlFile, 'w') as f:
           headerXML.export(f, 0, name_='eicss')