import os, re
import itertools
from EMICSS import EMICSS

class EmicssXML:
    "Writing annotations to output xml file according to the EMDB_EMICSS.xsd schema "

    def __init__(self, workDir, unip_map, cpx_map, lig_map, mw_map):
        self.workDir = workDir
        self.unip_map = unip_map
        self.cpx_map = cpx_map
        self.lig_map = lig_map
        self.mw_map = mw_map

    def execute(self):
        self.emicss_annotation = self.dict_emicss()
        self.writeXML_ligands()

    def dict_emicss(self):
        "Converts dictionary individual annotation to deeply nested dictionary of all the added annotations"
        emicss_dict = {}

        for mw in self.mw_map:
            if mw.emdb_id not in emicss_dict:
                emicss_dict[mw.emdb_id] = {}
            if mw.emdb_id not in emicss_dict[mw.emdb_id]:
                emicss_dict[mw.emdb_id][mw.pdb_id] = mw.__dict__
            else:
                emicss_dict[mw.emdb_id][mw.pdb_id] += mw.__dict__
        for unip in self.unip_map:
            if unip.emdb_id not in emicss_dict:
                emicss_dict[unip.emdb_id] = {}
            if unip.emdb_id not in emicss_dict[unip.emdb_id]:
                emicss_dict[unip.emdb_id][unip.uniprot_id] = unip.__dict__
            else:
                emicss_dict[unip.emdb_id][unip.uniprot_id] += unip.__dict__
        for emcpx in self.cpx_map:
            if emcpx:
                for cpx in emcpx.cpx_list:
                    lcpx = ["emdb_id", emcpx.emdb_id, "sample_id", emcpx.sample_id, "sample_name", emcpx.sample_name, "sample_copies",
                            emcpx.sample_copies, "cpx_id", cpx.cpx_id, "cpx_name", cpx.name, "provenance", emcpx.provenance, "score",
                            emcpx.score]
                    dcpx = dict(itertools.zip_longest(*[iter(lcpx)] * 2, fillvalue=""))
                    if emcpx.emdb_id not in emicss_dict:
                        emicss_dict[emcpx.emdb_id] = {}
                    if emcpx.emdb_id not in emicss_dict[emcpx.emdb_id]:
                        emicss_dict[emcpx.emdb_id][cpx.cpx_id] = dcpx
                    else:
                        emicss_dict[emcpx.emdb_id][cpx.cpx_id] += dcpx
        for ligand in self.lig_map:
            if ligand.emdb_id not in emicss_dict:
                emicss_dict[ligand.emdb_id] = {}
            if ligand.emdb_id not in emicss_dict[ligand.emdb_id]:
                emicss_dict[ligand.emdb_id][ligand.sample_id] = ligand.__dict__
            else:
                emicss_dict[ligand.emdb_id][ligand.sample_id] += ligand.__dict__

        return emicss_dict

    def writeXML_ligands(self):
        "Write added annotation to individual EMICSS file"
        # print(self.emicss_annotation)
        for em_id, val in self.emicss_annotation.items():
            all_DB = set()
            cp_id = set()
            headerXML = EMICSS.emicss()
            DBs_list = EMICSS.DBs_listType()
            models_list = EMICSS.models_listType()
            list_macro_molecules = EMICSS.list_macro_moleculesType()

            headerXML.set_EMDB_ID(em_id)

            for samp_id in val.keys():
                # print(em_id, samp_id)
                if samp_id is not None:
                    if (samp_id.isalnum() and not samp_id.isalpha() and not samp_id.isnumeric()):
                        if len(samp_id) == 4:
                            pdb_id = val.get(samp_id, {}).get('pdb_id')
                            assembly = val.get(samp_id, {}).get('assembly')
                            mw = val.get(samp_id, {}).get('molecular_weight')
                            if pdb_id:
                                if "PDBe" not in all_DB:
                                    DB = EMICSS.DBType()
                                    DB.set_DB_source("%s" % "PDBe")
                                    DB.set_DB_version("%s" % "2.0")
                                    DBs_list.add_DB(DB)
                            all_DB.add("PDBe")
                            model_annotation = EMICSS.model_annotationType()
                            model_annotation.set_PDBID("%s" % pdb_id)
                            model_annotation.set_assemblies(int(assembly))
                            model_annotation.set_weight(round(mw, 2))
                            model_annotation.set_units("%s" % "Da")
                            model_annotation.set_provenance("%s" % "PDBe")
                            models_list.add_model_annotation(model_annotation)

                        if len(samp_id) != 4:
                            list_crossRefDBs = EMICSS.list_crossRefDBsType()
                            sample_id = val.get(samp_id, {}).get('sample_id')
                            sample_copies = val.get(samp_id, {}).get('sample_copies')
                            sample_name = val.get(samp_id, {}).get('sample_name')
                            uniprot_id = val.get(samp_id, {}).get('uniprot_id')
                            uni_provenance = val.get(samp_id, {}).get('provenance')

                            macro_molecule_annotation = EMICSS.macro_molecule_annotationType()
                            macro_molecule_annotation.set_macro_kind("%s" % "protein")
                            macro_molecule_annotation.set_macro_ID(int(sample_id))
                            macro_molecule_annotation.set_macro_copies(int(sample_copies))
                            macro_molecule_annotation.set_macro_name("%s" % sample_name)
                            list_macro_molecules.add_macro_molecule_annotation(macro_molecule_annotation)
                            if uniprot_id:
                                if "UNIPROT" not in all_DB:
                                    DB = EMICSS.DBType()
                                    DB.set_DB_source("%s" % "UNIPROT")
                                    DB.set_DB_version("%s" % "2021.02")
                                    DBs_list.add_DB(DB)
                                crossRefDB = EMICSS.crossRefDBType()
                                crossRefDB.set_DB_source("%s" % "UNIPROT")
                                crossRefDB.set_provenance("%s" % uni_provenance)
                                crossRefDB.set_DB_accession_ID("%s" % uniprot_id)
                                list_crossRefDBs.add_crossRefDB(crossRefDB)
                                macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                            all_DB.add("UNIPROT")
                    if samp_id.isnumeric():
                        list_crossRefDBs = EMICSS.list_crossRefDBsType()
                        lig_copies = val.get(samp_id, {}).get('lig_copies')
                        lig_name = val.get(samp_id, {}).get('lig_name')
                        HET = val.get(samp_id, {}).get('HET')
                        chembl_id = val.get(samp_id, {}).get('chembl_id')
                        chebi_id = val.get(samp_id, {}).get('chebi_id')
                        drugbank_id = val.get(samp_id, {}).get('drugbank_id')
                        provenance = val.get(samp_id, {}).get('provenance')

                        macro_molecule_annotation = EMICSS.macro_molecule_annotationType()
                        macro_molecule_annotation.set_macro_kind("%s" % "ligand")
                        macro_molecule_annotation.set_macro_ID(int(samp_id))
                        macro_molecule_annotation.set_macro_CCD_ID("%s" % HET)
                        macro_molecule_annotation.set_macro_copies(int(lig_copies))
                        macro_molecule_annotation.set_macro_name("%s" % lig_name)
                        list_macro_molecules.add_macro_molecule_annotation(macro_molecule_annotation)

                        if chembl_id:
                            if "CHEMBL" not in all_DB:
                                DB = EMICSS.DBType()
                                DB.set_DB_source("%s" % "CHEMBL")
                                DB.set_DB_version("%s" % "4.2.0")
                                DBs_list.add_DB(DB)
                            crossRefDB = EMICSS.crossRefDBType()
                            crossRefDB.set_DB_source("%s" % "ChEMBL")
                            crossRefDB.set_provenance("%s" % provenance)
                            crossRefDB.set_DB_accession_ID("%s" % chembl_id)
                            list_crossRefDBs.add_crossRefDB(crossRefDB)
                            macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                        all_DB.add("CHEMBL")
                        if chebi_id:
                            if "CHEBI" not in all_DB:
                                DB = EMICSS.DBType()
                                DB.set_DB_source("%s" % "CHEBI")
                                DB.set_DB_version("%s" % "15.21")
                                DBs_list.add_DB(DB)
                            crossRefDB = EMICSS.crossRefDBType()
                            crossRefDB.set_DB_source("%s" % "ChEBI")
                            crossRefDB.set_provenance("%s" % provenance)
                            crossRefDB.set_DB_accession_ID("%s" % chebi_id)
                            list_crossRefDBs.add_crossRefDB(crossRefDB)
                            macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                        all_DB.add("CHEBI")
                        if drugbank_id:
                            if "DRUGBANK" not in all_DB:
                                DB = EMICSS.DBType()
                                DB.set_DB_source("%s" % "DRUGBANK")
                                DB.set_DB_version("%s" % "2021.03.30")
                                DBs_list.add_DB(DB)
                            crossRefDB = EMICSS.crossRefDBType()
                            crossRefDB.set_DB_source("%s" % "DrugBank")
                            crossRefDB.set_provenance("%s" % provenance)
                            crossRefDB.set_DB_accession_ID("%s" % drugbank_id)
                            list_crossRefDBs.add_crossRefDB(crossRefDB)
                            macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                        all_DB.add("DRUGBANK")

                    cpx_samp_id = re.search(r'%s\-\d+' % "CPX", samp_id)
                    cpx_sample_id = val.get(samp_id, {}).get('sample_id')
                    cpx_sample_copies = val.get(samp_id, {}).get('sample_copies')
                    cpx_sample_name = val.get(samp_id, {}).get('sample_name')
                    cpx_id = val.get(samp_id, {}).get('cpx_id')
                    cpx_name = val.get(samp_id, {}).get('cpx_name')
                    cpx_provenance = val.get(samp_id, {}).get('provenance')
                    cpx_score = val.get(samp_id, {}).get('score')
                    macro_molecule_annotation = EMICSS.macro_molecule_annotationType()
                    list_crossRefDBs = EMICSS.list_crossRefDBsType()
                    if cpx_samp_id is not None:
                        if cpx_sample_id not in cp_id:
                            crossRefDB = EMICSS.crossRefDBType()
                            macro_molecule_annotation.set_macro_kind("%s" % "complex")
                            macro_molecule_annotation.set_macro_ID(int(cpx_sample_id))
                            macro_molecule_annotation.set_macro_copies(int(cpx_sample_copies))
                            macro_molecule_annotation.set_macro_name("%s" % cpx_sample_name)
                            list_macro_molecules.add_macro_molecule_annotation(macro_molecule_annotation)
                            if "COMPLEX PORTAL" not in all_DB:
                                DB = EMICSS.DBType()
                                DB.set_DB_source("%s" % "COMPLEX PORTAL")
                                DB.set_DB_version("%s" % "1")
                                DBs_list.add_DB(DB)

                            crossRefDB.set_name("%s" % cpx_name)
                            crossRefDB.set_DB_source("%s" % "COMPLEX PORTAL")
                            crossRefDB.set_provenance("%s" % cpx_provenance)
                            crossRefDB.set_DB_accession_ID("%s" % cpx_id)
                            crossRefDB.set_score(float(cpx_score))
                            list_crossRefDBs.add_crossRefDB(crossRefDB)
                            macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                        if cpx_sample_id in cp_id:
                            crossRefDB = EMICSS.crossRefDBType()
                            crossRefDB.set_name("%s" % cpx_name)
                            crossRefDB.set_DB_source("%s" % "COMPLEX PORTAL")
                            crossRefDB.set_provenance("%s" % cpx_provenance)
                            crossRefDB.set_DB_accession_ID("%s" % cpx_id)
                            crossRefDB.set_score(float(cpx_score))
                            list_crossRefDBs.add_crossRefDB(crossRefDB)
                            macro_molecule_annotation.set_list_crossRefDBs(list_crossRefDBs)
                        cp_id.add(cpx_sample_id)
                        all_DB.add("COMPLEX PORTAL")

            headerXML.set_DBs_list(DBs_list)
            headerXML.set_sample_annotation(list_macro_molecules)
            headerXML.set_models_list(models_list)

            xmlFile = os.path.join(self.workDir, em_id + "_emicss.xml")
            with open(xmlFile, 'w') as f:
                headerXML.export(f, 0, name_='emicss')
