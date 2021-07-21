import os, re
import itertools
from EMICSS import EMICSS

class EmicssXML:
    """
    Writing annotations to output xml file according to the EMdb_EMICSS.xsd schema
    """

    def __init__(self, workDir, unip_map, cpx_map, lig_map, mw_map, sw_map, empiar_map, pmc_map):
        self.workDir = workDir
        self.unip_map = unip_map
        self.cpx_map = cpx_map
        self.lig_map = lig_map
        self.mw_map = mw_map
        self.sw_map = sw_map
        self.empiar_map = empiar_map
        self.pmc_map = pmc_map

    def execute(self):
        self.emicss_annotation = self.dict_emicss()
        self.writeXML_emicss()

    def dict_emicss(self):
        """
        Converts dictionary individual annotation to deeply nested dictionary of all the added annotations
        """

        emicss_dict = {}

        if self.mw_map:
            for mw in self.mw_map:
                if mw.emdb_id not in emicss_dict:
                    emicss_dict[mw.emdb_id] = {}
                if mw.emdb_id not in emicss_dict[mw.emdb_id]:
                    emicss_dict[mw.emdb_id][mw.pdb_id] = mw.__dict__
                else:
                    emicss_dict[mw.emdb_id][mw.pdb_id] += mw.__dict__

        if self.sw_map:
            for sw in self.sw_map:
                if sw.emdb_id not in emicss_dict:
                    emicss_dict[sw.emdb_id] = {}
                if sw.emdb_id not in emicss_dict[sw.emdb_id]:
                    emicss_dict[sw.emdb_id][sw.method] = sw.__dict__
                else:
                    emicss_dict[sw.emdb_id][sw.method] += sw.__dict__

        if self.empiar_map:
            for empiar in self.empiar_map:
                if empiar.emdb_id not in emicss_dict:
                    emicss_dict[empiar.emdb_id] = {}
                if empiar.emdb_id not in emicss_dict[empiar.emdb_id]:
                    emicss_dict[empiar.emdb_id][empiar.empiar_id] = empiar.__dict__
                else:
                    emicss_dict[empiar.emdb_id][empiar.empiar_id] += empiar.__dict__

        if self.unip_map:
            for unip in self.unip_map:
                if unip.emdb_id not in emicss_dict:
                    emicss_dict[unip.emdb_id] = {}
                if unip.emdb_id not in emicss_dict[unip.emdb_id]:
                    emicss_dict[unip.emdb_id][unip.uniprot_id] = unip.__dict__
                else:
                    emicss_dict[unip.emdb_id][unip.uniprot_id] += unip.__dict__

        if self.cpx_map:
            for emcpx in self.cpx_map:
                if emcpx:
                    for cpx in emcpx.cpx_list:
                        if emcpx.emdb_id not in emicss_dict.keys():
                            emicss_dict[emcpx.emdb_id] = {}
                        if emcpx.sample_id not in emicss_dict[emcpx.emdb_id].keys():
                            emicss_dict[emcpx.emdb_id][emcpx.sample_id] = {}
                            ind = 0
                        lcpx = ["supra_name", emcpx.supra_name, "sample_copies", emcpx.sample_copies, "cpx_id" + "_" + str(ind),
                                cpx.cpx_id, "cpx_name" + "_" + str(ind), cpx.name, "provenance" + "_" + str(ind),
                                emcpx.provenance, "score" + "_" + str(ind), emcpx.score]
                        dcpx = dict(itertools.zip_longest(*[iter(lcpx)] * 2, fillvalue=""))
                        for k in dcpx.keys():
                            emicss_dict[emcpx.emdb_id][emcpx.sample_id][k] = dcpx[k]
                        ind = ind + 1
                    emicss_dict[emcpx.emdb_id][emcpx.sample_id]["ind"] = ind

        if self.lig_map:
            for ligand in self.lig_map:
                if ligand.emdb_id not in emicss_dict:
                    emicss_dict[ligand.emdb_id] = {}
                if ligand.emdb_id not in emicss_dict[ligand.emdb_id]:
                    emicss_dict[ligand.emdb_id][ligand.sample_id] = ligand.__dict__
                else:
                    emicss_dict[ligand.emdb_id][ligand.sample_id] += ligand.__dict__

        if self.pmc_map:
            for pmc in self.pmc_map:
                if pmc.emdb_id not in emicss_dict:
                    emicss_dict[pmc.emdb_id] = {}
                if pmc.emdb_id not in emicss_dict[pmc.emdb_id]:
                    emicss_dict[pmc.emdb_id]["PMC"] = pmc.__dict__
                else:
                    emicss_dict[pmc.emdb_id]["PMC"] += pmc.__dict__

        return emicss_dict

    def writeXML_emicss(self):
        """
        Create and write added annotations to individual EMICSS file for every EMDB entry
        """
        # print(self.emicss_annotation)
        for em_id, val in self.emicss_annotation.items():
            all_db = set()
            headerXML = EMICSS.emicss()
            dbs = EMICSS.dbsType()
            molecular_weight = EMICSS.molecular_weightType()
            empiars = EMICSS.empiarsType()
            models = EMICSS.modelsType()
            weights = EMICSS.weightsType()
            europe_pmc = EMICSS.europe_pmcType()
            sample = EMICSS.sampleType()
            macromolecules = EMICSS.macromoleculesType()

            headerXML.set_emdb_id(em_id)
            supramolecules = None
            for samp_id in val.keys():
                if samp_id is not None:
                    if re.search(r'%s\-\d+' % "EMPIAR", samp_id):
                        self.EMICSS_empiar(val, samp_id, all_db, dbs, empiars)
                    if samp_id == "theoretical" or samp_id == "experimental":
                        self.EMICSS_weight(val, samp_id, weights)
                    if samp_id == "PMC":
                        self.EMICSS_PMC(val, samp_id, all_db, dbs, europe_pmc)
                    if (samp_id.isalnum() and not samp_id.isalpha() and not samp_id.isnumeric()):
                        if len(samp_id) == 4:
                            self.EMICSS_Pdbe(val, samp_id, all_db, dbs, models)
                        if len(samp_id) != 4:
                            self.EMICSS_uniprot(val, samp_id, all_db, dbs, macromolecules)
                    if samp_id.isnumeric():
                        self.EMICSS_ligands(val, samp_id, all_db, dbs, macromolecules)
                    if re.search(r'%s\_\d+' % em_id, samp_id):
                        supramolecules = self.EMICSS_CPX(val, samp_id, all_db, dbs)

            headerXML.set_dbs(dbs)
            headerXML.set_empiars(empiars)
            headerXML.set_europe_pmc(europe_pmc)
            molecular_weight.set_models(models)
            molecular_weight.set_weights(weights)
            headerXML.set_molecular_weight(molecular_weight)
            if supramolecules:
                sample.set_supramolecules(supramolecules)
            sample.set_macromolecules(macromolecules)
            headerXML.set_sample(sample)

            xmlFile = os.path.join(self.workDir, "emicss", em_id + "_emicss.xml")
            with open(xmlFile, 'w') as f:
                headerXML.export(f, 0, name_='emicss')

    def EMICSS_empiar(self, val, samp_id, all_db, dbs, empiars):
        "Adding EMPIAR_ID to EMICSS"

        empiar_id = val.get(samp_id, {}).get('empiar_id')
        if empiar_id:
            if "EMPIAR" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "EMPIAR")
                db.set_db_version("%s" % "2.0")
                dbs.add_db(db)
        all_db.add("EMPIAR")
        empiar = EMICSS.empiarType()
        empiar.set_empiar_id("%s" % empiar_id)
        empiar.set_provenance("%s" % "AUTHOR")
        empiars.add_empiar(empiar)

    def EMICSS_PMC(self, val, samp_id, all_db, dbs, europe_pmc):
        """
        Adding PUBMED_ID, DOI and ISSN to EMICSS
        """
        pmedid = val.get(samp_id, {}).get('pmedid')
        doi = val.get(samp_id, {}).get('doi')
        issn = val.get(samp_id, {}).get('issn')
        title = val.get(samp_id, {}).get('title')
        provenance = val.get(samp_id, {}).get('provenance')
        link = val.get(samp_id, {}).get('url')
        if pmedid or doi or issn:
            if "EuropePMC" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "EuropePMC")
                db.set_db_version("%s" % "2.0")
                dbs.add_db(db)
        all_db.add("EuropePMC")
        pmc = EMICSS.pmcType()
        if pmedid:
            pmc.set_pubmed_id("%s" % pmedid)
            europe_pmc.set_pmc_link("%s" % link)
        if doi:
            pmc.set_doi("%s" % doi)
        if issn:
            pmc.set_issn("%s" % issn)
        if pmedid or doi or issn:
            pmc.set_provenance("%s" % provenance)
        europe_pmc.add_pmc(pmc)

    def EMICSS_Pdbe(self, val, samp_id, all_db, dbs, models):
        """
        Adding pdb_id and calulated assembly weight annotations to EMICSS
        """
        pdb_id = val.get(samp_id, {}).get('pdb_id')
        assembly = val.get(samp_id, {}).get('assembly')
        mw = val.get(samp_id, {}).get('molecular_weight')
        if pdb_id:
            if "PDBe" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "PDBe")
                db.set_db_version("%s" % "2.0")
                dbs.add_db(db)
        all_db.add("PDBe")
        model = EMICSS.modelType()
        model.set_pdb_id("%s" % pdb_id)
        model.set_assemblies(int(assembly))
        model.set_weight(round(mw, 2))
        model.set_units("%s" % "Da")
        model.set_provenance("%s" % "PDBe")
        models.add_model(model)

    def EMICSS_weight(self, val, samp_id, weights):
        "Adding author provided calulated total sample weight annotations to EMICSS"

        kind = val.get(samp_id, {}).get('kind')
        th_weight = val.get(samp_id, {}).get('sample_th_weight')
        th_units = val.get(samp_id, {}).get('th_unit')
        exp_weight = val.get(samp_id, {}).get('sample_exp_weight')
        exp_units = val.get(samp_id, {}).get('exp_unit')
        if th_weight:
            weight = EMICSS.weightType()
            weight.set_kind("%s" % kind)
            weight.set_method("%s" % "theoretical")
            weight.set_weight(round(th_weight, 3))
            weight.set_units("%s" % th_units)
            weight.set_provenance("%s" % "AUTHOR")
            weights.add_weight(weight)
        if exp_weight:
            weight = EMICSS.weightType()
            weight.set_kind("%s" % kind)
            weight.set_method("%s" % "experimental")
            weight.set_weight(round(exp_weight, 3))
            weight.set_units("%s" % exp_units)
            weight.set_provenance("%s" % "AUTHOR")
            weights.add_weight(weight)

    def EMICSS_uniprot(self, val, samp_id, all_db, dbs, macromolecules):
        "Adding UNIPROT annotation to EMICSS"

        cross_ref_dbs = EMICSS.cross_ref_dbsType()
        sample_id = val.get(samp_id, {}).get('sample_id')
        sample_copies = val.get(samp_id, {}).get('sample_copies')
        name = val.get(samp_id, {}).get('sample_name')
        uniprot_id = val.get(samp_id, {}).get('uniprot_id')
        uni_provenance = val.get(samp_id, {}).get('provenance')

        macromolecule = EMICSS.macromoleculeType()
        macromolecule.set_kind("%s" % "protein")
        macromolecule.set_id(int(sample_id))
        macromolecule.set_copies(int(sample_copies))
        macromolecule.set_name("%s" % name)
        macromolecules.add_macromolecule(macromolecule)
        if uniprot_id:
            if "UNIPROT" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "UNIPROT")
                db.set_db_version("%s" % "2021.02")
                dbs.add_db(db)
            cross_ref_db = EMICSS.cross_ref_dbType()
            cross_ref_db.set_db_source("%s" % "UNIPROT")
            cross_ref_db.set_provenance("%s" % uni_provenance)
            cross_ref_db.set_db_accession_id("%s" % uniprot_id)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            macromolecule.set_cross_ref_dbs(cross_ref_dbs)
        all_db.add("UNIPROT")

    def EMICSS_ligands(self, val, samp_id, all_db, dbs, macromolecules):
        "Adding components annotation to EMICSS"

        cross_ref_dbs = EMICSS.cross_ref_dbsType()
        lig_copies = val.get(samp_id, {}).get('lig_copies')
        lig_name = val.get(samp_id, {}).get('lig_name')
        HET = val.get(samp_id, {}).get('HET')
        chembl_id = val.get(samp_id, {}).get('chembl_id')
        chebi_id = val.get(samp_id, {}).get('chebi_id')
        drugbank_id = val.get(samp_id, {}).get('drugbank_id')
        provenance = val.get(samp_id, {}).get('provenance')

        macromolecule = EMICSS.macromoleculeType()
        macromolecule.set_kind("%s" % "ligand")
        macromolecule.set_id(int(samp_id))
        macromolecule.set_ccd_id("%s" % HET)
        macromolecule.set_copies(int(lig_copies))
        macromolecule.set_name("%s" % lig_name)
        macromolecules.add_macromolecule(macromolecule)
        if chembl_id:
            if "CHEMBL" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "CHEMBL")
                db.set_db_version("%s" % "4.2.0")
                dbs.add_db(db)
            cross_ref_db = EMICSS.cross_ref_dbType()
            cross_ref_db.set_db_source("%s" % "ChEMBL")
            cross_ref_db.set_provenance("%s" % provenance)
            cross_ref_db.set_db_accession_id("%s" % chembl_id)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            macromolecule.set_cross_ref_dbs(cross_ref_dbs)
        all_db.add("CHEMBL")
        if chebi_id:
            if "CHEBI" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "CHEBI")
                db.set_db_version("%s" % "15.21")
                dbs.add_db(db)
            cross_ref_db = EMICSS.cross_ref_dbType()
            cross_ref_db.set_db_source("%s" % "ChEBI")
            cross_ref_db.set_provenance("%s" % provenance)
            cross_ref_db.set_db_accession_id("%s" % chebi_id)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            macromolecule.set_cross_ref_dbs(cross_ref_dbs)
        all_db.add("CHEBI")
        if drugbank_id:
            if "DRUGBANK" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "DRUGBANK")
                db.set_db_version("%s" % "2021.03.30")
                dbs.add_db(db)
            cross_ref_db = EMICSS.cross_ref_dbType()
            cross_ref_db.set_db_source("%s" % "DrugBank")
            cross_ref_db.set_provenance("%s" % provenance)
            cross_ref_db.set_db_accession_id("%s" % drugbank_id)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            macromolecule.set_cross_ref_dbs(cross_ref_dbs)
        all_db.add("DRUGBANK")

    def EMICSS_CPX(self, val, samp_id, all_db, dbs):
        """
        Adding complex ids to EMICSS
        """
        supramolecules = EMICSS.supramoleculesType()
        cp_id = set()
        cpx_samp_id = samp_id.split("_")[1]
        cpx_sample_copies = val.get(samp_id, {}).get('sample_copies')
        cpx_sample_name = val.get(samp_id, {}).get('supra_name')
        ind = val.get(samp_id, {}).get('ind')
        supramolecule = EMICSS.supramoleculeType()
        cross_ref_dbs = EMICSS.cross_ref_dbsType()
        for x in range(ind):
            c_id = "cpx_id_"+str(x)
            cpx_id = val.get(samp_id, {}).get(c_id)
            c_name = "cpx_name_"+str(x)
            c_provenance = "provenance_"+str(x)
            c_score = "score_"+str(x)
            cpx_name = val.get(samp_id, {}).get(c_name)
            cpx_provenance = val.get(samp_id, {}).get(c_provenance)
            cpx_score = val.get(samp_id, {}).get(c_score)
            if cpx_samp_id is not None:
                if cpx_samp_id not in cp_id:
                    cross_ref_db = EMICSS.cross_ref_dbType()
                    supramolecule.set_kind("%s" % "complex")
                    supramolecule.set_id(int(cpx_samp_id))
                    supramolecule.set_copies(int(cpx_sample_copies))
                    supramolecule.set_name("%s" % cpx_sample_name)
                    if "COMPLEX PORTAL" not in all_db:
                        db = EMICSS.dbType()
                        db.set_db_source("%s" % "COMPLEX PORTAL")
                        db.set_db_version("%s" % "1")
                        dbs.add_db(db)

                    cross_ref_db.set_name("%s" % cpx_name)
                    cross_ref_db.set_db_source("%s" % "COMPLEX PORTAL")
                    cross_ref_db.set_provenance("%s" % cpx_provenance)
                    cross_ref_db.set_db_accession_id("%s" % cpx_id)
                    cross_ref_db.set_score(float(cpx_score))
                    cross_ref_dbs.add_cross_ref_db(cross_ref_db)
                if cpx_samp_id in cp_id:
                    cross_ref_db = EMICSS.cross_ref_dbType()
                    cross_ref_db.set_name("%s" % cpx_name)
                    cross_ref_db.set_db_source("%s" % "COMPLEX PORTAL")
                    cross_ref_db.set_provenance("%s" % cpx_provenance)
                    cross_ref_db.set_db_accession_id("%s" % cpx_id)
                    cross_ref_db.set_score(float(cpx_score))
                    cross_ref_dbs.add_cross_ref_db(cross_ref_db)

                cp_id.add(cpx_samp_id)
                all_db.add("COMPLEX PORTAL")
        supramolecule.set_cross_ref_dbs(cross_ref_dbs)
        # print(f'Adding supramolecule {supramolecule}: {supramolecule.__dict__}')
        supramolecules.add_supramolecule(supramolecule)
        return supramolecules
