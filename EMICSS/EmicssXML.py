import os, re
from pathlib import Path
from EMICSS import EMICSS

class EmicssXML:
    """
    Writing annotations to output xml file according to the EMDB_EMICSS.xsd schema
    """

    def __init__(self, workDir, emicss_annotation, version_list):
        self.workDir = workDir
        self.emicss_annotation = emicss_annotation
        self.version_list = version_list

    def execute(self):
        self.writeXML_emicss(self.emicss_annotation, self.version_list)

    def writeXML_emicss(self, emicss_annotation, version_list):
        """
        Create and write added annotations to individual EMICSS file for every EMDB entry
        """
        # print(emicss_annotation)
        for em_id, val in emicss_annotation.items():
            all_db = set()
            for vers in range(0, len(version_list), 2):
                if version_list[vers] == "uniprot":
                    self.unip_vers = version_list[vers + 1]
                if version_list[vers] == "pdbe":
                    self.pdbe_vers = version_list[vers + 1]
                if version_list[vers] == "pdbekb":
                    self.pdbekb_vers = version_list[vers + 1]
                if version_list[vers] == "alphafold":
                    self.alphafold_vers = version_list[vers + 1]
                if version_list[vers] == "pfam":
                    self.pfam_vers = version_list[vers + 1]
                if version_list[vers] == "go":
                    self.go_vers = version_list[vers + 1]
                if version_list[vers] == "interpro":
                    self.interpro_vers = version_list[vers + 1]
                if version_list[vers] == "empiar":
                    self.empiar_vers = version_list[vers + 1]
                if version_list[vers] == "chembl":
                    self.chembl_vers = version_list[vers + 1]
                if version_list[vers] == "chebi":
                    self.chebi_vers = version_list[vers + 1]
                if version_list[vers] == "drugbank":
                    self.drugbank_vers = version_list[vers + 1]
                if version_list[vers] == "cpx":
                    self.cpx_vers = version_list[vers + 1]

            headerXML = EMICSS.emicss()
            dbs = EMICSS.dbsType()
            cross_ref_dbs = EMICSS.cross_ref_dbsType()
            citations = EMICSS.citationsType()
            weights = EMICSS.weightsType()
            sample = EMICSS.sampleType()
            macromolecules = EMICSS.macromoleculesType()
            supramolecules = EMICSS.supramoleculesType()

            emicss_supramolecules = None
            emicss_macromol = None
            emicss_ligand = None
            emicss_empiar = None
            emicss_model = None
            emicss_weight = None
            emicss_pmc = None
            emicss_citation = None

            for samp_id in val.keys():
                if samp_id is not None:
                    if re.search(r'%s\-\d+' % "EMPIAR", samp_id):
                        emicss_empiar = self.EMICSS_empiar(val, samp_id, all_db, dbs, cross_ref_dbs)
                    if samp_id == "theoretical" or samp_id == "experimental":
                        emicss_weight = self.EMICSS_weight(val, samp_id, weights)
                    if samp_id == "PMC":
                        emicss_pmc, emicss_citation = self.EMICSS_PMC(val, samp_id, all_db, dbs, cross_ref_dbs, citations)
                    if (samp_id.isalnum() and not samp_id.isalpha() and not samp_id.isnumeric()):
                        if len(samp_id) == 4:
                            emicss_model = self.EMICSS_Pdbe(val, samp_id, all_db, dbs, weights)
                        if len(samp_id) != 4:
                            emicss_macromol = self.EMICSS_proteins(val, samp_id, all_db, dbs, macromolecules)
                    if re.search(r'%s\_\d+' % "ligand", samp_id):
                        emicss_ligand = self.EMICSS_ligands(val, samp_id, all_db, dbs, macromolecules)
                    if re.search(r'%s\_\d+' % em_id, samp_id):
                        emicss_supramolecules = self.EMICSS_CPX(val, samp_id, all_db, dbs, supramolecules)

            entry_id = em_id.split("-")[1]
            entry = f"emd_{entry_id}"
            headerXML.set_emdb_id("%s" % entry)
            headerXML.set_schema_version("1.0.0")

            if len(all_db) != 0:
                headerXML.set_dbs(dbs)
            if emicss_empiar or emicss_model or emicss_pmc:
                headerXML.set_cross_ref_dbs(cross_ref_dbs)
            if emicss_citation:
                headerXML.set_citations(citations)
            if emicss_weight:
                headerXML.set_weights(weights)
            if emicss_supramolecules:
                sample.set_supramolecules(supramolecules)
            if emicss_macromol or emicss_ligand:
                sample.set_macromolecules(macromolecules)
            if emicss_supramolecules or emicss_macromol or emicss_ligand:
                headerXML.set_sample(sample)

            output_path = os.path.join(self.workDir, "emicss")
            Path(output_path).mkdir(parents=True, exist_ok=True)
            xmlFile = os.path.join(output_path, "emd-" + entry_id + "_emicss.xml")
            with open(xmlFile, 'w') as f:
                headerXML.export(f, 0, name_='emicss')
                # headerXML.export(f, 0, name_='emicss',
                #                  namespacedef_='xmlns="http://pdbe.org/empiar" '
                #                                'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                #                                'xsi:schemaLocation="https://ftp.ebi.ac.uk/pub/databases/emtest/empiar/schema/empiar.xsd"')

    def EMICSS_empiar(self, val, samp_id, all_db, dbs, cross_ref_dbs):
        "Adding EMPIAR_ID to EMICSS"

        empiar_id = val.get(samp_id, {}).get('empiar_id')
        if empiar_id:
            if "EMPIAR" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "EMPIAR")
                db.set_db_version("%s" % self.empiar_vers)
                dbs.add_db(db)
        all_db.add("EMPIAR")
        cross_ref_db = EMICSS.cross_ref_db()
        cross_ref_db.set_db_source("%s" % "EMPIAR")
        cross_ref_db.set_accession_id("%s" % empiar_id)
        cross_ref_db.set_provenance("%s" % "AUTHOR")
        cross_ref_dbs.add_cross_ref_db(cross_ref_db)
        return cross_ref_dbs

    def EMICSS_PMC(self, val, samp_id, all_db, dbs, cross_ref_dbs, citations):
        """
        Adding PUBMED_ID, DOI and ISSN to EMICSS
        """
        pmedid = val.get(samp_id, {}).get('pmedid')
        pmcid = val.get(samp_id, {}).get('pmcid')
        doi = val.get(samp_id, {}).get('doi')
        issn = val.get(samp_id, {}).get('issn')
        provenance_pm = val.get(samp_id, {}).get('provenance_pm')
        provenance_pmc = val.get(samp_id, {}).get('provenance_pmc')
        provenance_doi = val.get(samp_id, {}).get('provenance_doi')
        if pmedid:
            if "PUBMED" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "PUBMED")
                dbs.add_db(db)
            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_db_source("%s" % "PUBMED")
            cross_ref_db.set_accession_id("%s" % pmedid)
            cross_ref_db.set_provenance("%s" % provenance_pm)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            all_db.add("PUBMED")
        if pmcid:
            if "PUBMED CENTRAL" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "PUBMED CENTRAL")
                dbs.add_db(db)
            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_db_source("%s" % "PUBMED CENTRAL")
            cross_ref_db.set_accession_id("%s" % pmcid)
            cross_ref_db.set_provenance("%s" % provenance_pmc)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            all_db.add("PUBMED CENTRAL")
        if issn:
            if "ISSN" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "ISSN")
                dbs.add_db(db)
            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_db_source("%s" % "ISSN")
            cross_ref_db.set_accession_id("%s" % issn)
            cross_ref_db.set_provenance("%s" % "AUTHOR")
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            all_db.add("ISSN")
        if doi:
            citation = EMICSS.citationType()
            citation.set_doi("%s" % doi)
            citation.set_provenance("%s" % provenance_doi)
            citations.add_citation(citation)
        return cross_ref_dbs, citations

    def EMICSS_Pdbe(self, val, samp_id, all_db, dbs, weights):
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
                db.set_db_version("%s" % "38.21")
                dbs.add_db(db)
            all_db.add("PDBe")
            weight = EMICSS.weightType()
            weight.set_pdb_id("%s" % pdb_id)
            weight.set_assemblies(int(assembly))
            weight.set_weight("%s" % mw)
            weight.set_unit("%s" % "Da")
            weight.set_provenance("%s" % "PDBe")
            weights.add_weight(weight)
        return weights

    def EMICSS_weight(self, val, samp_id, weights):
        "Adding author provided calulated total sample weight annotations to EMICSS"

        th_weight = val.get(samp_id, {}).get('sample_th_weight')
        th_units = val.get(samp_id, {}).get('th_unit')
        exp_weight = val.get(samp_id, {}).get('sample_exp_weight')
        exp_units = val.get(samp_id, {}).get('exp_unit')
        if th_weight:
            weight = EMICSS.weightType()
            weight.set_method("%s" % "theoretical")
            weight.set_weight("%s" % th_weight)
            weight.set_unit("%s" % th_units)
            weight.set_provenance("%s" % "AUTHOR")
            weights.add_weight(weight)
        if exp_weight:
            weight = EMICSS.weightType()
            weight.set_method("%s" % "experimental")
            weight.set_weight("%s" % exp_weight)
            weight.set_unit("%s" % exp_units)
            weight.set_provenance("%s" % "AUTHOR")
            weights.add_weight(weight)
        return weights

    def EMICSS_proteins(self, val, samp_id, all_db, dbs, macromolecules):
        "Adding UNIPROT, GO, INTERPRO and PFAM annotations to EMICSS"

        macromolecule = EMICSS.macromoleculeType()
        cross_ref_dbs = EMICSS.cross_ref_dbsType()
        sample_id = val.get(samp_id, {}).get('sample_id')
        sample_copies = val.get(samp_id, {}).get('sample_copies')
        name = val.get(samp_id, {}).get('sample_name')
        uniprot_id = val.get(samp_id, {}).get('uniprot_id')
        uni_provenance = val.get(samp_id, {}).get('provenance')
        macromolecule.set_type("%s" % "protein")
        macromolecule.set_id(int(sample_id))
        macromolecule.set_copies(int(sample_copies))
        macromolecule.set_name("%s" % name)
        if uniprot_id:
            if "UNIPROT" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "UNIPROT")
                db.set_db_version("%s" % self.unip_vers)
                dbs.add_db(db)
            all_db.add("UNIPROT")
            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_name("%s" % name.lower())
            cross_ref_db.set_db_source("%s" % "UNIPROT")
            if not "+" in uni_provenance:
                cross_ref_db.set_provenance("%s" % uni_provenance)
            if "+" in uni_provenance:
                uni_provenance1 = uni_provenance.split(" + ")[0]
                uni_provenance2 = uni_provenance.split(" + ")[1]
                cross_ref_db.set_provenance1("%s" % uni_provenance1.strip())
                cross_ref_db.set_provenance2("%s" % uni_provenance2.strip())
            cross_ref_db.set_accession_id("%s" % uniprot_id)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)

        ind = val.get(samp_id, {}).get('ind')
        for x in range(ind):
            go_id = "id_" + str(x)
            GO_id = val.get(samp_id, {}).get(go_id)
            go_namespace = "namespace_" + str(x)
            GO_namespace = val.get(samp_id, {}).get(go_namespace)
            go_type = "type_" + str(x)
            GO_type = val.get(samp_id, {}).get(go_type)
            go_provenance = "provenance_" + str(x)
            GO_provenance = val.get(samp_id, {}).get(go_provenance)
            if GO_id:
                if "GO" not in all_db:
                    db = EMICSS.dbType()
                    db.set_db_source("%s" % "GO")
                    db.set_db_version("%s" % self.go_vers)
                    dbs.add_db(db)
                all_db.add("GO")

                cross_ref_db = EMICSS.cross_ref_db()
                cross_ref_db.set_db_source("%s" % "GO")
                cross_ref_db.set_accession_id("%s" % GO_id)
                cross_ref_db.set_name("%s" % GO_namespace)
                if GO_type == "P":
                    cross_ref_db.set_type("%s" % "biological_process")
                if GO_type == "C":
                    cross_ref_db.set_type("%s" % "cellular_component")
                if GO_type == "F":
                    cross_ref_db.set_type("%s" % "molecular_function")
                cross_ref_db.set_provenance("%s" % GO_provenance)
                cross_ref_dbs.add_cross_ref_db(cross_ref_db)

            ipr_id = "ipr_id_" + str(x)
            IPR_id = val.get(samp_id, {}).get(ipr_id)
            ipr_namespace = "ipr_namespace_" + str(x)
            IPR_namespace = val.get(samp_id, {}).get(ipr_namespace)
            ipr_provenance = "ipr_provenance_" + str(x)
            IPR_provenance = val.get(samp_id, {}).get(ipr_provenance)
            if IPR_id:
                if "INTERPRO" not in all_db:
                    db = EMICSS.dbType()
                    db.set_db_source("%s" % "INTERPRO")
                    db.set_db_version("%s" % self.interpro_vers)
                    dbs.add_db(db)
                all_db.add("INTERPRO")

                cross_ref_db = EMICSS.cross_ref_db()
                cross_ref_db.set_db_source("%s" % "INTERPRO")
                cross_ref_db.set_accession_id("%s" % IPR_id)
                cross_ref_db.set_name("%s" % IPR_namespace)
                cross_ref_db.set_provenance("%s" % IPR_provenance)
                cross_ref_dbs.add_cross_ref_db(cross_ref_db)

            pfam_id = "pfam_id_" + str(x)
            PFAM_id = val.get(samp_id, {}).get(pfam_id)
            pfam_namespace = "pfam_namespace_" + str(x)
            PFAM_namespace = val.get(samp_id, {}).get(pfam_namespace)
            pfam_provenance = "pfam_provenance_" + str(x)
            PFAM_provenance = val.get(samp_id, {}).get(pfam_provenance)
            if PFAM_id:
                if "PFAM" not in all_db:
                    db = EMICSS.dbType()
                    db.set_db_source("%s" % "PFAM")
                    db.set_db_version("%s" % self.pfam_vers)
                    dbs.add_db(db)
                all_db.add("PFAM")

                cross_ref_db = EMICSS.cross_ref_db()
                cross_ref_db.set_db_source("%s" % "PFAM")
                cross_ref_db.set_accession_id("%s" % PFAM_id)
                cross_ref_db.set_name("%s" % PFAM_namespace)
                cross_ref_db.set_provenance("%s" % PFAM_provenance)
                cross_ref_dbs.add_cross_ref_db(cross_ref_db)

            kb_link = "kb_link_" + str(x)
            KB_link = val.get(samp_id, {}).get(kb_link)
            kb_provenance = "kb_provenance_" + str(x)
            KB_provenance = val.get(samp_id, {}).get(kb_provenance)
            if KB_link:
                if "PDBe-KB" not in all_db:
                    db = EMICSS.dbType()
                    db.set_db_source("%s" % "PDBe-KB")
                    db.set_db_version("%s" % self.pdbekb_vers)
                    dbs.add_db(db)
                all_db.add("PDBe-KB")

                cross_ref_db = EMICSS.cross_ref_db()
                cross_ref_db.set_name("%s" % name.lower())
                cross_ref_db.set_db_source("%s" % "PDBe-KB")
                cross_ref_db.set_accession_id("%s" % uniprot_id)
                # cross_ref_db.set_link("%s" % KB_link)
                cross_ref_db.set_provenance("%s" % KB_provenance)
                cross_ref_dbs.add_cross_ref_db(cross_ref_db)

            af_link = "af_link_" + str(x)
            AF_link = val.get(samp_id, {}).get(af_link)
            af_provenance = "af_provenance_" + str(x)
            AF_provenance = val.get(samp_id, {}).get(af_provenance)
            if AF_link:
                if "ALPHAFOLD DB" not in all_db:
                    db = EMICSS.dbType()
                    db.set_db_source("%s" % "ALPHAFOLD DB")
                    db.set_db_version("%s" % self.alphafold_vers)
                    dbs.add_db(db)
                all_db.add("ALPHAFOLD DB")

                cross_ref_db = EMICSS.cross_ref_db()
                cross_ref_db.set_name("%s" % name.lower())
                cross_ref_db.set_db_source("%s" % "ALPHAFOLD DB")
                cross_ref_db.set_accession_id("%s" % uniprot_id)
                # cross_ref_db.set_link("%s" % AF_link)
                cross_ref_db.set_provenance("%s" % AF_provenance)
                cross_ref_dbs.add_cross_ref_db(cross_ref_db)

        macromolecule.set_cross_ref_dbs(cross_ref_dbs)
        macromolecules.add_macromolecule(macromolecule)
        return macromolecules

    def EMICSS_ligands(self, val, samp_id, all_db, dbs, macromolecules):
        "Adding components annotation to EMICSS"

        cross_ref_dbs = EMICSS.cross_ref_dbsType()
        lig_copies = val.get(samp_id, {}).get('lig_copies')
        lig_name = val.get(samp_id, {}).get('lig_name')
        sample_id = val.get(samp_id, {}).get('sample_id')
        HET = val.get(samp_id, {}).get('HET')
        chembl_id = val.get(samp_id, {}).get('chembl_id')
        chebi_id = val.get(samp_id, {}).get('chebi_id')
        drugbank_id = val.get(samp_id, {}).get('drugbank_id')
        provenance_chembl = val.get(samp_id, {}).get('provenance_chembl')
        provenance_chebi = val.get(samp_id, {}).get('provenance_chebi')
        provenance_drugbank = val.get(samp_id, {}).get('provenance_drugbank')

        macromolecule = EMICSS.macromoleculeType()
        macromolecule.set_type("%s" % "ligand")
        macromolecule.set_id(int(sample_id))
        macromolecule.set_ccd_id("%s" % HET)
        macromolecule.set_copies(int(lig_copies))
        macromolecule.set_name("%s" % lig_name)
        macromolecules.add_macromolecule(macromolecule)
        if chembl_id:
            if "CHEMBL" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "CHEMBL")
                db.set_db_version("%s" % self.chembl_vers)
                dbs.add_db(db)
            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_name("%s" % lig_name)
            cross_ref_db.set_db_source("%s" % "CHEMBL")
            cross_ref_db.set_provenance("%s" % provenance_chembl)
            cross_ref_db.set_accession_id("%s" % chembl_id)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            macromolecule.set_cross_ref_dbs(cross_ref_dbs)
        all_db.add("CHEMBL")
        if chebi_id:
            if "CHEBI" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "CHEBI")
                db.set_db_version("%s" % self.chebi_vers)
                dbs.add_db(db)
            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_name("%s" % lig_name)
            cross_ref_db.set_db_source("%s" % "CHEBI")
            cross_ref_db.set_provenance("%s" % provenance_chebi)
            cross_ref_db.set_accession_id("%s" % chebi_id)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            macromolecule.set_cross_ref_dbs(cross_ref_dbs)
        all_db.add("CHEBI")
        if drugbank_id:
            if "DRUGBANK" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "DRUGBANK")
                db.set_db_version("%s" % self.drugbank_vers)
                dbs.add_db(db)
            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_name("%s" % lig_name)
            cross_ref_db.set_db_source("%s" % "DRUGBANK")
            cross_ref_db.set_provenance("%s" % provenance_drugbank)
            cross_ref_db.set_accession_id("%s" % drugbank_id)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            macromolecule.set_cross_ref_dbs(cross_ref_dbs)
        all_db.add("DRUGBANK")
        return macromolecules

    def EMICSS_CPX(self, val, samp_id, all_db, dbs, supramolecules):
        """
        Adding complex ids to EMICSS
        """

        cp_id = set()
        cpx_samp_id = samp_id.split("_")[1]
        cpx_sample_copies = val.get(samp_id, {}).get('sample_copies')
        cpx_sample_name = val.get(samp_id, {}).get('supra_name')
        ind = val.get(samp_id, {}).get('ind')
        supramolecule = EMICSS.supramoleculeType()
        cross_ref_dbs = EMICSS.cross_ref_dbsType()
        for x in range(ind):
            c_id = "cpx_id_" + str(x)
            cpx_id = val.get(samp_id, {}).get(c_id)
            c_name = "cpx_name_" + str(x)
            c_provenance = "provenance_" + str(x)
            c_score = "score_" + str(x)
            cpx_name = val.get(samp_id, {}).get(c_name)
            cpx_provenance = val.get(samp_id, {}).get(c_provenance)
            cpx_score = val.get(samp_id, {}).get(c_score)
            if cpx_samp_id is not None:
                if cpx_samp_id not in cp_id:
                    cross_ref_db = EMICSS.cross_ref_db()
                    supramolecule.set_type("%s" % "complex")
                    supramolecule.set_id(int(cpx_samp_id))
                    supramolecule.set_copies(int(cpx_sample_copies))
                    supramolecule.set_name("%s" % cpx_sample_name)
                    if "COMPLEX PORTAL" not in all_db:
                        db = EMICSS.dbType()
                        db.set_db_source("%s" % "COMPLEX PORTAL")
                        db.set_db_version("%s" % self.cpx_vers)
                        dbs.add_db(db)

                    cross_ref_db.set_name("%s" % cpx_name)
                    cross_ref_db.set_db_source("%s" % "COMPLEX PORTAL")
                    if "+" in cpx_provenance:
                        cpx_provenance1 = cpx_provenance.split(" + ")[0]
                        cpx_provenance2 = cpx_provenance.split(" + ")[1]
                        cross_ref_db.set_provenance1("%s" % cpx_provenance1.strip())
                        cross_ref_db.set_provenance2("%s" % cpx_provenance2.strip())
                    if not "+" in cpx_provenance:
                        cross_ref_db.set_provenance("%s" % cpx_provenance)
                    cross_ref_db.set_accession_id("%s" % cpx_id)
                    cross_ref_db.set_score(float(cpx_score))
                    cross_ref_dbs.add_cross_ref_db(cross_ref_db)
                if cpx_samp_id in cp_id:
                    cross_ref_db = EMICSS.cross_ref_db()
                    cross_ref_db.set_name("%s" % cpx_name)
                    cross_ref_db.set_db_source("%s" % "COMPLEX PORTAL")
                    if "+" in cpx_provenance:
                        cpx_provenance1 = cpx_provenance.split("+")[0]
                        cpx_provenance2 = cpx_provenance.split("+")[1]
                        cross_ref_db.set_provenance1("%s" % cpx_provenance1.strip())
                        cross_ref_db.set_provenance2("%s" % cpx_provenance2.strip())
                    if not "+" in cpx_provenance:
                        cross_ref_db.set_provenance("%s" % cpx_provenance)
                    cross_ref_db.set_accession_id("%s" % cpx_id)
                    cross_ref_db.set_score(float(cpx_score))
                    cross_ref_dbs.add_cross_ref_db(cross_ref_db)

                cp_id.add(cpx_samp_id)
                all_db.add("COMPLEX PORTAL")
        supramolecule.set_cross_ref_dbs(cross_ref_dbs)
        supramolecules.add_supramolecule(supramolecule)
        return supramolecules