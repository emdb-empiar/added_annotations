import os, re
from pathlib import Path
from EMICSS.EMICSS import *

class EmicssXML:
    """
    Writing annotations to output xml file according to the EMDB_EMICSS.xsd schema
    # TODO: db_source should be a enum field
    # TODO: weight_info has an attribute named method that seems to be not used anymore
    """

    def __init__(self, workDir, version_list):
        self.workDir = workDir
        self.version_list = version_list

    def write(self, mapping_list):
        """
        Create and write added annotations to individual EMICSS file for every EMDB entry
        """
        emdb_id = mapping_list['HEADER'].emdb_id
        #TODO: Schema version can not be hard coded here. It must follow the version of the xsd file used to generate pymodels
        headerXML = EMICSS.emicss(emdb_id=emdb_id, schema_version="0.9.0")
        dbs = EMICSS.dbsType()
        entry_ref_dbs = EMICSS.entry_ref_dbsType()
        weights = EMICSS.weightsType()
        sample = EMICSS.sampleType()
        macromolecules = EMICSS.macromoleculesType()
        supramolecules = EMICSS.supramoleculesType()
        all_db = set()
        uniq_id = set()

    
        if "EMPIAR" in mapping_list:
            empiar_objects = mapping_list['EMPIAR']
            if len(empiar_objects) > 0:
                all_db.add("EMPIAR")
                for empiar in empiar_objects:
                    emp_ref_db_obj = entry_ref_dbType(db_source="EMPIAR", accession_id=empiar.empiar_id, provenance="EMDB")
                    entry_ref_dbs.add_entry_ref_db(emp_ref_db_obj)
        # MW calculated from header
        if "WEIGHT" in mapping_list: 
            mw = mapping_list["WEIGHT"]
            if mw.overall_mw > 0:
                all_db.add("EMDB")
                mw_info_obj = weight_infoType(weight=mw.overall_mw, unit=mw.units, provenance=mw.provenance)
                weights.add_weight_info(mw_info_obj)
        # MW calculated from assemblies
        if "MODEL" in mapping_list:
            pdb_objects = mapping_list["MODEL"]
            if len(pdb_objects) > 0:
                all_db.add("PDBe")
                for pdb in pdb_objects:
                    pdb_info_obj = weight_infoType(pdb_id=pdb.pdb_id, assemblies=pdb.assembly, weight=pdb.molecular_weight, unit="Da", provenance="PDBe")
                    weights.add_weight_info(pdb_info_obj)
        if "CITATION" in mapping_list:
            citation = mapping_list["CITATION"]
            if citation.pmedid:
                primary_citation = primary_citationType()
                authors_obj = authorsType()
                all_db.add("PUBMED")
                pm_citation_obj = ref_citationType(db_source="PUBMED", accession_id=citation.pmedid, provenance=citation.provenance_pm)
                primary_citation.add_ref_citation(pm_citation_obj)
                if citation.pmcid:
                    all_db.add("PUBMED CENTRAL")
                    pmc_citation_obj = ref_citationType(db_source="PUBMED CENTRAL", accession_id=citation.pmcid, provenance=citation.provenance_pmc)
                    primary_citation.add_ref_citation(pmc_citation_obj)
                if citation.issn:
                    all_db.add("ISSN")
                    issn_citation_obj = ref_citationType(db_source="ISSN", accession_id=citation.issn, provenance=citation.provenance_issn)
                    primary_citation.add_ref_citation(issn_citation_obj)
                if citation.doi:
                    primary_citation.set_doi(citation.doi)
                    primary_citation.set_provenance(citation.provenance_doi)
                for author in citation.authors:
                    orcid = author.orcid if author.orcid else None
                    author_obj = authorType(name=author.name, orcid_id=orcid, order=author.order, provenance=author.provenance)
                    authors_obj.add_author(author_obj)
                primary_citation.set_authors(authors)






        # TODO: Use all_db to create dbType and DBVersion




            





        for vers in range(0, len(self.version_list), 2):
            if self.version_list[vers] == "uniprot":
                self.unip_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "pdbe":
                self.pdbe_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "pdbekb":
                self.pdbekb_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "alphafold":
                self.alphafold_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "pfam":
                self.pfam_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "go":
                self.go_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "interpro":
                self.interpro_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "cath":
                self.cath_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "scop":
                self.scop_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "scop2":
                self.scop2_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "empiar":
                self.empiar_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "chembl":
                self.chembl_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "chebi":
                self.chebi_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "drugbank":
                self.drugbank_vers = self.version_list[vers + 1]
            if self.version_list[vers] == "cpx":
                self.cpx_vers = self.version_list[vers + 1]



        emicss_supramolecules = None
        emicss_macromol = None
        emicss_ligand = None
        emicss_empiar = None
        emicss_model = None
        emicss_weight = None
        emicss_authors = None
        emicss_primary_citation = None

        for samp_id in val.keys():
            if samp_id is not None:
                if samp_id == "PMC":
                    emicss_primary_citation = self.EMICSS_PMC(val, samp_id, all_db, dbs, primary_citation)
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
        entry = f"EMD_{entry_id}"
        headerXML.set_emdb_id("%s" % entry)
        headerXML.set_schema_version("0.9.0")

        if len(all_db) != 0:
            headerXML.set_dbs(dbs)
        if emicss_empiar:
            headerXML.set_entry_ref_dbs(entry_ref_dbs)
        if emicss_primary_citation:
            headerXML.set_primary_citation(primary_citation)
        if emicss_model or emicss_weight:
            headerXML.set_weights(weights)
        if emicss_supramolecules:
            sample.set_supramolecules(supramolecules)
        if emicss_macromol or emicss_ligand:
            sample.set_macromolecules(macromolecules)
        if emicss_supramolecules or emicss_macromol or emicss_ligand:
            headerXML.set_sample(sample)

        output_path = os.path.join(self.workDir, "emicss_xml")
        Path(output_path).mkdir(parents=True, exist_ok=True)
        xmlFile = os.path.join(output_path, "emd-" + entry_id + "_emicss.xml")
        with open(xmlFile, 'w') as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            headerXML.export(f, 0, name_='emicss',  namespacedef_='schemaLocation=""')
            # headerXML.export(f, 0, name_='emicss',
            #                  namespacedef_='xmlns="http://pdbe.org/empiar" '
            #                                'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
            #                                'xsi:schemaLocation="https://ftp.ebi.ac.uk/pub/databases/emtest/empiar/schema/empiar.xsd"')


def parse_empiar(self, objects):
    """
    Parse EMPIAR objects into EMICSS data model
    """
    if len(objects) > 0:
        for empiar in objects:
            emp_db_type = entry_ref_dbType(db_source="EMPIAR", accession_id=empiar.empiar_id, provenance="EMDB")


    empiar_id = val.get(samp_id, {}).get('empiar_id')
    if "EMDB" not in all_db:
        db = EMICSS.dbType()
        db.set_db_source("%s" % "EMDB")
        dbs.add_db(db)
    all_db.add("EMDB")
    if empiar_id:
        if "EMPIAR" not in all_db:
            db = EMICSS.dbType()
            db.set_db_source("%s" % "EMPIAR")
            db.set_db_version("%s" % self.empiar_vers)
            dbs.add_db(db)
    all_db.add("EMPIAR")
    entry_ref_db = EMICSS.entry_ref_dbType()
    entry_ref_db.set_db_source("%s" % "EMPIAR")
    entry_ref_db.set_accession_id("%s" % empiar_id)
    entry_ref_db.set_provenance("%s" % "EMDB")
    entry_ref_dbs.add_entry_ref_db(entry_ref_db)
    return entry_ref_dbs

def EMICSS_PMC(self, val, samp_id, all_db, dbs, primary_citation):
    """
    Adding PUBMED_ID, DOI and ISSN to EMICSS
    """
    authors = EMICSS.authorsType()
    pmedid = val.get(samp_id, {}).get('pmedid')
    pmcid = val.get(samp_id, {}).get('pmcid')
    pub_doi = val.get(samp_id, {}).get('doi')
    issn = val.get(samp_id, {}).get('issn')
    orcid_ids = val.get(samp_id, {}).get('orcid_ids')
    provenance_pm = val.get(samp_id, {}).get('provenance_pm')
    provenance_pmc = val.get(samp_id, {}).get('provenance_pmc')
    provenance_doi = val.get(samp_id, {}).get('provenance_doi')
    provenance_issn = val.get(samp_id, {}).get('provenance_issn')
    if "EMDB" not in all_db:
        db = EMICSS.dbType()
        db.set_db_source("%s" % "EMDB")
        dbs.add_db(db)
    all_db.add("EMDB")
    if pmedid:
        if "PUBMED" not in all_db:
            db = EMICSS.dbType()
            db.set_db_source("%s" % "PUBMED")
            dbs.add_db(db)
        ref_citation = EMICSS.ref_citationType()
        ref_citation.set_db_source("%s" % "PUBMED")
        ref_citation.set_accession_id("%s" % pmedid)
        ref_citation.set_provenance("%s" % provenance_pm)
        primary_citation.add_ref_citation(ref_citation)
        all_db.add("PUBMED")
    if pmcid:
        if "PUBMED CENTRAL" not in all_db:
            db = EMICSS.dbType()
            db.set_db_source("%s" % "PUBMED CENTRAL")
            dbs.add_db(db)
        ref_citation = EMICSS.ref_citationType()
        ref_citation.set_db_source("%s" % "PUBMED CENTRAL")
        ref_citation.set_accession_id("%s" % pmcid)
        ref_citation.set_provenance("%s" % provenance_pmc)
        primary_citation.add_ref_citation(ref_citation)
        all_db.add("PUBMED CENTRAL")
    if issn:
        if "ISSN" not in all_db:
            db = EMICSS.dbType()
            db.set_db_source("%s" % "ISSN")
            dbs.add_db(db)
        ref_citation = EMICSS.ref_citationType()
        ref_citation.set_db_source("%s" % "ISSN")
        ref_citation.set_accession_id("%s" % issn)
        ref_citation.set_provenance("%s" % provenance_issn)
        primary_citation.add_ref_citation(ref_citation)
        all_db.add("ISSN")
    if pub_doi:
        primary_citation.set_doi("%s" % pub_doi)
        primary_citation.set_provenance("%s" % provenance_doi)
    if orcid_ids:
        ind = orcid_ids.get('ind')
        if "ORCID" not in all_db:
            db = EMICSS.dbType()
            db.set_db_source("%s" % "ORCID")
            dbs.add_db(db)
        all_db.add("ORCID")
        for x in range(int(ind)):
            auth_name = "name_" + str(x)
            name = orcid_ids.get(auth_name)
            auth_id = "id_" + str(x)
            id = orcid_ids.get(auth_id)
            auth_ord = "order_" + str(x)
            order = orcid_ids.get(auth_ord)
            prov = "provenance_orcid_" + str(x)
            provenance_orcid = orcid_ids.get(prov)
            author = EMICSS.authorType()
            author.set_name("%s" % name)
            if id != "N/A":
                author.set_orcid_id("%s" % id)
            if int(order) != 0:
                author.set_order(int(order))
            author.set_provenance("%s" % provenance_orcid)
            authors.add_author(author)
        primary_citation.set_authors(authors)

    return primary_citation

def EMICSS_Pdbe(self, val, samp_id, all_db, dbs, weights):
    """
    Adding pdb_id and calulated assembly weight annotations to EMICSS
    """
    pdb_id = val.get(samp_id, {}).get('pdb_id')
    assembly = val.get(samp_id, {}).get('assembly')
    mw = val.get(samp_id, {}).get('molecular_weight')
    if "EMDB" not in all_db:
        db = EMICSS.dbType()
        db.set_db_source("%s" % "EMDB")
        dbs.add_db(db)
    all_db.add("EMDB")
    if pdb_id:
        if "PDBe" not in all_db:
            db = EMICSS.dbType()
            db.set_db_source("%s" % "PDBe")
            dbs.add_db(db)
        all_db.add("PDBe")
        weight_info = EMICSS.weight_infoType()
        weight_info.set_pdb_id("%s" % pdb_id)
        weight_info.set_assemblies(int(assembly))
        weight_info.set_weight("%s" % mw)
        weight_info.set_unit("%s" % "Da")
        weight_info.set_provenance("%s" % "PDBe")
        weights.add_weight_info(weight_info)
    return weights

def EMICSS_weight(self, val, all_db, dbs, samp_id, weights):
    "Adding author provided calulated total sample weight annotations to EMICSS"

    overall_mw = val.get(samp_id, {}).get('overall_mw')
    units = val.get(samp_id, {}).get('units')
    provenance = val.get(samp_id, {}).get('provenance')
    if overall_mw:
        if "EMDB" not in all_db:
            db = EMICSS.dbType()
            db.set_db_source("%s" % "EMDB")
            dbs.add_db(db)
        weight_info = EMICSS.weight_infoType()
        weight_info.set_weight("%s" % overall_mw)
        weight_info.set_unit("%s" % units)
        weight_info.set_provenance("%s" % provenance)
        weights.add_weight_info(weight_info)
    all_db.add("EMDB")
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
    macromolecule.set_provenance("%s" % "EMDB")
    macromolecule.set_name("%s" % name)
    if "EMDB" not in all_db:
        db = EMICSS.dbType()
        db.set_db_source("%s" % "EMDB")
        dbs.add_db(db)
    all_db.add("EMDB")
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
        cross_ref_db.set_provenance("%s" % uni_provenance)
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
        ipr_uniprot_start = "ipr_start_" + str(x)
        IPR_uniprot_start = val.get(samp_id, {}).get(ipr_uniprot_start)
        ipr_uniprot_end = "ipr_unp_end_" + str(x)
        IPR_uniprot_end = val.get(samp_id, {}).get(ipr_uniprot_end)
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
            if IPR_uniprot_start and IPR_uniprot_end is not None:
                cross_ref_db.set_uniprot_start(int(IPR_uniprot_start))
                cross_ref_db.set_uniprot_end(int(IPR_uniprot_end))
            cross_ref_db.set_provenance("%s" % IPR_provenance)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)

        pfam_id = "pfam_id_" + str(x)
        PFAM_id = val.get(samp_id, {}).get(pfam_id)
        pfam_namespace = "pfam_namespace_" + str(x)
        PFAM_namespace = val.get(samp_id, {}).get(pfam_namespace)
        pfam_uniprot_start = "pfam_unp_start_" + str(x)
        PFAM_uniprot_start = val.get(samp_id, {}).get(pfam_uniprot_start)
        pfam_uniprot_end = "pfam_unp_end_" + str(x)
        PFAM_uniprot_end = val.get(samp_id, {}).get(pfam_uniprot_end)
        pfam_provenance = "pfam_provenance_" + str(x)
        PFAM_provenance = val.get(samp_id, {}).get(pfam_provenance)
        if PFAM_id:
            if "PFAM" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "PFAM")
                try:
                    if self.pfam_vers:
                        db.set_db_version("%s" % self.pfam_vers)
                except AttributeError as e:
                    print("Pfam version not collected.", e)
                dbs.add_db(db)
            all_db.add("PFAM")

            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_db_source("%s" % "PFAM")
            cross_ref_db.set_accession_id("%s" % PFAM_id)
            cross_ref_db.set_name("%s" % PFAM_namespace)
            if PFAM_uniprot_start and PFAM_uniprot_end is not None:
                cross_ref_db.set_uniprot_start(int(PFAM_uniprot_start))
                cross_ref_db.set_uniprot_end(int(PFAM_uniprot_end))
            cross_ref_db.set_provenance("%s" % PFAM_provenance)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)

        cath_id = "cath_id_" + str(x)
        CATH_id = val.get(samp_id, {}).get(cath_id)
        cath_uniprot_start = "cath_unp_start_" + str(x)
        CATH_uniprot_start = val.get(samp_id, {}).get(cath_uniprot_start)
        cath_uniprot_end = "cath_unp_end_" + str(x)
        CATH_uniprot_end = val.get(samp_id, {}).get(cath_uniprot_end)
        cath_provenance = "cath_provenance_" + str(x)
        CATH_provenance = val.get(samp_id, {}).get(cath_provenance)
        if CATH_id:
            if "CATH" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "CATH")
                db.set_db_version("%s" % self.cath_vers)
                dbs.add_db(db)
            all_db.add("CATH")

            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_db_source("%s" % "CATH")
            cross_ref_db.set_accession_id("%s" % CATH_id)
            if CATH_uniprot_start and CATH_uniprot_end is not None:
                cross_ref_db.set_uniprot_start(int(CATH_uniprot_start))
                cross_ref_db.set_uniprot_end(int(CATH_uniprot_end))
                if CATH_uniprot_start > CATH_uniprot_end:
                    print("Something wrong as UniProt start is greater than UniProt end")
            cross_ref_db.set_provenance("%s" % CATH_provenance)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)

        scop_id = "scop_id_" + str(x)
        SCOP_id = val.get(samp_id, {}).get(scop_id)
        scop_uniprot_start = "scop_unp_start_" + str(x)
        SCOP_uniprot_start = val.get(samp_id, {}).get(scop_uniprot_start)
        scop_uniprot_end = "scop_unp_end_" + str(x)
        SCOP_uniprot_end = val.get(samp_id, {}).get(scop_uniprot_end)
        scop_provenance = "scop_provenance_" + str(x)
        SCOP_provenance = val.get(samp_id, {}).get(scop_provenance)
        if SCOP_id:
            if "SCOP" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "SCOP")
                db.set_db_version("%s" % self.scop_vers)
                dbs.add_db(db)
            all_db.add("SCOP")

            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_db_source("%s" % "SCOP")
            cross_ref_db.set_accession_id("%s" % SCOP_id)
            if SCOP_uniprot_start and SCOP_uniprot_end is not None:
                cross_ref_db.set_uniprot_start(int(SCOP_uniprot_start))
                cross_ref_db.set_uniprot_end(int(SCOP_uniprot_end))
                if SCOP_uniprot_start > SCOP_uniprot_end:
                    print("Something wrong as UniProt start is greater than UniProt end")
            cross_ref_db.set_provenance("%s" % SCOP_provenance)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)

        scop2_id = "scop2_id_" + str(x)
        SCOP2_id = val.get(samp_id, {}).get(scop2_id)
        scop2_uniprot_start = "scop2_unp_start_" + str(x)
        SCOP2_uniprot_start = val.get(samp_id, {}).get(scop2_uniprot_start)
        scop2_uniprot_end = "scop2_unp_end_" + str(x)
        SCOP2_uniprot_end = val.get(samp_id, {}).get(scop2_uniprot_end)
        scop2_provenance = "scop2_provenance_" + str(x)
        SCOP2_provenance = val.get(samp_id, {}).get(scop2_provenance)
        if SCOP2_id:
            if "SCOP2" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "SCOP2")
                db.set_db_version("%s" % self.scop2_vers)
                dbs.add_db(db)
            all_db.add("SCOP2")

            cross_ref_db = EMICSS.cross_ref_db()
            cross_ref_db.set_db_source("%s" % "SCOP2")
            cross_ref_db.set_accession_id("%s" % SCOP2_id)
            if SCOP2_uniprot_start and SCOP2_uniprot_end is not None:
                cross_ref_db.set_uniprot_start(int(SCOP2_uniprot_start))
                cross_ref_db.set_uniprot_end(int(SCOP2_uniprot_end))
                if SCOP2_uniprot_start > SCOP2_uniprot_end:
                    print("Something wrong as UniProt start is greater than UniProt end")
            cross_ref_db.set_provenance("%s" % SCOP2_provenance)
            cross_ref_dbs.add_cross_ref_db(cross_ref_db)

        kb_link = "kb_link_" + str(x)
        KB_link = val.get(samp_id, {}).get(kb_link)
        kb_provenance = "kb_provenance_" + str(x)
        KB_provenance = val.get(samp_id, {}).get(kb_provenance)
        if KB_link:
            if "PDBe-KB" not in all_db:
                db = EMICSS.dbType()
                db.set_db_source("%s" % "PDBe-KB")
                # db.set_db_version("%s" % self.pdbekb_vers)
                dbs.add_db(db)
            all_db.add("PDBe-KB")

            if uniprot_id not in uniq_id:
                cross_ref_db = EMICSS.cross_ref_db()
                cross_ref_db.set_name("%s" % name.lower())
                cross_ref_db.set_db_source("%s" % "PDBe-KB")
                cross_ref_db.set_accession_id("%s" % uniprot_id)
                # cross_ref_db.set_link("%s" % KB_link)
                cross_ref_db.set_provenance("%s" % KB_provenance)
                cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            uniq_id.add(uniprot_id)

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
    macromolecule.set_provenance("%s" % "EMDB")
    macromolecule.set_name("%s" % lig_name)
    macromolecules.add_macromolecule(macromolecule)
    if "EMDB" not in all_db:
        db = EMICSS.dbType()
        db.set_db_source("%s" % "EMDB")
        dbs.add_db(db)
    all_db.add("EMDB")
    if "PDBe-CCD" not in all_db:
        db = EMICSS.dbType()
        db.set_db_source("%s" % "PDBe-CCD")
        dbs.add_db(db)
    all_db.add("PDBe-CCD")
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
    if "EMDB" not in all_db:
        db = EMICSS.dbType()
        db.set_db_source("%s" % "EMDB")
        dbs.add_db(db)
    all_db.add("EMDB")
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
                supramolecule.set_provenance("%s" % "EMDB")
                supramolecule.set_name("%s" % cpx_sample_name)
                if "COMPLEX PORTAL" not in all_db:
                    db = EMICSS.dbType()
                    db.set_db_source("%s" % "COMPLEX PORTAL")
                    dbs.add_db(db)

                cross_ref_db.set_name("%s" % cpx_name)
                cross_ref_db.set_db_source("%s" % "COMPLEX PORTAL")
                cross_ref_db.set_accession_id("%s" % cpx_id)
                cross_ref_db.set_score(float(cpx_score))
                cross_ref_dbs.add_cross_ref_db(cross_ref_db)
            if cpx_samp_id in cp_id:
                cross_ref_db = EMICSS.cross_ref_db()
                cross_ref_db.set_name("%s" % cpx_name)
                cross_ref_db.set_db_source("%s" % "COMPLEX PORTAL")
                cross_ref_db.set_provenance("%s" % cpx_provenance)
                cross_ref_db.set_accession_id("%s" % cpx_id)
                cross_ref_db.set_score(float(cpx_score))
                cross_ref_dbs.add_cross_ref_db(cross_ref_db)

            cp_id.add(cpx_samp_id)
            all_db.add("COMPLEX PORTAL")
    supramolecule.set_cross_ref_dbs(cross_ref_dbs)
    supramolecules.add_supramolecule(supramolecule)
    return supramolecules