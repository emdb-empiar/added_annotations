import os
from pathlib import Path
from EMICSS.EMICSSGDS import *

class EmicssXML:
    """
    Writing annotations to output xml file according to the EMDB_EMICSS.xsd schema
    """

    def __init__(self, emdb_id, workDir, version_list):
        self.emdb_id = emdb_id
        self.workDir = workDir
        self.version_list = version_list

    def write(self, packed_models):
        """
        Create and write added annotations to individual EMICSS file for every EMDB entry
        """
        emdb_id = self.emdb_id
        headerXML = emicss(emdb_id=emdb_id)
        dbs = dbsType(collection_date=self.version_list['date'])
        entry_ref_dbs = entry_ref_dbsType()
        weights = weightsType()
        sample = sampleType()
        macromolecules = macromoleculesType()
        supramolecules = supramoleculesType()
        primary_citation = primary_citationType()
        all_db = set()

        if "EMPIAR" in packed_models:
            empiar_objects = packed_models['EMPIAR']
            if len(empiar_objects) > 0:
                all_db.add("EMPIAR")
                for empiar in empiar_objects:
                    emp_ref_db_obj = entry_ref_dbType(source="EMPIAR", accession_id=empiar.empiar_id, provenance="EMDB")
                    entry_ref_dbs.add_entry_ref_db(emp_ref_db_obj)
        # MW calculated from header
        if "WEIGHT" in packed_models: 
            mweight = packed_models["WEIGHT"]
            if len(mweight) > 0:
                all_db.add("EMDB")
                for mw in mweight:
                    mw_info_obj = weight_infoType(weight=round(mw.overall_mw,2), unit=mw.units, provenance=mw.provenance)
                    weights.add_weight_info(mw_info_obj)
        # MW calculated from assemblies
        if "MODEL" in packed_models:
            pdb_objects = packed_models["MODEL"]
            if len(pdb_objects) > 0:
                all_db.add("EMDB")
                for pdb in pdb_objects:
                    pdb_info_obj = weight_infoType(pdb_id=pdb.pdb_id, assemblies=pdb.assembly, weight=pdb.molecular_weight, unit="Da", provenance="PDBe")
                    weights.add_weight_info(pdb_info_obj)
        if "CITATION" in packed_models:
            citations = packed_models["CITATION"]
            if len(citations) > 0:
                for citation in citations:
                    authors_obj = authorsType()
                    for author in citation.authors:
                        orcid = author.orcid if author.orcid else None
                        author_obj = authorType(name=author.name, orcid_id=orcid, order=author.order, provenance=author.provenance)
                        authors_obj.add_author(author_obj)
                    primary_citation.set_authors(authors_obj)
                    if citation.pmedid:
                        all_db.add("PubMed")
                        pm_citation_obj = ref_citationType(source="PubMed", accession_id=citation.pmedid, provenance=citation.provenance_pm)
                        primary_citation.add_ref_citation(pm_citation_obj)
                        if citation.pmcid:
                            all_db.add("PubMed Central")
                            pmc_citation_obj = ref_citationType(source="PubMed Central", accession_id=citation.pmcid, provenance=citation.provenance_pmc)
                            primary_citation.add_ref_citation(pmc_citation_obj)
                        if citation.issn:
                            all_db.add("ISSN")
                            issn_citation_obj = ref_citationType(source="ISSN", accession_id=citation.issn, provenance=citation.provenance_issn)
                            primary_citation.add_ref_citation(issn_citation_obj)
                        if citation.doi:
                            primary_citation.set_doi(citation.doi)
                            primary_citation.set_provenance(citation.provenance_doi)
        if "PROTEIN-TERMS" in packed_models:
            proteins = packed_models["PROTEIN-TERMS"]
            if len(proteins) > 0:
                all_db.add("EMDB")
                for protein in proteins:
                    cross_ref_dbs = cross_ref_dbsType()
                    if len(protein.uniprot_id) > 0:
                        all_db.add("UniProt")
                        for uni in protein.uniprot_id:
                            unp_xref_obj = cross_ref_db(source="UniProt", accession_id=uni.uniprot_id, provenance=uni.provenance)
                            cross_ref_dbs.add_cross_ref_db(unp_xref_obj)
                    if len(protein.go) > 0:
                        all_db.add("GO")
                        for go in protein.go:
                            go_xref_obj = cross_ref_db(name=go.namespace, source="GO", accession_id=go.id, provenance=go.provenance)
                            if go.type == "P": go_xref_obj.set_type("biological process")
                            elif go.type == "C": go_xref_obj.set_type("cellular component")
                            elif go.type == "F": go_xref_obj.set_type("molecular function")
                            cross_ref_dbs.add_cross_ref_db(go_xref_obj)
                    if len(protein.interpro) > 0:
                        all_db.add("InterPro")
                        for ipr in protein.interpro:
                            ipr_xref_obj = cross_ref_db(name=ipr.namespace, source="InterPro", accession_id=ipr.id, uniprot_start=ipr.start, uniprot_end=ipr.end, provenance=ipr.provenance)
                            cross_ref_dbs.add_cross_ref_db(ipr_xref_obj)
                    if len(protein.pfam) > 0:
                        all_db.add("Pfam")
                        for pfam in protein.pfam:
                            pfam_xref_obj = cross_ref_db(name=pfam.namespace, source="Pfam", accession_id=pfam.id, uniprot_start=pfam.start, uniprot_end=pfam.end, provenance=pfam.provenance)
                            cross_ref_dbs.add_cross_ref_db(pfam_xref_obj)
                    if len(protein.cath) > 0:
                        all_db.add("CATH")
                        for cath in protein.cath:
                            cath_xref_obj = cross_ref_db(source="CATH", accession_id=cath.id, uniprot_start=cath.start, uniprot_end=cath.end, provenance=cath.provenance)
                            cross_ref_dbs.add_cross_ref_db(cath_xref_obj)
                    if len(protein.scop) > 0:
                        all_db.add("SCOP")
                        for scop in protein.scop:
                            scop_xref_obj = cross_ref_db(source="SCOP", accession_id=scop.id, uniprot_start=scop.start, uniprot_end=scop.end, provenance=scop.provenance)
                            cross_ref_dbs.add_cross_ref_db(scop_xref_obj)
                    if len(protein.scop2) > 0:
                        all_db.add("SCOP2")
                        for scop2 in protein.scop2:
                            scop2_xref_obj = cross_ref_db(source="SCOP2", accession_id=scop2.id, uniprot_start=scop2.start, uniprot_end=scop2.end, provenance=scop2.provenance)
                            cross_ref_dbs.add_cross_ref_db(scop2_xref_obj)
                    if len(protein.scop2B) > 0:
                        all_db.add("SCOP2")
                        for scop2B in protein.scop2B:
                            scop2B_xref_obj = cross_ref_db(source="SCOP2B", accession_id=scop2B.id, uniprot_start=scop2B.start, uniprot_end=scop2B.end, provenance=scop2B.provenance)
                            cross_ref_dbs.add_cross_ref_db(scop2B_xref_obj)
                    if len(protein.pdbekb) > 0:
                        all_db.add("PDBe-KB")
                        for pdbekb in protein.pdbekb:
                            pdbekb_xref_obj = cross_ref_db(source="PDBe-KB", accession_id=pdbekb.unip_id, provenance=pdbekb.provenance)
                            cross_ref_dbs.add_cross_ref_db(pdbekb_xref_obj)
                    if len(protein.alphafold) > 0:
                        all_db.add("AlphaFold DB")
                        for alphafold in protein.alphafold:
                            afdb_xref_obj = cross_ref_db(source="AlphaFold DB", accession_id=alphafold.unip_id, provenance=alphafold.provenance)
                            cross_ref_dbs.add_cross_ref_db(afdb_xref_obj)
                    if cross_ref_dbs.hasContent_():
                        macromolecule = macromoleculeType(type_="protein", id=int(protein.sample_id), copies=int(protein.sample_copies), provenance="EMDB", name=protein.sample_name, cross_ref_dbs=cross_ref_dbs)
                        macromolecules.add_macromolecule(macromolecule)
        if "LIGANDS" in packed_models:
            ligand_obj = packed_models["LIGANDS"]
            if len(ligand_obj) > 0:
                for ligand in ligand_obj:
                    cross_ref_dbs = cross_ref_dbsType()
                    if ligand.chembl_id:
                        all_db.add("ChEMBL")
                        chembl_obj = cross_ref_db(source="ChEMBL", accession_id=ligand.chembl_id, provenance=ligand.provenance_chembl)
                        cross_ref_dbs.add_cross_ref_db(chembl_obj)
                    if ligand.chebi_id:
                        all_db.add("ChEBI")
                        chebi_obj = cross_ref_db(source="ChEBI", accession_id=ligand.chebi_id, provenance=ligand.provenance_chebi)
                        cross_ref_dbs.add_cross_ref_db(chebi_obj)
                    if ligand.drugbank_id:
                        all_db.add("DrugBank")
                        drugbank_obj = cross_ref_db(source="DrugBank", accession_id=ligand.drugbank_id, provenance=ligand.provenance_drugbank)
                        cross_ref_dbs.add_cross_ref_db(drugbank_obj)
                    if cross_ref_dbs.hasContent_():
                        macromolecule = macromoleculeType(type_="ligand", id=ligand.sample_id, copies=ligand.lig_copies,
                            provenance="EMDB", name=ligand.lig_name, ccd_id=ligand.HET, cross_ref_dbs=cross_ref_dbs)
                        macromolecules.add_macromolecule(macromolecule)
        if "COMPLEX" in packed_models:
            complex_objects = packed_models['COMPLEX']
            if len(complex_objects) > 0:
                all_db.add("Complex Portal")
                for emdb_complex in complex_objects:
                    if emdb_complex:
                        cross_ref_dbs = cross_ref_dbsType()
                        if emdb_complex.cpx_list:
                            for cpx in emdb_complex.cpx_list:
                                cpx_obj = cross_ref_db(name=cpx[1], source="Complex Portal", accession_id=cpx[0], provenance=emdb_complex.provenance, score=round(emdb_complex.score,2))
                                cross_ref_dbs.add_cross_ref_db(cpx_obj)
                            if cross_ref_dbs.hasContent_():
                                sample_id = emdb_complex.sample_id.split('_')[1]
                                supramolecule = supramoleculeType(type_="complex", id=sample_id, copies=emdb_complex.sample_copies,
                                    provenance=emdb_complex.provenance, name=emdb_complex.supra_name, cross_ref_dbs=cross_ref_dbs)
                                supramolecules.add_supramolecule(supramolecule)

        for database in all_db:
            if database in self.version_list:
                version = self.version_list[database]
                db_ver = dbType(source=database, version=version)
                dbs.add_db(db_ver)

        if all_db:
            headerXML.set_dbs(dbs)
        if entry_ref_dbs.hasContent_():
            headerXML.set_entry_ref_dbs(entry_ref_dbs)
        if primary_citation.hasContent_():
            headerXML.set_primary_citation(primary_citation)
        if weights.hasContent_():
            headerXML.set_weights(weights)
        if supramolecules.hasContent_():
            sample.set_supramolecules(supramolecules)
        if macromolecules.hasContent_():
            sample.set_macromolecules(macromolecules)
        if sample.hasContent_():
            headerXML.set_sample(sample)

        output_path = os.path.join(self.workDir, "emicss")
        Path(output_path).mkdir(parents=True, exist_ok=True)
        xmlFile = os.path.join(output_path, f"emd_{emdb_id[4:]}_emicss.xml")
        with open(xmlFile, 'w') as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            headerXML.export(f, 0, name_='emicss', namespacedef_=f'version="{headerXML.version}" '
                                                                 f'schema_location="https://github.com/emdb-empiar/emicss-schema/tree/main/versions/emdb_emicss_{headerXML.version}.xsd" ')
