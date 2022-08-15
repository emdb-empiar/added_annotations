from EMICSS.EMICSS import *
from models import *
import os

class Parser:
    """
    Read log files and create emicss objects
    """
    def __init__(self, version):
        self.emdb_ids = set()
        self.proteins = {}
        self.empiars = {}
        self.weights = {}
        self.models = {}
        self.citations = {}
        self.ligands = {}
        self.complexes = {}
        self.version = version
        self.__parse_uniprot()
        self.__parse_empiar()
        self.__parse_mw()
        self.__parse_models()
        self.__parse_citations()
        self.__parse_chembl()
        self.__parse_chebi()
        self.__parse_drugbank()
        self.__parse_complex()
        self.__parse_go()
        self.__parse_interpro()
        self.__parse_pfam()
        self.__parse_cath()
        self.__parse_scop()
        self.__parse_scop2()
        self.__parse_scop2B()
        self.__parse_pdbekb()
        self.__parse_afdb()

    def __parse_uniprot(self):
        fileuniprot = os.path.join(self.workDir, 'emdb_uniprot.log')
        with open(fileuniprot, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, sample_name, copies, ncbi_id, uniprot_id, provenance, sample_complex_ids = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                protein = Protein(emdb_id=emdb_id, sample_id=sample_id, sample_name=sample_name, sample_copies=copies, sample_organism=ncbi_id, uniprot_id=uniprot_id, provenance=provenance)
                self.entries[emdb_sample_id] = protein

    def __parse_empiar(self):
        fileempiar = os.path.join(self.workDir, "emdb_empiar.log")
        with open(fileempiar, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, empiar_id, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                empiar = Empiar(emdb_id=emdb_id, empiar_id=empiar_id)
                self.empiars[emdb_id] = empiar

    def __parse_mw(self):
        filemw = os.path.join(self.workDir, "overall_mw.log")
        with open(filemw, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, mw = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                weight = Weight(emdb_id=emdb_id, overall_mw=float(mw), units="MDa", provenance="EMDB")
                self.weights[emdb_id] = weight

    def __parse_models(self):
        filemodel = os.path.join(self.workDir, "emdb_model.log")
        with open(filemodel, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, pdb_id, assembly, mw = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                model = Model(emdb_id=emdb_id, pdb_id=pdb_id, assembly=assembly, molecular_weight=float(mw))
                self.models[emdb_id] = model

    def __parse_citations(self):
        author_dict = {}
        fileauthor = os.path.join(self.workDir, 'emdb_author.log')
        with open(fileauthor, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, author_name, orcid_id, order, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                author = Author(name=author_name, orcid=orcid_id, order=order, provenance=provenance)
                if emdb_id in author_dict:
                    author_dict[emdb_id].append(author)
                else:
                    author_dict[emdb_id] = [author]
        filepubmed = os.path.join(self.workDir, "emdb_pubmed.log")
        with open(filepubmed, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, pubmed_id, pubmed_provenance, pmc_id, pmc_provenance, issn, issn_provenance, doi, doi_provenance, journal_name, journal_abbv = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                authors_list = author_dict[emdb_id]
                citation = Citation(emdb_id=emdb_id, pmedid=pubmed_id, provenance_pm=pubmed_provenance, pmcid=pmc_id, 
                                    provenance_pmc=pmc_provenance, issn=issn, provenance_issn=issn_provenance, doi=doi, 
                                    provenance_doi=doi_provenance, authors=authors_list)
                self.citations[emdb_id] = citation

    def __parse_chembl(self):
        filechembl = os.path.join(self.workDir, 'emdb_chembl.log')
        with open(filechembl, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, het_code, ligand_name, copies, chembl_id, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                if emdb_sample_id in self.ligands:
                    self.ligands[emdb_sample_id].chembl_id = chembl_id
                    self.ligands[emdb_sample_id].provenance_chembl = provenance
                else:
                    ligand = Ligand(emdb_id=emdb_id, sample_id=sample_id, chembl_id=chembl_id, 
                        provenance_chembl=provenance, HET=het_code, lig_name=ligand_name, lig_copies=copies)
                    self.ligands[emdb_sample_id] = ligand

    def __parse_chebi(self):
        filechebi = os.path.join(self.workDir, 'emdb_chebi.log')
        with open(filechebi, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, het_code, ligand_name, copies, chebi_id, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                if emdb_sample_id in self.ligands:
                    self.ligands[emdb_sample_id].chebi_id = chebi_id
                    self.ligands[emdb_sample_id].provenance_chebi = provenance
                else:
                    ligand = Ligand(emdb_id=emdb_id, sample_id=sample_id, chebi_id=chebi_id,
                        provenance_chebi=provenance, HET=het_code, lig_name=ligand_name, lig_copies=copies)
                    self.ligands[emdb_sample_id] = ligand

    def __parse_drugbank(self):
        filedb = os.path.join(self.workDir, 'emdb_drugbank.log')
        with open(filedb, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, het_code, ligand_name, copies, db_id, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                if emdb_sample_id in self.ligands:
                    self.ligands[emdb_sample_id].drugbank_id = db_id
                    self.ligands[emdb_sample_id].provenance_drugbank = provenance
                else:
                    ligand = Ligand(emdb_id=emdb_id, sample_id=sample_id, drugbank_id=db_id,
                        provenance_drugbank=provenance, HET=het_code, lig_name=ligand_name, lig_copies=copies)
                    self.ligands[emdb_sample_id] = ligand

    def __parse_complex(self):
        filecpx = os.path.join(self.workDir, "emdb_cpx.log")
        cpx_dict = {}
        with open(filecpx, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, sample_name, copies, cpx_id, cpx_title, provenance, score = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                cpx = CPX([cpx_id, cpx_title, "", "", "", "", "", "", ""])
                if emdb_sample_id in self.complexes:
                    self.complexes[emdb_sample_id].cpx_list.append(cpx)
                else:
                    self.complexes[emdb_sample_id] = EMDB_complex(emdb_id=emdb_id, sample_id=emdb_sample_id, supra_name=sample_name,
                                                                  sample_copies=copies, complex_sample_id=sample_id, cpx_list=[cpx],
                                                                  proteins=None,provenance=provenance, score=float(score))

    def __parse_go(self):
        filego = os.path.join(self.workDir, 'emdb_go.log')
        with open(filego, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, go_id, go_namespace, go_type, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                go = GO(id=go_id, namespace=go_namespace, type=go_type, provenance=provenance)
                self.proteins[emdb_sample_id].go.add(go)

    def __parse_interpro(self):
        fileinterpro = os.path.join(self.workDir, 'emdb_interpro.log')
        with open(fileinterpro, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, interpro_id, interpro_namespace, start, end, uniprot_start, uniprot_end, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                interpro = Interpro(id=interpro_id, namespace=interpro_namespace, start=int(start), end=int(end),
                                    unp_start=int(uniprot_start), unp_end=int(uniprot_end), provenance=provenance)
                self.proteins[emdb_sample_id].interpro.add(interpro)

    def __parse_pfam(self):
        filepfam = os.path.join(self.workDir, 'emdb_pfam.log')
        with open(filepfam, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, pfam_id, pfam_namespace, start, end, uniprot_start, uniprot_end, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                pfam = Pfam(id=pfam_id, namespace=pfam_namespace, start=int(start), end=int(end), unp_start=int(uniprot_start),
                            unp_end=int(uniprot_end), provenance=provenance)
                self.proteins[emdb_sample_id].pfam.add(pfam)

    def __parse_cath(self):
        filecath = os.path.join(self.workDir, 'emdb_cath.log')
        with open(filecath, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, cath_id, start, end, unipro_start, uniprot_end, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                cath = Cath(id=cath_id, start=int(start), end=int(end), unp_start=int(uniprot_start), unp_end=int(uniprot_end),
                            provenance=provenance)
                self.proteins[emdb_sample_id].cath.add(cath)

    def __parse_scop(self):
        filescop = os.path.join(self.workDir, 'emdb_scop.log')
        with open(filescop, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, scop_id, start, end, uniprot_start, uniprot_end, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                scop = SCOP(id=scop_id, start=int(start), end=int(end), unp_start=int(uniprot_start), unp_end=int(uniprot_end),
                            provenance=provenance)
                self.proteins[emdb_sample_id].scop.add(scop)

    def __parse_scop2(self):
        filescop = os.path.join(self.workDir, 'emdb_scop2.log')
        with open(filescop, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, scop_id, start, end, uniprot_start, uniprot_end, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                scop2 = SCOP2(id=scop_id, start=int(start), end=int(end), unp_start=int(uniprot_start), unp_end=int(uniprot_end),
                            provenance=provenance)
                self.proteins[emdb_sample_id].scop2.add(scop2)

    def __parse_scop2B(self):
        filescop = os.path.join(self.workDir, 'emdb_scop2B.log')
        with open(filescop, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, scop_id, start, end, uniprot_start, uniprot_end, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                scop2b = SCOP2B(id=scop_id, start=int(start), end=int(end), unp_start=int(uniprot_start), unp_end=int(uniprot_end),
                            provenance=provenance)
                self.proteins[emdb_sample_id].scop2b.add(scop2b)

    def __parse_pdbekb(self):
        filepdbekb = os.path.join(self.workDir, 'emdb_pdbekb.log')
        with open(filepdbekb, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, pdbekb_id, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                pdbekb = Pdbekb(uniprot_id=pdbekb_id, provenance=provenance)
                self.proteins[emdb_sample_id].pdbekb = pdbekb

    def __parse_afdb(self):
        fileafdb = os.path.join(self.workDir, 'emdb_alphafold.log')
        with open(fileafdb, 'r') as filereader:
            for line in filereader.readlines()[1:]:
                emdb_id, sample_id, afdb_id, provenance = line.strip().split('\t')
                self.emdb_ids.add(emdb_id)
                emdb_sample_id = f"{emdb_id}_{sample_id}"
                afdb = Alphafold(uniprot_id=afdb_id, provenance=provenance)
                self.proteins[emdb_sample_id].alphafold = afdb



class EmicssXML:
    """
    Store and Write EMICSS data
    """
    def __init__(self, emdb_id, version):
        self.emdb_id = emdb_id
        self.headerXML = emicss(emdb_id=emdb_id)
        self.dbs = dbsType(collection_date=version.today)
        self.entry_ref_dbs = entry_ref_dbsType()
        self.weights = weightsType()
        self.sample = sampleType()
        self.macromolecules = macromoleculesType()
        self.supramolecules = supramoleculesType()
        self.primary_citation = primary_citationType()
        self.all_db = set()




