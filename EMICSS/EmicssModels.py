import os
import models

class EmicssModels:
    """
     Convert the annotations to dictionary from log_files to write the XML files
    """

    def __init__(self, workDir):
        self.workDir = workDir
        self.packed_models = {}
        self.PT_dict = {}

    def worker(self):
        #### EMPIAR
        fileemp = os.path.join(self.workDir, "emdb_empiar.log")
        with open(fileemp, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                row = line.strip('\n').split('\t')
                emp = models.Empiar(emdb_id=row[0], empiar_id=row[1])
                self.generating_dictionary(row[0], "EMPIAR", emp)
        #### WEIGHT
        filemw = os.path.join(self.workDir, "emdb_overall_mw.log")
        with open(filemw, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                row = line.strip('\n').split('\t')
                wgt = models.Weight(emdb_id=row[0], overall_mw=float(row[1]), units="MDa", provenance="EMDB")
                self.generating_dictionary(row[0], "WEIGHT", wgt)
        #### MODEL
        filemod = os.path.join(self.workDir, "emdb_model.log")
        with open(filemod, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                row = line.strip('\n').split('\t')
                model = models.Model(emdb_id=row[0], pdb_id=row[1], assembly=row[2], molecular_weight=row[3])
                self.generating_dictionary(row[0], "MODEL", model)
        #### CITATION
        filepub = os.path.join(self.workDir, "emdb_pubmed.log")
        author_dict = {}
        fileauth = os.path.join(self.workDir, 'emdb_author.log')
        with open(fileauth, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                authors = models.Author(name=rw[1], orcid=rw[2], order=rw[3], provenance=rw[4])
                if rw[0] not in author_dict:
                    author_dict[rw[0]] = [authors]
                else:
                    author_dict[rw[0]].append(authors)
        with open(filepub, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                row = line.strip('\n').split('\t')
                authors_list = author_dict[row[0]]
                pubm = models.Citation(emdb_id=row[0], pmedid=row[1], pmcid=row[3], issn=row[5], doi=row[7], authors=authors_list,
                                       provenance_pm=row[2], provenance_pmc=row[4], provenance_issn=row[6], provenance_doi=row[8])
                self.generating_dictionary(row[0], "CITATION", pubm)
        #### LIGANDS
        het_dict = {}
        filechembl = os.path.join(self.workDir, 'emdb_chembl.log')
        with open(filechembl, 'r') as fileconts:
            for ln in fileconts.readlines()[1:]:
                rw = ln.strip('\n').split('\t')
                hetchembl = {"sample_id": rw[1], "lig_name": rw[3], "lig_copies": rw[4], "chembl_id": rw[5], "provenance_chembl": rw[6]}
                if rw[0] not in het_dict:
                    het_dict[rw[0]] = {}
                if rw[2] not in het_dict[rw[0]]:
                    het_dict[rw[0]][rw[2]] = hetchembl
        filechebi = os.path.join(self.workDir, 'emdb_chebi.log')
        with open(filechebi, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                row = line.strip('\n').split('\t')
                hetchebi = {"sample_id": row[1], "lig_name": row[3], "lig_copies": row[4], "chebi_id": row[5], "provenance_chebi": row[6]}
                if row[0] not in het_dict:
                    het_dict[row[0]] = {}
                if row[2] not in het_dict[row[0]]:
                    het_dict[row[0]][row[2]] = hetchebi
                else:
                    het_dict[row[0]][row[2]].update(hetchebi)
        filedb = os.path.join(self.workDir, 'emdb_drugbank.log')
        with open(filedb, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                rows = line.strip('\n').split('\t')
                hetdb = {"sample_id": rows[1], "lig_name": rows[3], "lig_copies": rows[4], "drugbank_id": rows[5], "provenance_drugbank": rows[6]}
                if rows[0] not in het_dict:
                    het_dict[rows[0]] = {}
                if rows[2] not in het_dict[rows[0]]:
                    het_dict[rows[0]][rows[2]] = hetdb
                else:
                    het_dict[rows[0]][rows[2]].update(hetdb)
        for em_key in list(het_dict.keys()):
            for het in list(het_dict[em_key].keys()):
                allkeys = (["emdb_id", "sample_id", "lig_name", "lig_copies", "chembl_id", "provenance_chembl", "chebi_id",
                           "provenance_chebi", "drugbank_id", "provenance_drugbank"])
                het_values = het_dict[em_key][het]
                for index in allkeys:
                    if index not in het_values.keys():
                        het_values[index] = ''
                lig_db = models.Ligand(emdb_id=het_values['emdb_id'], sample_id=het_values['sample_id'], HET=het, lig_name=het_values['lig_name'],
                                       lig_copies=het_values['lig_copies'], chembl_id=het_values['chembl_id'], provenance_chembl=het_values['provenance_chembl'],
                                       chebi_id=het_values['chebi_id'], provenance_chebi=het_values['provenance_chebi'],
                                       drugbank_id=het_values['drugbank_id'], provenance_drugbank=het_values['provenance_drugbank'])
                self.generating_dictionary(em_key, "LIGANDS", lig_db)
        # #### COMPLEX
        filecpx = os.path.join(self.workDir, "emdb_cpx.log")
        cpx_dict = {}
        with open(filecpx, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                row = line.strip('\n').split('\t')
                cpxl = [row[4], row[5]]
                cpx_dict = self.generating_PT_dictionary(row[0], row[1], "cpx_list", cpxl)
                if "supra_name" not in cpx_dict[row[0]][row[1]]:
                    cpx_dict[row[0]][row[1]]["supra_name"] = row[2]
                if "supra_copies" not in cpx_dict[row[0]][row[1]]:
                    cpx_dict[row[0]][row[1]]["supra_copies"] = row[3]
                if "score" not in cpx_dict[row[0]][row[1]]:
                    cpx_dict[row[0]][row[1]]["score"] = row[7]
                if "provenance" not in cpx_dict[row[0]][row[1]]:
                    cpx_dict[row[0]][row[1]]["provenance"] = row[6]
        for em_key in list(cpx_dict.keys()):
            for samp_key in list(cpx_dict[em_key].keys()):
                cpx_val = cpx_dict[em_key][samp_key]
                cpx = models.EMDB_complex(emdb_id=em_key, sample_id=samp_key, supra_name=cpx_val["supra_name"], sample_copies=cpx_val["supra_copies"],
                                          complex_sample_id="", cpx_list=cpx_val["cpx_list"], provenance=cpx_val["provenance"], score=float(cpx_val["score"]))
                self.generating_dictionary(em_key, "COMPLEX", cpx)
        #### PROTEIN_TERMS
        self.PT_dict = {}
        fileuni = os.path.join(self.workDir, 'emdb_uniprot.log')
        with open(fileuni, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                uni = models.Protein(emdb_id=rw[0], sample_id=rw[1], sample_name=rw[2], sample_copies=rw[3], uniprot_id=rw[5], provenance=rw[6])
                uni_dict = self.generating_PT_dictionary(rw[0], rw[1], "uni", uni)
                if "sample_name" not in uni_dict[rw[0]][rw[1]]:
                    uni_dict[rw[0]][rw[1]]["sample_name"] = rw[2]
                if "sample_copies" not in uni_dict[rw[0]][rw[1]]:
                    uni_dict[rw[0]][rw[1]]["sample_copies"] = rw[3]
        filego = os.path.join(self.workDir, 'emdb_go.log')
        with open(filego, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                go = models.GO(id=rw[2], namespace=rw[3], type=rw[4], provenance=rw[5])
                self.generating_PT_dictionary(rw[0], rw[1], "go", go)
        fileinter = os.path.join(self.workDir, 'emdb_interpro.log')
        with open(fileinter, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                interpro = models.Interpro(id=rw[2], namespace=rw[3], start=rw[4], end=rw[5], unp_start=rw[6], unp_end=rw[7], provenance=rw[8])
                self.generating_PT_dictionary(rw[0], rw[1], "inter", interpro)
        filepfam = os.path.join(self.workDir, 'emdb_pfam.log')
        with open(filepfam, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                pfam = models.Pfam(id=rw[2], namespace=rw[3], start=rw[4], end=rw[5], unp_start=rw[6], unp_end=rw[7], provenance=rw[8])
                self.generating_PT_dictionary(rw[0], rw[1], "pfam", pfam)
        filecath = os.path.join(self.workDir, 'emdb_cath.log')
        with open(filecath, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                cath = models.Cath(id=rw[2], start=rw[3], end=rw[4], unp_start=rw[5], unp_end=rw[6], provenance=rw[7])
                self.generating_PT_dictionary(rw[0], rw[1], "cath", cath)
        filescop = os.path.join(self.workDir, 'emdb_scop.log')
        with open(filescop, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                scop = models.SCOP(id=rw[2], start=rw[3], end=rw[4], unp_start=rw[5], unp_end=rw[6], provenance=rw[7])
                self.generating_PT_dictionary(rw[0], rw[1], "scop", scop)
        filescop2 = os.path.join(self.workDir, 'emdb_scop2.log')
        with open(filescop2, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                scop2 = models.SCOP2(id=rw[2], start=rw[3], end=rw[4], unp_start=rw[5], unp_end=rw[6], provenance=rw[7])
                self.generating_PT_dictionary(rw[0], rw[1], "scop2", scop2)
        filescop2B = os.path.join(self.workDir, 'emdb_scop2B.log')
        with open(filescop2B, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                scop2B = models.SCOP(id=rw[2], start=rw[3], end=rw[4], unp_start=rw[5], unp_end=rw[6], provenance=rw[7])
                self.generating_PT_dictionary(rw[0], rw[1], "scop2B", scop2B)
        filekb = os.path.join(self.workDir, 'emdb_pdbekb.log')
        with open(filekb, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                pdbekb = models.Pdbekb(unip_id=rw[2], provenance=rw[3])
                self.generating_PT_dictionary(rw[0], rw[1], "kb", pdbekb)
        fileaf = os.path.join(self.workDir, 'emdb_alphafold.log')
        with open(fileaf, 'r') as filecontents:
            for lin in filecontents.readlines()[1:]:
                rw = lin.strip('\n').split('\t')
                af = models.Alphafold(unip_id=rw[2], provenance=rw[3])
                self.generating_PT_dictionary(rw[0], rw[1], "af", af)
        for PT in list(self.PT_dict.keys()):
            PT_samp = self.PT_dict[PT]
            for PTs in list(PT_samp.keys()):
                PTres = self.PT_dict[PT][PTs]
                allkeys = (["uni", "go", "inter", "pfam", "cath", "scop", "scop2", "scop2B", "kb", "af"])
                for index in allkeys:
                    if index not in PTres.keys():
                        PTres[index] = []
                PT_db = models.Protein(emdb_id=PT, sample_id=PTs, sample_copies=PTres["sample_copies"], sample_name=PTres["sample_name"],
                                       uniprot_id=PTres["uni"], go=PTres["go"], interpro=PTres["inter"], pfam=PTres["pfam"], cath=PTres["cath"],
                                       scop=PTres["scop"], scop2=PTres["scop2"], scop2B=PTres["scop2B"], pdbekb=PTres["kb"],
                                       alphafold=PTres["af"], provenance=PTres["uni"])
                self.generating_dictionary(PT, "PROTEIN-TERMS", PT_db)
        return self.packed_models

    def generating_dictionary(self, emdb_id, db, emicss_values):
        if emdb_id not in self.packed_models:
            self.packed_models[emdb_id] = {}
        if db not in self.packed_models[emdb_id]:
            self.packed_models[emdb_id][db] = [emicss_values]
        else:
            self.packed_models[emdb_id][db].append(emicss_values)
        return self.packed_models

    def generating_PT_dictionary(self, em_id, samp_id, db, db_values):
        if em_id not in self.PT_dict:
            self.PT_dict[em_id] = {}
        if samp_id not in self.PT_dict[em_id]:
            self.PT_dict[em_id][samp_id] = {}
        if db not in self.PT_dict[em_id][samp_id]:
            self.PT_dict[em_id][samp_id][db] = [db_values]
        else:
            self.PT_dict[em_id][samp_id][db].append(db_values)
        return self.PT_dict