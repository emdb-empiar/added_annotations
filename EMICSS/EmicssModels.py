import os
import models

class EmicssModels:
    """
     Convert the annotations to dictionary from log_files to write the XML files
    """

    def __init__(self, workDir):
        self.workDir = workDir
        self.packed_models = {}

    def worker(self):
        #### EMPIAR
        fileemp = os.path.join(self.workDir, "emdb_empiar.log")
        with open(fileemp, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                row = line.strip('\n').split('\t')
                emp = models.Empiar(emdb_id=row[0], empiar_id=row[1])
                self.generating_dictionary(row, "EMPIAR", emp)
        #### WEIGHT
        filemw = os.path.join(self.workDir, "emdb_overall_mw.log")
        with open(filemw, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                row = line.strip('\n').split('\t')
                wgt = models.Weight(emdb_id=row[0], overall_mw=float(row[1]), units="MDa", provenance="EMDB")
                self.generating_dictionary(row, "WEIGHT", wgt)
        #### MODEL
        filemod = os.path.join(self.workDir, "emdb_model.log")
        with open(filemod, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                row = line.strip('\n').split('\t')
                model = models.Model(emdb_id=row[0], pdb_id=row[1], assembly=row[2], molecular_weight=row[3])
                self.generating_dictionary(row, "MODEL", model)
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
                self.generating_dictionary(row, "CITATION", pubm)
        #### LIGANDS
        het_dict = {}
        filechembl = os.path.join(self.workDir, 'emdb_chembl.log')
        with open(filechembl, 'r') as fileconts:
            for ln in fileconts.readlines()[1:]:
                rw = ln.strip('\n').split('\t')
                hetchembl = {"emdb_id": rw[0], "sample_id": rw[1], "lig_name": rw[3], "lig_copies": rw[4], "chembl_id": rw[5], "provenance_chembl": rw[6]}
                             # "chebi_id": "", "provenance_chebi": "", "drugbank_id": "", "provenance_drugbank": ""}
                if rw[2] not in het_dict:
                    het_dict[rw[2]] = hetchembl
        filechebi = os.path.join(self.workDir, 'emdb_chebi.log')
        with open(filechebi, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                row = line.strip('\n').split('\t')
                hetchebi = {"emdb_id": row[0], "sample_id": row[1], "lig_name": row[3], "lig_copies": row[4], "chebi_id": row[5], "provenance_chebi": row[6]}
                if row[2] not in het_dict:
                    het_dict[row[2]] = hetchebi
                else:
                    het_dict[row[2]].update(hetchebi)
        filedb = os.path.join(self.workDir, 'emdb_drugbank.log')
        with open(filedb, 'r') as filecont:
            for line in filecont.readlines()[1:]:
                rows = line.strip('\n').split('\t')
                hetdb = {"emdb_id": rows[0], "sample_id": rows[1], "lig_name": rows[3], "lig_copies": rows[4], "drugbank_id": rows[5], "provenance_drugbank": rows[6]}
                if rows[2] not in het_dict:
                    het_dict[rows[2]] = hetdb
                else:
                    het_dict[rows[2]].update(hetdb)
        for het in list(het_dict.keys()):
            allkeys = (["emdb_id", "sample_id", "lig_name", "lig_copies", "chembl_id", "provenance_chembl", "chebi_id",
                       "provenance_chebi", "drugbank_id", "provenance_drugbank"])
            het_values = het_dict[het]
            for index in allkeys:
                if index not in het_values.keys():
                    het_values[index] = ''
            lig_db = models.Ligand(emdb_id=het_values['emdb_id'], sample_id=het_values['sample_id'], HET=het, lig_name=het_values['lig_name'],
                                   lig_copies=het_values['lig_copies'], chembl_id=het_values['chembl_id'], provenance_chembl=het_values['provenance_chembl'],
                                   chebi_id=het_values['chebi_id'], provenance_chebi=het_values['provenance_chebi'],
                                   drugbank_id=het_values['drugbank_id'], provenance_drugbank=het_values['provenance_drugbank'])
            self.generating_dictionary(row, "LIGANDS", lig_db)

        return self.packed_models

    def generating_dictionary(self, row, db, emicss_values):
        if row[0] not in self.packed_models:
            self.packed_models[row[0]] = {}
        if db not in self.packed_models[row[0]]:
            self.packed_models[row[0]][db] = [emicss_values]
        else:
            self.packed_models[row[0]][db].append(emicss_values)
        return self.packed_models
