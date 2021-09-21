import requests
from models import Pdbekb, Alphafold

uni_api = r'https://www.uniprot.org/uniprot/?query=id:'
uni_format = r'%20database:(type:pdb)&sort=score&columns=id&format=tab'
alphafold_api = r'https://alphafold.ebi.ac.uk/api/prediction/'

class PdbeKB_AlphafoldMapping:
    """
    Validates the PDBe-KB and Alphafold links for the UniProt ID
    """

    def __init__(self, proteins, is_pdbekb=True, is_alphafold=True):
        self.proteins = proteins
        self.is_pdbekb = is_pdbekb
        self.is_alphafold = is_alphafold

    def execute(self):
        for protein in self.proteins:
            if protein.uniprot_id:
                protein = self.worker(protein)
        return self.proteins

    def worker(self, protein):
        return self.link_validate(protein)

    def link_validate(self, protein):
        uid = protein.uniprot_id
        pdbekb_url = uni_api + uid + uni_format
        response = requests.get(pdbekb_url)
        if response.status_code == 200 and response.content:
            if self.is_pdbekb:
                pdbekb = Pdbekb()
                pdbekb.emdb_id = protein.emdb_id
                pdbekb.sample_id = protein.sample_id
                pdbekb.link = "https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/" + uid
                pdbekb.provenance = "PDBe-KB"
                protein.pdbekb.append(pdbekb)

        alphafold_url = alphafold_api + uid
        response = requests.get(alphafold_url)
        if response.status_code == 200 and response.content:
            if self.is_alphafold:
                alphafold = Alphafold()
                alphafold.emdb_id = protein.emdb_id
                alphafold.sample_id = protein.sample_id
                alphafold.link = "https://alphafold.ebi.ac.uk/entry/" + uid
                alphafold.provenance = "ALPHAFOLD"
                protein.alphafold.append(alphafold)
        return protein

    def export_tsv(self, pdbekb_logger, alphafold_logger):
        for protein in self.proteins:
            if self.is_pdbekb and protein.pdbekb:
                for pdbekb in protein.pdbekb:
                    row = f"{protein.emdb_id}\t{protein.sample_id}\t{pdbekb.link}\t{pdbekb.provenance}"
                    pdbekb_logger.info(row)
            if self.is_alphafold and protein.alphafold:
                for af in protein.alphafold:
                    row = f"{protein.emdb_id}\t{protein.sample_id}\t{af.link}\t{af.provenance}"
                    alphafold_logger.info(row)
