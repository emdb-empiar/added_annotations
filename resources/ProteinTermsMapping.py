import requests
from models import GO, Interpro, Pfam
import lxml.etree as ET

uni_api = r'https://www.uniprot.org/uniprot/'

class ProteinTermsMapping:
    """
    Extract GO, InterPro and Pfam terms from the UniProt ID mapping web API
    """

    def __init__(self, proteins, is_go=True, is_interpro=True, is_pfam=True):
        self.proteins = proteins
        self.is_interpro = is_interpro
        self.is_go = is_go
        self.is_pfam = is_pfam

    def execute(self):
        for protein in self.proteins:
            if protein.uniprot_id:
                protein = self.worker(protein)
        return self.proteins

    def worker(self, protein):
        return self.uniprot_api(protein)

    def uniprot_api(self, protein):
        uid = protein.uniprot_id
        url = uni_api + uid + ".xml"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            root = ET.fromstring(response.content)

            if self.is_go:
                go_elements = root.findall(".//{http://uniprot.org/uniprot}dbReference[@type='GO']")
                for element in go_elements:
                    go = GO()
                    go.id = element.get("id")
                    terms = element.findall("{http://uniprot.org/uniprot}property[@type='term']")
                    if terms:
                        term_text = terms[0].get("value")
                        go.type = term_text[0]
                        go.namespace = term_text[2:]
                        go.unip_id = uid
                        go.provenance = "UNIPROT"
                        protein.go.append(go)

            if self.is_interpro:
                interpro_elements = root.findall(".//{http://uniprot.org/uniprot}dbReference[@type='InterPro']")
                for element in interpro_elements:
                    interpro = Interpro()
                    interpro.id = element.get("id")
                    terms = element.findall("{http://uniprot.org/uniprot}property[@type='entry name']")
                    if terms:
                        interpro.namespace = terms[0].get("value")
                        interpro.provenance = "UNIPROT"
                        interpro.unip_id = uid
                        protein.interpro.append(interpro)

            if self.is_pfam:
                pfam_elements = root.findall(".//{http://uniprot.org/uniprot}dbReference[@type='Pfam']")
                for element in pfam_elements:
                    pfam = Pfam()
                    pfam.id = element.get("id")
                    terms = element.findall("{http://uniprot.org/uniprot}property[@type='entry name']")
                    if terms:
                        pfam.namespace = terms[0].get("value")
                        pfam.provenance = "UNIPROT"
                        pfam.unip_id = uid
                        protein.pfam.append(pfam)
        return protein

    def export_tsv(self, go_logger, interpro_logger, pfam_logger):
        for protein in self.proteins:
            if self.is_go and protein.go:
                for go in protein.go:
                    row = f"{protein.emdb_id}\t{protein.sample_id}\t{go.id}\t{go.namespace}\t{go.type}\t{go.provenance}"
                    go_logger.info(row)
            if self.is_interpro and protein.interpro:
                for ipr in protein.interpro:
                    row = f"{protein.emdb_id}\t{protein.sample_id}\t{ipr.id}\t{ipr.namespace}\t{ipr.provenance}"
                    interpro_logger.info(row)
            if self.is_pfam and protein.pfam:
                for pfam in protein.pfam:
                    row = f"{protein.emdb_id}\t{protein.sample_id}\t{pfam.id}\t{pfam.namespace}\t{pfam.provenance}"
                    pfam_logger.info(row)


