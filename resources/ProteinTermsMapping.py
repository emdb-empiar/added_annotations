import requests
from models import GO, Interpro
import lxml.etree as ET

uni_api = r'https://www.uniprot.org/uniprot/'

class ProteinTermsMapping:
    """
    Extract GO and InterPro terms from the UniProt ID mapping web API
    """

    def __init__(self, proteins, is_go=True, is_interpro=True):
        self.proteins = proteins
        self.is_interpro = is_interpro
        self.is_go = is_go


    def execute(self):
        for protein in self.proteins:
            if protein.uniprot_id:
                protein = self.worker(protein)
        return self.proteins

    def worker(self, protein):
        #TODO: Add translator to obtain namespace and type from author submited GO and InterPro
        return self.uniprot_api(protein) 

    def uniprot_api(self, protein):
        ids = set()
        uid = protein.uniprot_id
        url = uni_api + uid + ".xml"
        response = requests.get(url)
        if response.status_code == 200:   
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
                    protein.interpro.append(interpro)

        return protein

    def export_tsv(self, go_logger, interpro_logger):
        for protein in self.proteins:
            if self.is_go and protein.go:
                for go in protein.go:
                    row = f"{protein.emdb_id}\t{protein.sample_id}\t{go.id}\t{go.namespace}\t{go.type}\t{go.provenance}"
                    go_logger.info(row)
            if self.is_interpro and protein.interpro:
                for ipr in protein.interpro:
                    row = f"{protein.emdb_id}\t{protein.sample_id}\t{ipr.id}\t{ipr.namespace}\t{ipr.provenance}"
                    interpro_logger.info(row)

