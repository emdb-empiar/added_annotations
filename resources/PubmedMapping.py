import json, csv
import requests

class PubmedMapping:
    """
    Author provided publication IDs (PUBMED, DOI) and querying PMC API for title if any publication IDs available in
    EuropePMC if not provided by author.
    """

    def __init__(self, citations, pmc_api, is_orcid=True):
        self.citations = citations
        self.api = pmc_api
        self.is_orcid = is_orcid

    def execute(self):
        for citation in self.citations:
            citation = self.worker(citation)
        return self.citations

    def worker(self, citation):
        if citation.pmcid:
            citation.provenance_pmc = "AUTHOR"

        if citation.doi:
            citation.provenance_doi = "AUTHOR"

        if citation.pmedid:
            citation.provenance_pm = "AUTHOR"
            if citation.doi is None or citation.pmcid is None:
                webAPI = self.pmc_api_query(("ext_id:" + citation.pmedid))
                if webAPI[1]:
                    citation.pmcid = webAPI[1]
                    citation.provenance_pmc = "EuropePMC"
                if webAPI[2]:
                    citation.doi = webAPI[2]
                    citation.provenance_doi = "EuropePMC"

        else:
            if citation.doi:
                webAPI = self.pmc_api_query(("DOI:" + citation.doi))
                if webAPI[0]:
                    citation.pmedid = webAPI[0]
                    citation.provenance_pm = "EuropePMC"
                if webAPI[1]:
                    citation.pmcid = webAPI[1]
                    citation.provenance_pmc = "EuropePMC"

        if not citation.pmedid and not citation.doi:
            if citation.title:
                webAPI = self.pmc_api_query((f'TITLE:"{citation.title}"'))
                if webAPI[0]:
                    citation.pmedid = webAPI[0]
                    citation.provenance_pm = "EuropePMC"
                if webAPI[1]:
                    citation.pmcid = webAPI[1]
                    citation.provenance_pmc = "EuropePMC"
                if webAPI[2]:
                    citation.doi = webAPI[2]
                    citation.provenance_doi = "EuropePMC"

        if self.is_orcid:
            orcid_dict = {}
            with open('/Users/amudha/project/emdb_orcid.log', 'r') as f:
                for line in f.readlines():
                    if citation.emdb_id in line:
                        name = line.split('\t')[1]
                        id = line.split('\t')[2]
                        pvn = line.split('\t')[3]
                        orcid_dict[id] = name
                        citation.orcid_ids = orcid_dict
                        citation.provenance_orcid = pvn.strip()
        return citation

    def pmc_api_query(self, queryString):
        pm_id = pmc_id = doi = None
        data = {'query': queryString,
                'format': 'json',
                'resultType': 'lite'}
        response = requests.post(self.api, data=data)
        res_text = response.text
        pmcjdata = json.loads(res_text)
        result = pmcjdata['resultList']['result']
        if result:
            source = result[0]['source']
            if source == 'MED':
                pm_id = result[0]['id']
                pmc_id = result[0]['pmcid'] if 'pmcid' in result[0] else ""
                doi = result[0]['doi'] if 'doi' in result[0] else ""
        return pm_id, pmc_id, doi

    def export_tsv(self, pubmed_logger):
        for citation in self.citations:
            row = f"{citation.emdb_id}\t{citation.pmedid}\t{citation.pmcid}\t{citation.doi}"
            pubmed_logger.info(row)

