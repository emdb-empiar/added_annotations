import json
import lxml.etree as ET
import xmltodict
import requests

class PubmedMapping:
    """
    Author provided publication IDs (PUBMED, DOI) and querying PMC API for title if any publication IDs available in
    EuropePMC if not provided by author.
    """

    def __init__(self, citations, pmc_api):
        self.citations = citations
        self.api = pmc_api

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
            self.oricid_for_pubmed(citation.pmedid)
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

    def oricid_for_pubmed(self, pubmed_id):
        orcid_ids = []
        url = f"https://pub.orcid.org/v3.0/search/?q=pmid:{pubmed_id}"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            data = xmltodict.parse(response.text)
            datas = data['search:search']
            num = datas['@num-found']
            if 'search:result' in datas:
                for x in range(int(num)):
                    orcid_id = datas['search:result'][x]['common:orcid-identifier']['common:path']
                    orcid_ids.append(orcid_id)

        return orcid_ids


