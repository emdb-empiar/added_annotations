import json
import requests
import urllib3
from multiprocessing import Pool

pmc_geturl = r'https://www.ebi.ac.uk/europepmc/webservices/rest/search?'
pmc_posturl = r'https://www.ebi.ac.uk/europepmc/webservices/rest/searchPOST'
pmc_append = r'%22&resultType=lite&pageSize=25&format=json'

class PubmedMapping:
    """
    Author provided publication IDs (PUBMED, DOI) and querying PMC API for title if any publication IDs available in
    EuropePMC if not provided by author.
    """

    def __init__(self, workDir, citations):
        self.workDir = workDir
        self.citations = citations

    def execute(self, threads):
        with Pool(processes=threads) as pool:
            self.citations = pool.map(self.worker, self.citations)
        return self.citations

    def worker(self, citation):
        if citation.pmcid is not None:
            citation.provenance_pmc = "AUTHOR"

        if citation.doi is not None:
            citation.provenance_doi = "AUTHOR"

        if citation.pmedid is not None:
            citation.provenance_pm = "AUTHOR"
            if citation.doi is None or citation.pmcid is None:
                webAPI = self.pmc_api_query(("ext_id:" + citation.pmedid))
                citation.pmcid = webAPI[1]
                citation.provenance_pmc = "EuropePMC"
                citation.doi = webAPI[2]
                citation.provenance_doi = "EuropePMC"

        if citation.pmedid is None:
            if citation.doi is not None:
                webAPI = self.pmc_api_query(("DOI:" + citation.doi))
                if webAPI[0]:
                    citation.pmedid = webAPI[0]
                    citation.provenance_pm = "EuropePMC"
                if webAPI[1]:
                    citation.pmcid = webAPI[1]
                    citation.provenance_pmc = "EuropePMC"

        if not citation.pmedid and not citation.doi:
            if citation.title:
                queryString = (citation.title).replace("%", "%25")
                queryString = queryString.replace(' ', '%20')
                queryString = queryString.replace("\n", "%0A")
                queryString = queryString.replace("=", "%3D")
                queryString = queryString.replace("(", "%28")
                queryString = queryString.replace(")", "%29")
                queryString = queryString.replace(",", "%2C")
                queryString = queryString.replace("-", "%2D")
                queryString = queryString.replace("&#183;", "%2D")
                queryString = queryString.replace("&#966;", "%CF%86")
                queryString = queryString.replace("/", "%2F")
                url = pmc_geturl + "query=%22" + (queryString) + pmc_append
                http = urllib3.PoolManager()
                response = http.request('GET', url)
                data = response.data
                pmcjdata = json.loads(data)
                id = pmcjdata['resultList']['result']
                if id:
                    pm_id = pmcjdata['resultList']['result'][0]['id']
                    citation.pmedid = pm_id
                    citation.provenance_pm = "EuropePMC"
                    pmc_id = pmcjdata['resultList']['result'][0]['pmcid']
                    citation.pmcid = pmc_id
                    citation.provenance_pmc = "EuropePMC"
                    doi_id = pmcjdata['resultList']['result'][0]['doi']
                    citation.doi= doi_id
                    citation.provenance_doi = "EuropePMC"
        # print(citation.__dict__)
        return citation

    def pmc_api_query(self, queryString):
        pm_id = pmc_id = doi = None
        data = {'query': queryString,
                'format': 'json',
                'resultType': 'lite'}
        response = requests.post(pmc_posturl, data=data)
        res_text = response.text
        pmcjdata = json.loads(res_text)
        id = pmcjdata['resultList']['result']
        if id:
            source = pmcjdata['resultList']['result'][0]['source']
            if source == 'MED':
                pm_id = pmcjdata['resultList']['result'][0]['id']
                pmc_id = pmcjdata['resultList']['result'][0]['pmcid']
                doi = pmcjdata['resultList']['result'][0]['doi']
        return pm_id, pmc_id, doi
