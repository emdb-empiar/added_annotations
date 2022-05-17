import json
import requests

class PublicationMapping:
    """
    Author provided publication IDs (PUBMED, DOI) and querying PMC API for PubMed, PubMED central, DOI and ORCID using
    any publication IDs or title available in EuropePMC if not provided by author.
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
            citation.provenance_pmc = "EMDB"

        if citation.doi:
            citation.provenance_doi = "EMDB"

        if citation.pmedid:
            citation.provenance_pm = "EMDB"
            if not citation.doi or not citation.pmcid:
                webAPI = self.pmc_api_query(("ext_id:" + citation.pmedid))
                if not citation.pmcid:
                    if webAPI[1]:
                        citation.pmcid = webAPI[1]
                        citation.provenance_pmc = "EuropePMC"
                if not citation.doi:
                    if webAPI[2]:
                        citation.doi = webAPI[2]
                        citation.provenance_doi = "EuropePMC"
                if self.is_orcid:
                    if webAPI[3]:
                        citation.orcid_ids = webAPI[3]

        else:
            if citation.doi:
                webAPI = self.pmc_api_query(("DOI:" + citation.doi))
                if not citation.pmedid:
                    if webAPI[0]:
                        citation.pmedid = webAPI[0]
                        citation.provenance_pm = "EuropePMC"
                if not citation.pmcid:
                    if webAPI[1]:
                        citation.pmcid = webAPI[1]
                        citation.provenance_pmc = "EuropePMC"
                if self.is_orcid:
                    if webAPI[3]:
                        citation.orcid_ids = webAPI[3]

        if not citation.pmedid and not citation.doi:
            if citation.title:
                webAPI = self.pmc_api_query((f'TITLE:"{citation.title}"'))
                if not citation.pmedid:
                    if webAPI[0]:
                        citation.pmedid = webAPI[0]
                        citation.provenance_pm = "EuropePMC"
                if not citation.pmcid:
                    if webAPI[1]:
                        citation.pmcid = webAPI[1]
                        citation.provenance_pmc = "EuropePMC"
                if not citation.doi:
                    if webAPI[2]:
                        citation.doi = webAPI[2]
                        citation.provenance_doi = "EuropePMC"
                if self.is_orcid:
                    if webAPI[3]:
                        citation.orcid_ids = webAPI[3]

        return citation

    def pmc_api_query(self, queryString):
        orcid_ids = {}
        pm_id = pmc_id = doi = None
        data = {'query': queryString,
                'format': 'json',
                'resultType': 'core'}
        response = requests.post(self.api, data=data)
        res_text = response.text
        try:
            pmcjdata = json.loads(res_text)
        except json.JSONDecodeError:
            return "", "", ""
        if 'resultList' not in pmcjdata:
            return "", "", ""
        result = pmcjdata['resultList']['result']
        if result:
            source = result[0]['source']
            if source == 'MED':
                pm_id = result[0]['id']
                pmc_id = result[0]['pmcid'] if 'pmcid' in result[0] else ""
                doi = result[0]['doi'] if 'doi' in result[0] else ""
            res = pmcjdata['resultList']['result'][0]
            if 'authorList' in res:
                auth = res['authorList']['author']
                size = len(auth)
                for ind in range(size):
                    auth_list = auth[ind]
                    firstname = auth_list['firstName']
                    lastname = auth_list['lastName']
                    author_name = firstname + " " + lastname
                    orcid_ids["name_" + str(ind)] = author_name
                    orcid_ids["order_" + str(ind)] = ind + 1
                    orcid_ids["provenance_orcid_" + str(ind)] = "EuropePMC"
                    if 'authorId' in auth_list:
                        type_id = auth_list['authorId']['type']
                        if type_id == "ORCID":
                            orcid_id = auth_list['authorId']['value']
                            orcid_ids["id_" + str(ind)] = orcid_id
                    else:
                        orcid_ids["id_" + str(ind)] = "N/A"
                orcid_ids["ind"] = size

        return pm_id, pmc_id, doi, orcid_ids

    def export_tsv(self, pubmed_logger, orcid_logger):
        for citation in self.citations:
            row = f"{citation.emdb_id}\t{citation.pmedid}\t{citation.pmcid}\t{citation.issn}\t{citation.doi}"
            pubmed_logger.info(row)
            ind = citation.orcid_ids["ind"]
            for i in range(int(ind)):
                name = citation.orcid_ids["name_" + str(i)]
                id = citation.orcid_ids["id_" + str(i)]
                order = citation.orcid_ids["order_" + str(i)]
                provenance = citation.orcid_ids["provenance_orcid_" + str(i)]
                row = f"{citation.emdb_id}\t{name}\t{id}\t{order}\t{provenance}"
                orcid_logger.info(row)