import json, os
import requests
import itertools

def generate_orcid_dictionary(workDir):
    orcid_dict = {}
    emdb_orcid = os.path.join(workDir, "emdb_orcid.log")
    with open(emdb_orcid, 'r') as f:
        for line in f.readlines()[1:]:
            line = line.strip('\n')
            emdb_id = line.split('\t')[0]
            name = line.split('\t')[1]
            id = line.split('\t')[2]
            pvn = line.split('\t')[3]
            if emdb_id not in orcid_dict:
                orcid_dict[emdb_id] = {}
            orc_list = ["name", name, "id", id, "provenance", pvn]
            list_dict = dict(itertools.zip_longest(*[iter(orc_list)] * 2, fillvalue=""))
            for k in list_dict.keys():
                orcid_dict[emdb_id][k] = list_dict[k]
    return orcid_dict

class PubmedMapping:
    """
    Author provided publication IDs (PUBMED, DOI) and querying PMC API for title if any publication IDs available in
    EuropePMC if not provided by author.
    """

    def __init__(self, citations, pmc_api, orcid_dict={}, is_orcid=True):
        self.citations = citations
        self.api = pmc_api
        self.is_orcid = is_orcid
        self.orcid_dict = orcid_dict

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
            orcid_id = {}
            orc = self.orcid_dict.get(citation.emdb_id)
            if orc:
                name = orc.get('name')
                id = orc.get('id')
                orcid_id[id] = name
                citation.orcid_ids = orcid_id
                citation.provenance_orcid = orc.get('provenance')
        print(citation)
        return citation

    def pmc_api_query(self, queryString):
        pm_id = pmc_id = doi = None
        data = {'query': queryString,
                'format': 'json',
                'resultType': 'lite'}
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
        return pm_id, pmc_id, doi

    def export_tsv(self, pubmed_logger):
        for citation in self.citations:
            row = f"{citation.emdb_id}\t{citation.pmedid}\t{citation.pmcid}\t{citation.doi}"
            pubmed_logger.info(row)

