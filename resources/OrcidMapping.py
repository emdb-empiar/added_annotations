import json, re
import lxml.etree as ET
import xmltodict
import requests
import time

class OrcidMapping:
    """
    Orcid annotations both from author and API.
    """

    def __init__(self, citations, pmc_api, emdb_pubmed):
        self.citations = citations
        self.api = pmc_api
        self.emdb_pubmed = emdb_pubmed

    def execute(self):
        for citation in self.citations:
            citation = self.worker(citation)
        return self.citations

    def worker(self, citation):
        ############# TEST ORCID API ACCESS #######################
        pubmed_file = open(self.emdb_pubmed)
        for line in pubmed_file.readlines():
            if citation.emdb_id in line:
                pubmedid = line.split('\t')[1]
                if pubmedid != "None":
                    time.sleep(0.5)
                    citation.orcid_ids = self.orcid_from_Epmc(pubmedid)
        return citation

    def pmc_api_query(self, queryString):
        pm_id = None
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

    def orcid_from_Epmc(self, pubmed_id):
        for citation in self.citations:
            orcid_ids = {}
            ids = set()
            url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={pubmed_id}&resultType=core&cursorMark=*&pageSize=25&fromSearchPost=false"
            response = requests.get(url)
            if response.status_code == 200 and response.content:
                root = ET.fromstring(response.content)
                if list(root.iter('authorIdList')):
                    for x in list(root.iter('authorIdList')):
                        for auth in x.iter('authorId'):
                            if auth.attrib['type'] == 'ORCID':
                                orcid_id = auth.text
                                if orcid_id not in ids:
                                    if list(root.iter('authorList')):
                                        for y in list(root.iter('author')):
                                            for id in y.iter('authorId'):
                                                if id.attrib['type'] == 'ORCID':
                                                    if orcid_id == id.text:
                                                        if y.find('firstName') is not None:
                                                            firstname = y.find('firstName').text
                                                            if y.find('lastName') is not None:
                                                                lastname = y.find('lastName').text
                                                                author_name = firstname + " " + lastname
                                                                orcid_ids[author_name] = orcid_id
                                                                citation.provenance_orcid = "EuropePMC"
                                    else:
                                        author_name = self.name_for_orcid_id(orcid_id)
                                        orcid_ids[author_name] = orcid_id
                                        citation.provenance_orcid = "ORCID"
                                ids.add(orcid_id)
                                citation.orcid_ids = ids
            return orcid_ids

    def name_for_orcid_id(self, orcid_id):
        author_name = ""
        url_name = f"https://pub.orcid.org/v3.0/{orcid_id}"
        response = requests.get(url_name)
        if response.status_code == 200 and response.content:
            xmld = xmltodict.parse(response.text)
            xmld = xmld['record:record']
            if 'person:person' in xmld:
                if 'person:name' in xmld['person:person']:
                    given_name = xmld['person:person']['person:name']['personal-details:given-names']
                    family_name = xmld['person:person']['person:name']['personal-details:family-name']
                    author_name = given_name + " " + family_name
        return author_name

    def export_tsv(self, orcid_logger):
        for citation in self.citations:
            for name, id in (citation.orcid_ids).items():
                row = f"{citation.emdb_id}\t{name}\t{id}\t{citation.provenance_orcid}"
                orcid_logger.info(row)
