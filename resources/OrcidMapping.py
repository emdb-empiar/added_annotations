import json, os
import lxml.etree as ET
import xmltodict
import requests
import time
import collections

class OrcidMapping:
    """
    Orcid annotations both from author and API.
    """

    def __init__(self, citations, pmc_api, workDir):
        self.citations = citations
        self.api = pmc_api
        self.workDir = workDir

    def execute(self):
        for citation in self.citations:
            citation = self.worker(citation)
        return self.citations

    def worker(self, citation):
        ############# TEST ORCID API ACCESS #######################
        emdb_pubmed = os.path.join(self.workDir, "emdb_pubmed.log")
        pubmed_file = open(emdb_pubmed)
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
                if list(root.iter('result')):
                    for a in list(root.iter('result')):
                        for pmid in a.iter('pmid'):
                            if pmid.text == pubmed_id:
                                if list(a.iter('authorIdList')):
                                    for x in list(a.iter('authorIdList')):
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
                                                                                if lastname in citation.author_order:
                                                                                    order = citation.author_order[lastname]
                                                                                    name_order = f'{order} [{author_name}]'
                                                                                    orcid_ids[name_order] = orcid_id
                                                                                else:
                                                                                    name_order = f'0 [{author_name}]'
                                                                                    orcid_ids[name_order] = orcid_id
                                                                                citation.provenance_orcid = "EuropePMC"
                                                    else:
                                                        author_name, lastname = self.name_for_orcid_id(orcid_id)
                                                        if lastname in citation.author_order:
                                                            order = citation.author_order[lastname]
                                                            name_order = f'{order} [{author_name}]'
                                                            orcid_ids[name_order] = orcid_id
                                                        else:
                                                            name_order = f'0 [{author_name}]'
                                                            orcid_ids[name_order] = orcid_id
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
        return author_name, family_name

    def export_tsv(self, orcid_logger):
        for citation in self.citations:
            sort_orcid_ids = collections.OrderedDict(sorted((citation.orcid_ids).items(), key=lambda item: int(item[0].split("[")[0])))
            for auth_order, id in (sort_orcid_ids).items():
                order = auth_order.split('[')[0]
                if '[' in auth_order:
                    name = auth_order.split('[', 1)[1].split(']')[0]

                row = f"{citation.emdb_id}\t{name}\t{id}\t{order}\t{citation.provenance_orcid}"
                orcid_logger.info(row)
