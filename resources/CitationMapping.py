import json, os
import requests
from glob import glob
import lxml.etree as ET
from models import Citation

class CitationMapping:
    """
    Author provided publication IDs (PUBMED, DOI) and querying PMC API for PubMed, PubMED central, DOI and ORCID using
    any publication IDs or title available in EuropePMC if not provided by author.
    """

    def __init__(self, headerDir, pmc_api, pubmed_logger, orcid_logger):
        self.headerDir = headerDir
        self.citations = []
        self.emdb_pubmed_ids = {}
        self.pmc_api = pmc_api
        self.pubmed_logger = pubmed_logger
        self.orcid_logger = orcid_logger

    def execute(self):
        pubmed_set = self.get_pubmed_ids()[0]
        for i in self.citations:
            self.emdb_pubmed_ids.update(i.pubmed_ids)
        pubmed_list = list(pubmed_set)
        n = 3
        pm = [pubmed_list[i:i + n] for i in range(0, len(pubmed_list), n)]
        for i in pm:
            pms = ' '.join([str(elem) for elem in i])
            pubmeds = pms.replace(" ", " OR ")
            webAPI = self.pmc_api_query(pubmeds)
            for cite in webAPI:
                for emdb_id, value in self.emdb_pubmed_ids.items():
                    if cite[0] == value:
                        row = f"{emdb_id}\t{cite[0]}\t{cite[1]}\t{cite[3]}\t{cite[2]}"
                        self.pubmed_logger.info(row)
                        ind = cite[4]["ind"]
                        for i in range(int(ind)):
                            name = cite[4]["name_" + str(i)]
                            id = cite[4]["id_" + str(i)]
                            order = cite[4]["order_" + str(i)]
                            provenance = cite[4]["provenance_orcid_" + str(i)]
                            row = f"{emdb_id}\t{name}\t{id}\t{order}\t{provenance}"
                            self.orcid_logger.info(row)
        return

    def get_pubmed_ids(self):
        pubmed_ids = set()
        for fn in glob(os.path.join(str(self.headerDir), '*')):
            id_num = fn.split('-')[1]
            xml_filename = "emd-" + id_num + "-v30.xml"
            xml_dirpath = os.path.join(str(self.headerDir),fn, "header")
            xml_filepath = os.path.join(xml_dirpath, xml_filename)
            with open(xml_filepath, 'r') as filexml:
                tree = ET.parse(filexml)
                root = tree.getroot()
                a = root.attrib
                self.emdb_id = a.get('emdb_id')
            if list(root.iter('primary_citation')):
                for y in list(root.iter('primary_citation')):
                    citation = Citation(self.emdb_id)
                    pub = y.find('journal_citation')
                    ind = 0
                    for auth in y.iter('author'):
                        author = auth.text
                        author_lastname = author.split(" ")[0]
                        citation.orcid_ids["name_" + str(ind)] = author
                        citation.orcid_ids["provenance_orcid_" + str(ind)] = "EMDB"
                        if 'order' in auth.attrib:
                            auth_order = auth.attrib['order']
                            citation.orcid_ids["order_" + str(ind)] = auth_order
                        if 'ORCID' in auth.attrib:
                            orcid_id = auth.attrib['ORCID']
                            citation.orcid_ids["id_" + str(ind)] = orcid_id
                        else:
                            citation.orcid_ids["id_" + str(ind)] = "N/A"
                        ind = ind + 1

                    citation.orcid_ids["ind"] = ind
                    nas = pub.find('title').text
                    title = nas.split('\n\n', 1)[0]
                    citation.title = title
                    pubStatus = pub.attrib['published']
                    if pubStatus == 'true':
                        citation.status = "published"
                    if pubStatus is None:
                        continue
                    for child in pub:
                        pmedid = child.text
                        pmedty = (child.attrib).get('type')
                        if pmedty is not None:
                            if pmedty == 'PUBMED':
                                citation.pmedid = pmedid
                                pubmed_ids.add(pmedid)
                                citation.pubmed_ids[self.emdb_id] = pmedid
                            if pmedty == 'DOI':
                                doi = pmedid.split(":")[1]
                                citation.doi = doi
                            if pmedty == 'ISSN':
                                citation.issn = pmedid
                    self.citations.append(citation)
        return pubmed_ids, self.citations

    def pmc_api_query(self, pubmeds):
        orcid_ids = {}
        pm_id = pmc_id = doi = None
        queryString = "SRC:MED AND EXT_ID:(" + str(pubmeds) + ")"
        data = {'query': queryString,
                'format': 'json',
                'resultType': 'core'}
        response = requests.post(self.pmc_api, data=data)
        res_text = response.text
        try:
            pmcjdata = json.loads(res_text)
        except json.JSONDecodeError:
            yield "", "", "", "", ""
        if 'resultList' not in pmcjdata:
            yield "", "", "", "", ""
        hit = pmcjdata['hitCount']
        result = pmcjdata['resultList']['result']
        if result:
            for i in range(hit):
                source = result[i]['source']
                if source == 'MED':
                    pm_id = result[i]['id']
                    pmc_id = result[i]['pmcid'] if 'pmcid' in result[i] else ""
                    doi = result[i]['doi'] if 'doi' in result[i] else ""
                res = pmcjdata['resultList']['result'][i]
                if 'journalInfo' in res:
                    ji = res['journalInfo']
                    if 'journal' in ji:
                        issn = ji['journal']['issn'] if 'issn' in ji['journal'] else ""
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
                yield pm_id, pmc_id, doi, issn, orcid_ids
        return