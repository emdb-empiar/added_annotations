import os, csv
from multiprocessing import Pool
import json
import urllib3
import configparser
import gzip
import shutil

pmc_baseurl = r'https://www.ebi.ac.uk/europepmc/webservices/rest/search?'
pmc_append = r'%22&resultType=lite&pageSize=25&format=json'

class PubmedMapping:
    """
    Author provided publication IDs (PUBMED, DOI) and querying PMC API for title if any publication IDs available in 
    EuropePMC if not provided by author.
    """

    def __init__(self, workDir, citations):
        self.workDir = workDir
        self.citations = citations

        self.pm_doi, self.pm_pmc = self.pm_doi_dict()

    def execute(self, threads):
        with Pool(processes=threads) as pool:
            self.citations = pool.map(self.worker, self.citations)
        return self.citations

    def worker(self, citation):
        if citation.pmcid:
            citation.provenance_pmc = "AUTHOR"

        if citation.doi:
            citation.provenance_doi = "AUTHOR"

        if citation.pmedid:
            citation.provenance_pm = "AUTHOR"
            citation.url = "https://pubmed.ncbi.nlm.nih.gov/" + citation.pmedid + "/"
            if citation.pmedid in self.pm_pmc:
                citation.pmcid = self.pm_pmc[citation.pmedid]
                citation.provenance_pmc = "EuropePMC"

        if not citation.pmedid:
            if citation.doi:
                citation.provenance_doi = "AUTHOR"
                doi = "https://doi.org/" + citation.doi
                if doi in self.pm_doi:
                    citation.pmedid = self.pm_doi[doi]
                    citation.provenance_pm = "EuropePMC"
                    if citation.pmedid in self.pm_pmc:
                        citation.pmcid = self.pm_pmc[citation.pmedid]
                        citation.provenance_pmc = "EuropePMC"

        if not citation.doi:
            if citation.pmedid:
                if citation.pmedid in self.pm_doi:
                    citation.doi = self.pm_doi[citation.pmedid]
                    citation.provenance_doi = "EuropePMC"

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
                url = pmc_baseurl + "query=%22" + (queryString) + pmc_append

                http = urllib3.PoolManager()
                response = http.request('GET', url)
                data = response.data
                pmcjdata = json.loads(data)
                id = pmcjdata['resultList']['result']
                if id:
                    pm_id = pmcjdata['resultList']['result'][0]['id']
                    citation.pmedid = pm_id
                    pmc_id = pmcjdata['resultList']['result'][0]['pmcid']
                    citation.pmcid = pmc_id
                citation.provenance_pm = "EuropePMC"
                citation.provenance_pmc = "EuropePMC"
        # print(citation.__dict__)
        return citation

    def pm_doi_dict(self):
        """
        Extract if PMID, PMCID, if DOI exists for a publication from PMC's ftp file and convert it to dictionary
        """

        config = configparser.ConfigParser()
        config.read(os.path.join(self.workDir, "git_code/added_annotations/config.ini"))
        pmc_ftp_gz = config.get("file_paths", "pmc_ftp_gz")

        if not os.path.exists(config.get("file_paths", "pmc_ftp")):
            with gzip.open(pmc_ftp_gz, 'rb') as f_in:
                with open(config.get("file_paths", "pmc_ftp"), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        pmc_ftp = config.get("file_paths", "pmc_ftp")
        pm_doi = {}
        pm_pmc = {}
        with open(pmc_ftp, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            next(reader, None)
            for row in reader:
                if row[0] and row[2]:
                    pmid = row[0]
                    doi = row[2]
                    pm_doi[doi] = pmid
                    pm_doi[pmid] = doi
                if (row[0] and row[1]) or (row[0] and row[1] and row[2]):
                    pmid = row[0]
                    pmcid = row[1]
                    pm_pmc[pmid] = pmcid
        return pm_doi, pm_pmc