import os, csv
from multiprocessing import Pool
import json
import urllib3
import gzip
import shutil

pmc_baseurl = r'https://www.ebi.ac.uk/europepmc/webservices/rest/search?'
pmc_append = r'%22&resultType=lite&pageSize=25&format=json'

class PubmedMapping:
    """
    Author provided publication IDs (PUBMED, DOI) and querying PMC API for title if any publication IDs available in 
    EuropePMC if not provided by author.
    """

    def __init__(self, workDir, citations, pmc_ftp, pmc_ftp_gz):
        self.workDir = workDir
        self.citations = citations
        self.pmc_ftp = pmc_ftp
        self.pmc_ftp_gz = pmc_ftp_gz

        self.pm_doi, self.pm_pmc = self.pm_doi_dict()

    def execute(self, threads):
        with Pool(processes=threads) as pool:
            self.citations = pool.map(self.worker, self.citations)
        return self.citations

    def worker(self, citation):
        if citation.pmcid:
            citation.provenance_pmc = "AUTHOR"

        if citation.pmedid:
            citation.provenance_pm = "AUTHOR"
            citation.url = "https://pubmed.ncbi.nlm.nih.gov/" + citation.pmedid + "/"
            if citation.pmedid in self.pm_pmc:
                citation.pmcid = self.pm_pmc[citation.pmedid]
                citation.provenance_pmc = "EuropePMC"
        else:
            if citation.doi:
                doi = "https://doi.org/" + citation.doi
                if doi in self.pm_doi:
                    citation.pmedid = self.pm_doi[doi]
                    citation.provenance_pm = "EuropePMC"
                    if citation.pmedid in self.pm_pmc:
                        citation.pmcid = self.pm_pmc[citation.pmedid]
                        citation.provenance_pmc = "EuropePMC"

        if citation.doi:
            citation.provenance_doi = "AUTHOR"
        else:
            if citation.pmedid:
                if citation.pmedid in self.pm_doi:
                    citation.doi = self.pm_doi[citation.pmedid]
                    citation.provenance_doi = "EuropePMC"

        return citation

    def pm_doi_dict(self):
        """
        Extract if PMID, PMCID, if DOI exists for a publication from PMC's ftp file and convert it to dictionary
        """
        if not os.path.exists(self.pmc_ftp):
            with gzip.open(self.pmc_ftp_gz, 'rb') as f_in:
                with open(self.pmc_ftp, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        pm_doi = {}
        pm_pmc = {}
        with open(self.pmc_ftp, 'r') as f:
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
        