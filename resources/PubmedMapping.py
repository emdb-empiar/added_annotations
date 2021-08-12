import os, csv
from multiprocessing import Pool
import json
import urllib3
import gzip
import shutil

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

        self.pmdic = self.pm_doi_dict()

    def execute(self, threads):
        with Pool(processes=threads) as pool:
            self.citations = pool.map(self.worker, self.citations)
        return self.citations

    def worker(self, citation):
        print(str(citation))
        if not citation.pmcid and not citation.doi:
            if citation.pmedid in self.pmdic:
                pmc, doi = self.pmdic[citation.pmedid]
                if pmc:
                    citation.pmcid = pmc
                    citation.provenance_pmc = "EuropePMC"
                if doi:
                    citation.doi = doi
                    citation.provenance_doi = "EuropePMC"
        else:
            if citation.pmcid:
                citation.provenance_pmc = "AUTHOR"
            else:
                if citation.pmedid in self.pmdic:
                    pmc, doi = self.pmdic[citation.pmedid]
                    if pmc:
                        citation.pmcid = pmc
                        citation.provenance_pmc = "EuropePMC"

            if citation.doi:
                citation.provenance_doi = "AUTHOR"
            else:
                if citation.pmedid in self.pmdic:
                    pmc, doi = self.pmdic[citation.pmedid]
                    if doi:
                        citation.doi = doi
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

        pmdic = {}
        with open(self.pmc_ftp, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            next(reader, None)
            for row in reader:
                pmdic[row[0]] = (row[1], row[2])
        return pmdic
