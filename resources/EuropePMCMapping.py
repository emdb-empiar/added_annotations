import csv
from multiprocessing import Pool
import pickle, copyreg
import json
import ssl
from urllib.request import urlopen

pmc_ftp = r'/Users/amudha/project/ftp_data/PMC/PMID_PMCID_DOI.csv'
# pmc_ftp = r'/nfs/ftp/pub/databases/pmc/DOI/PMID_PMCID_DOI.csv.gz'
pmc_baseurl = r'https://www.ebi.ac.uk/europepmc/webservices/rest/search?'
pmc_append = r'%22&resultType=lite&pageSize=25&format=json'

class EuropePMCMapping:
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
        if citation.pmedid or citation.doi or citation.issn:
            citation.provenance = "AUTHOR"
        else:
            citation.provenance = "EuropePMC"
        if not citation.pmedid:
            if citation.doi:
                with open(pmc_ftp, 'r') as f:
                    reader = csv.reader(f, delimiter=',')
                    next(reader, None)  # skip the headers
                    batch_data = list(reader)
                    for line in batch_data:
                        pmid = line[0]
                        doi = line[2]
                        if citation.doi == doi:
                            citation.pmedid = pmid
                            citation.provenance = "EuropePMC"
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

                def save_sslcontext(obj):
                    return obj.__class__, (obj.protocol,)
                copyreg.pickle(ssl.SSLContext, save_sslcontext)
                context = ssl.create_default_context()
                foo = pickle.dumps(context)
                gcontext = pickle.loads(foo)
                pmcjson = urlopen(url, context=gcontext).read()
                pmcjdata = json.loads(pmcjson.decode('utf-8'))
                print(pmcjdata)

                id = pmcjdata['resultList']['result']
                if id:
                    pm_id = pmcjdata['resultList']['result'][0]['id']
                    citation.pmedid = pm_id
                citation.provenance = "EuropePMC"
        print(citation.__dict__)
        return citation



