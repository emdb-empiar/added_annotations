import os, csv, re
import configparser
import urllib3
import json
from multiprocessing import Pool

pmc_annotation_api = r'https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?articleIds='
pmc_baseurl = r'https://www.ebi.ac.uk/europepmc/webservices/rest/search?'
pmc_append = r'%22&resultType=lite&pageSize=25&format=json'

class GOMapping:
    """
    If map+model entry extract the GO ids from PDBe sifts file. If map only entry query in EuropePMC for 
    GO ids by the publication
    """

    def __init__(self, workDir, GOs):
        self.workDir = workDir
        self.GOs = GOs

        self.pdbe_GO = self.extract_resources_from_sifts()
        self.obo_dict = self.go_namespace()

    def execute(self, threads):
        with Pool(processes=threads) as pool:
            self.GOs = pool.map(self.worker, self.GOs)
        return self.GOs

    def worker(self, GO):
        if GO.provenance != "AUTHOR":
            if GO.pdb_id in self.pdbe_GO:
                GO.GO_id = self.pdbe_GO[GO.pdb_id]
                for go_id in GO.GO_id:
                    if go_id in self.obo_dict:
                        go_namespace = self.obo_dict[go_id]
                        (GO.GO_namespace).append(go_namespace)
                GO.provenance = "PDBe"

            if not GO.pdb_id:
                if GO.pubmed_id is not None:
                    ids = self.pmc_annotation(GO.pubmed_id)
                    GO.GO_id = ids
                    GO.provenance = "EuropePMC"
                if GO.pubmed_id is None:
                    if GO.title:
                        queryString = (GO.title).replace("%", "%25")
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
                            GO.pubmed_id = pm_id
                        ids = self.pmc_annotation(GO.pubmed_id)
                        GO.GO_id = ids
                        GO.provenance = "EuropePMC"
        # print(GO.__dict__)
        return GO

    def extract_resources_from_sifts(self):
        """
        Extract the GO ids for every pdb_id from pdbe sifts file
        """

        config = configparser.ConfigParser()
        config.read(os.path.join(self.workDir, "git_code/added_annotations/config.ini"))
        shifts_GO = config.get("file_paths", "sifts_GO")

        pdb_GO = {}
        with open(shifts_GO, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            next(reader, None)
            for row in reader:
                if row[0] not in pdb_GO:
                    pdb_GO[row[0]] = []
                else:
                    pdb_GO[row[0]].append(row[-1])
        pdbe_GO = {a: list(set(b)) for a, b in pdb_GO.items()}
        return pdbe_GO

    def pmc_annotation(self, id):
        """
        Query for GO_id if pubmed id exists in EuropePMC
        """

        ids = set()
        url = pmc_annotation_api + 'MED%3A' + str(id) + '&format=JSON'
        http = urllib3.PoolManager()
        response = http.request('GET', url)
        data = response.data
        pmcjdata = json.loads(data)
        if pmcjdata:
            for result in re.findall("go/(.*?)'", str(pmcjdata)):
                ident = ",".join(re.findall(r'GO:\d*', result))
                if ident:
                    ids.add(ident)
        return ids

    def go_namespace(self):
        """
        For every GO_ID get the corresponding GO_namespace
        """

        config = configparser.ConfigParser()
        config.read(os.path.join(self.workDir, "git_code/added_annotations/config.ini"))
        GO_obo = config.get("file_paths", "GO_obo")

        obo_dict = {}
        with open(GO_obo, 'r') as f:
            for line in f:
                if "[Term]" in line:
                    id = (next(f).strip()).split('id: ')[1]
                    name = (next(f).strip())
                    name_space = (next(f).strip()).split(': ')[1]
                    obo_dict[str(id)] = str(name_space)
        return obo_dict



