import re
import requests

uni_api = r'https://www.uniprot.org/uniprot/'

class GOMapping:
    """
    If map+model entry extract the GO ids from PDBe sifts file. If map only entry query in EuropePMC for 
    GO ids by the publication
    """

    def __init__(self, workDir, GOs, GO_obo, uniprot_ids):
        self.workDir = workDir
        self.GOs = GOs
        self.GO_obo = GO_obo
        self.uniprot_ids = uniprot_ids

        self.obo_dict = self.go_namespace()

    def execute(self):
        for GO in self.GOs:
            GO = self.worker(GO)
        return self.GOs

    def worker(self, GO):
        if GO.provenance == "AUTHOR":
            for go_id in GO.GO_id:
                if go_id in self.obo_dict:
                    go_namespace = self.obo_dict[go_id]
                    (GO.GO_namespace).append(go_namespace)
                else:
                    (GO.GO_namespace).append("N/A")
        else:
            if self.uniprot_ids:
                for uid in self.uniprot_ids:
                    ids = self.uniprot_api(uid)
                    GO.GO_id = ids
                    for go_id in GO.GO_id:
                        if go_id in self.obo_dict:
                            go_namespace = self.obo_dict[go_id]
                            (GO.GO_namespace).append(go_namespace)
                        else:
                            (GO.GO_namespace).append("N/A")
                    GO.provenance = "UNIPROT"
        print(GO.__dict__)
        return GO

    def uniprot_api(self, uid):
        """
        Query for GO_id if Uniprot_id exists for map only entries
        """

        ids = set()
        url = uni_api + uid + ".txt"
        response = requests.get(url)
        data = response.text
        if data:
            ids = re.findall(r'GO:\d+', response.text)
        return ids

    def go_namespace(self):
        """
        For every GO_ID get the corresponding GO_namespace
        """
        obo_dict = {}
        with open(self.GO_obo, 'r') as f:
            for line in f:
                if "[Term]" in line:
                    id = (next(f).strip()).split('id: ')[1]
                    name = (next(f).strip())
                    name_space = (next(f).strip()).split(': ')[1]
                    obo_dict[str(id)] = str(name_space)
        return obo_dict



