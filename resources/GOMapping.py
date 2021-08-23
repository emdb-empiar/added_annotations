import csv
import urllib3
import lxml.etree as ET

uni_api = r'https://www.uniprot.org/uniprot/'

class GOMapping:
    """
    If map+model entry extract the GO ids from PDBe sifts file. If map only entry query in EuropePMC for 
    GO ids by the publication
    """

    def __init__(self, workDir, GOs, shifts_GO, GO_obo, uniprot_ids):
        self.workDir = workDir
        self.GOs = GOs
        self.shifts_GO = shifts_GO
        self.GO_obo = GO_obo
        self.uniprot_ids = uniprot_ids

        self.pdbe_GO = self.extract_resources_from_sifts()
        self.obo_dict = self.go_namespace()

    def execute(self):
        for GO in self.GOs:
            GO = self.worker(GO)
        return self.GOs

    def worker(self, GO):
        if GO.provenance != "AUTHOR":
            if GO.pdb_id:
                if GO.pdb_id in self.pdbe_GO:
                    GO.GO_id = self.pdbe_GO[GO.pdb_id]
                    for go_id in GO.GO_id:
                        if go_id in self.obo_dict:
                            go_namespace = self.obo_dict[go_id]
                            (GO.GO_namespace).append(go_namespace)
                    GO.provenance = "PDBe"
            else:
                if self.uniprot_ids:
                    for uid in self.uniprot_ids:
                        ids = self.uniprot_api(uid)
                        GO.GO_id = ids
                        for go_id in GO.GO_id:
                            if go_id in self.obo_dict:
                                go_namespace = self.obo_dict[go_id]
                                (GO.GO_namespace).append(go_namespace)
                        GO.provenance = "UNIPROT"
        # print(GO.__dict__)
        return GO

    def extract_resources_from_sifts(self):
        """
        Extract the GO ids for every pdb_id from pdbe sifts file
        """
        pdb_GO = {}
        with open(self.shifts_GO, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            next(reader, None)
            for row in reader:
                if row[0] not in pdb_GO:
                    pdb_GO[row[0]] = []
                else:
                    pdb_GO[row[0]].append(row[-1])
        pdbe_GO = {a: list(set(b)) for a, b in pdb_GO.items()}
        return pdbe_GO

    def uniprot_api(self, uid):
        """
        Query for GO_id if Uniprot_id exists for map only entries
        """

        ids = set()
        url = uni_api + uid + ".xml"
        http = urllib3.PoolManager()
        response = http.request('GET', url)
        data = response.data
        if data:
            tree = ET.parse(data)
            root = tree.getroot()
            if list(root.iter('dbReference')):
                for x in list(root.iter('dbReference')):
                    if x.attrib == "GO":
                        go_id = x.get('id').text
                        if go_id:
                            ids.add(go_id)
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



