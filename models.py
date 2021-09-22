import re
import requests
import json

class Protein:
    """
    Defines the attributes of a protein sample in a EMDB entry
    """
    def __init__(self, emdb_id, sample_id):
        self.emdb_id = emdb_id
        self.sample_id = sample_id
        self.sample_name = ""
        self.sample_organism = None
        self.pdb = []
        self.sample_complexes = []
        self.uniprot_id = None
        self.provenance = None
        self.sequence = ""
        self.sample_copies = ""
        self.go = []
        self.interpro = []
        self.pfam = []
        self.pdbekb = []
        self.alphafold = []

    def __str__(self):
        return "%s (%s)\n%s (%s) %s - %s [%s]\nComplexes: %s\nPDB: \n%s\n%s\n%s\n  %s\n%s\n" % (self.sample_name, self.sample_organism,
                                                                              self.emdb_id, self.sample_id, self.sample_copies,
                                                                              self.uniprot_id, self.provenance, str(self.sample_complexes),
                                                                              str(self.go), str(self.interpro), str(self.pfam),
                                                                              str(self.pdbekb), str(self.alphafold))

    def get_tsv(self):
        complex_str = ';'.join([str(elem) for elem in self.sample_complexes])
        if self.provenance:
            return ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.emdb_id, self.sample_id, self.sample_name, self.sample_copies,
                                                          self.sample_organism, self.uniprot_id, self.provenance, complex_str))
        else:
            return ""

class CPX:
    """
    Complex Portal entry obtained from their FTP area
    """

    def __init__(self, row):
        self.cpx_id = row[0]
        self.name = row[1]
        self.taxonomy = row[3]
        self.uniprot = set()
        self.identifiers = re.sub(r'\(\d+\)', '', row[4]).split('|')
        self.confidence = row[5]
        self.GO = re.sub(r'\(.+?\)', '', row[7]).split('|')
        self.cross_ref = re.findall(r':(.*?)\(', row[8], re.S)

        for idt in self.identifiers:
            if 'CHEBI:' in idt:
                continue
            if '-PRO_' in idt:
                self.uniprot.add(idt.split('-')[0])
                continue
            if  '_' in idt:
                continue
            self.uniprot.add(idt)

class Sample:
    """
    Unique sample along its parents and childrens
    """
    def __init__(self, sample_id, mw=None, copies=1):
        self.id = sample_id
        self.mw = mw
        self.parent = []
        self.children = []
        self.copies = copies

    def add_parent(self, node):
        self.parent.append(node)

    def add_child(self, node):
        self.children.append(node)

    def __str__(self):
        return f"{self.id}: {self.mw} ({self.copies})"

class Supra:
    """
    Defines the attributes of a supra_molecules in a EMDB entry
    """
    def __init__(self, emdb_id, supra_id):
        self.emdb_id = emdb_id
        self.supra_id = supra_id
        self.supra_name = ""
        self.kind = ""

    def __str__(self):
        return "%s\t%s\t%s\t%s" % (self.emdb_id, self.supra_id, self.supra_name, self.kind)

class EMDB_complex:
    """
    EMDB complex sample obtained from the header files in the Uniprot mapping
    """

    def __init__(self, emdb_id, sample_id, supra_name, sample_copies, complex_sample_id):
        self.emdb_id = emdb_id
        self.sample_id = emdb_id+"_"+sample_id
        self.supra_name = supra_name
        self.sample_copies = sample_copies
        self.complex_sample_id = complex_sample_id
        self.cpx_list = []
        self.proteins = set()
        self.provenance = ""
        self.score = 0.0

    def add_protein(self, uniprot_id):
        self.proteins.add(uniprot_id)

class Ligand:
    """
    Defines the attributes of a ligands sample in a EMDB entry
    """
    def __init__(self, emdb_id, sample_id):
        self.emdb_id = emdb_id
        self.sample_id = sample_id
        self.provenance_chebi = ""
        self.provenance_chembl = ""
        self.provenance_drugbank = ""
        self.HET = ""
        self.lig_name = ""
        self.chembl_id = ""
        self.chebi_id = ""
        self.drugbank_id = ""
        self.lig_copies = ""

    def get_chembl_tsv(self):
        if self.chembl_id:
            return "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.emdb_id, self.sample_id, self.HET, self.lig_name,
                                                     self.lig_copies, self.chembl_id, self.provenance_chembl)
        return ""

    def get_chebi_tsv(self):
        if self.chebi_id:
            return "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.emdb_id, self.sample_id, self.HET, self.lig_name,
                                                     self.lig_copies, self.chebi_id, self.provenance_chebi)
        return ""

    def get_drugbank_tsv(self):
        if self.drugbank_id:
            return "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.emdb_id, self.sample_id, self.HET, self.lig_name,
                                                     self.lig_copies, self.drugbank_id, self.provenance_drugbank)
        return ""

class Model:
    """
    Define the PDB model, preferred assembly and molecular weight
    """
    def __init__(self, emdb_id, pdb_id):
        self.emdb_id = emdb_id
        self.pdb_id = pdb_id
        self.assembly = 1
        self.molecular_weight = 0.0 #Dalton

    def __str__(self):
        return ("{}\t{}\t{}\t{:0.3f}".format(self.emdb_id, self.pdb_id, self.assembly, self.molecular_weight))

class Weight:
    """
    Total weight of the sample provided by author
    """
    def __init__(self, emdb_id):
        self.emdb_id = emdb_id
        self.provenance = ""
        self.kind = ""
        self.type = None
        self.method = ""
        self.sup_th_weight = []
        self.sup_th_unit = ""
        self.sup_exp_weight = []
        self.sup_exp_unit = ""
        self.macro_th_weight = []
        self.macro_th_unit = ""
        self.macro_exp_weight = []
        self.macro_exp_unit = ""
        self.sample_th_weight = 0.0
        self.th_unit = ""
        self.sample_exp_weight = 0.0
        self.exp_unit = ""

    def __str__(self):
            return ("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.emdb_id, self.sample_th_weight, self.sample_exp_weight, self.sup_th_weight,
                                                      self.sup_exp_weight, self.macro_th_weight, self.macro_exp_weight))

class Empiar:
    """
    Defines the EMPIAR ID in a EMDB entry
    """
    def __init__(self, emdb_id):
        self.emdb_id = emdb_id
        self.empiar_id = ""

    def __str__(self):
        return "%s\t%s\n" % (self.emdb_id, self.empiar_id)

class Citation:
    """
    Defines the attributes of a publication in a EMDB entry
    """
    def __init__(self, emdb_id):
        self.emdb_id = emdb_id
        self.pmedid = None
        self.pmcid = None
        self.doi = None
        self.issn = None
        self.status = ""
        self.title = ""
        self.provenance_pm = ""
        self.provenance_pmc = ""
        self.provenance_doi = ""
        self.url = ""

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\n" % (self.emdb_id, self.pmedid, self.pmcid, self.doi, self.issn)

class GO:
    """
    Define the GO terms for the sample in the EMDB entry
    """
    def __init__(self):
        self.id = ""
        self.namespace = ""
        self.type = ""
        self.unip_id = ""
        self.provenance = ""

    def __str__(self):
        return f"{self.id}\t{self.namespace}\t{self.type}\t{self.provenance}"

    def add_from_author(self, go_text, unip_id):
        self.provenance = "AUTHOR"
        self.unip_id = unip_id
        if "GO:" in go_text:
            self.id = go_text
        elif go_text.isdigit():
            self.id = f"GO:{go_text}"

        if self.id and not self.namespace:
            url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{self.id}"
            response = requests.get(url)
            if response.status_code == 200:
                res_text = response.text
                data = json.loads(res_text)
                hits = data['numberOfHits']
                if hits > 0:
                    result = data['results'][0]
                    self.namespace = result['name']
                    aspect = result['aspect']
                    if aspect == 'biological_process':
                        self.type = "P"
                    elif aspect == 'cellular_component':
                        self.type = "C"
                    elif aspect == 'molecular_function':
                        self.type = "F"


class Interpro:
    """
    Define the InterPro terms for the sample in the EMDB entry
    """
    def __init__(self):
        self.id = ""
        self.namespace = ""
        self.unip_id = ""
        self.provenance = ""

    def __str__(self):
        return f"{self.id}\t{self.namespace}\t{self.provenance}"

    def add_from_author(self, ipr_text, unip_id):
        self.provenance = "AUTHOR"
        self.unip_id = unip_id
        if "IPR" in ipr_text:
            self.id = ipr_text

        if self.id and not self.namespace:
            url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/{self.id}"
            response = requests.get(url)
            if response.status_code == 200:
                res_text = response.text
                data = json.loads(res_text)
                if 'metadata' in data:
                    result = data['metadata']
                    if 'hierarchy' in result:
                        if 'name' in result['hierarchy']:
                            self.namespace = result['hierarchy']['name']

class Pfam:
    """
    Define the InterPro terms for the sample in the EMDB entry
    """
    def __init__(self):
        self.id = ""
        self.namespace = ""
        self.unip_id = ""
        self.sample_id = ""
        self.provenance = ""

    def __str__(self):
        return f"{self.id}\t{self.namespace}\t{self.provenance}"

    def add_from_author(self, pfam_text, unip_id):
        self.provenance = "AUTHOR"
        self.unip_id = unip_id
        if "PF" in pfam_text:
            self.id = pfam_text

        if self.id and not self.namespace:
            url = f"https://pfam.xfam.org/family/{self.id}?output=xml"
            response = requests.get(url)
            if response.status_code == 200:
                res_text = response.text
                data = json.loads(res_text)
                if 'description' in data:
                    result = data['description']
                    self.namespace = result['description']

class Pdbekb:
    """
    Define the PDBeKB terms for the sample in the EMDB entry
    """
    def __init__(self):
        self.link = ""
        self.unip_id = ""
        self.provenance = ""

    def __str__(self):
        return f"{self.unip_id}\t{self.link}\t{self.provenance}"

class Alphafold:
    """
    Define the Alphafold terms for the sample in the EMDB entry
    """
    def __init__(self):
        self.link = ""
        self.unip_id = ""
        self.provenance = ""

    def __str__(self):
        return f"{self.unip_id}\t{self.link}\t{self.provenance}"

