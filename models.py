import re

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
        self.method = None
        self.sequence = ""

    def __str__(self):
        return "%s (%s)\n%s (%s) - %s [%s]\nComplexes: %s\nPDB: %s\n%s" % (self.sample_name, self.sample_organism, self.emdb_id, self.sample_id, self.uniprot_id, self.method, str(self.sample_complexes), str(self.pdb_ids), self.sequence)

    def get_tsv(self):
        complex_str = ';'.join([str(elem) for elem in self.sample_complexes])
        if self.method:
            return ("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.emdb_id, self.sample_id, self.sample_name, self.sample_organism, self.uniprot_id, self.method, complex_str))
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

class EMDB_complex:
    """
    EMDB complex sample obtained from the header files in the Uniprot mapping
    """

    def __init__(self, emdb_id, sample_id, complex_sample_id):
        self.emdb_id = emdb_id
        self.sample_id = sample_id
        self.complex_sample_id = complex_sample_id
        self.cpx_list = []
        self.proteins = set()
        self.method = ""
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
        self.method = ""
        self.HET  = ""
        self.lig_name = ""
        self.chembl_id = ""
        self.chebi_id = ""
        self.drugbank_id = ""

    def get_chembl_tsv(self):
        if self.chembl_id:
            return "%s\t%s\t%s\t%s\t%s\t%s\n" % (self.emdb_id, self.sample_id, self.HET, self.lig_name, self.chembl_id, self.method)
        return ""

    def get_chebi_tsv(self):
        if self.chebi_id:
            return "%s\t%s\t%s\t%s\t%s\t%s\n" % (self.emdb_id, self.sample_id, self.HET, self.lig_name, self.chebi_id, self.method)
        return ""

    def get_drugbank_tsv(self):
        if self.drugbank_id:
            return "%s\t%s\t%s\t%s\t%s\t%s\n" % (self.emdb_id, self.sample_id, self.HET, self.lig_name, self.drugbank_id, self.method)
        return ""

    def __str__(self):
        return "%s (%s) %s [%s]\nHET: %s" % (self.emdb_id, self.sample_id, self.lig_name, self.method, str(self.HET))

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
        return ("%s\t%s\t%d\t%f\n" % (self.emdb_id, self.pdb_id, self.assembly, self.molecular_weight))
