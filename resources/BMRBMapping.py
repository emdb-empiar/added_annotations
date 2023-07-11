from models import BMRB

BMRB_PDB = "https://bmrb.io/ftp/pub/bmrb/nmr_pdb_integrated_data/adit_nmr_matched_pdb_bmrb_entry_ids.csv"
BMRB_UNIPROT = "http://api.bmrb.io/v2/mappings/bmrb/uniprot?format=text&match_type=all"

def generate_bmrb_dictionary(bmrb_pdb_file):
    bmrb_dictionary = {}
    bmrb_pdb = open(bmrb_pdb_file, 'r')
    for line in bmrb_pdb:
        line = line.strip()
        row = line.split(",")
        if len(row) > 0:
            bmrb_id = row[0]
            pdb_id = row[1]
            if bmrb_id in bmrb_dictionary:
                bmrb_dictionary[bmrb_id].append(pdb_id.lower())
            else:
                bmrb_dictionary[bmrb_id] = [pdb_id.lower()]
    return bmrb_dictionary

def generate_bmrb_uni_dictionary(bmrb_uniprot_file):
    bmrb_uni_dictionary = {}
    bmrb_uni = open(bmrb_uniprot_file, 'r')
    for line in bmrb_uni:
        line = line.strip()
        row = line.split()
        if len(row) > 0:
            bmrb_id = row[0]
            uni_id = row[1]
            bmrb_uni_dictionary[bmrb_id] = uni_id
    return bmrb_uni_dictionary

class BMRBMapping:
    """
    Mapping BMRB ID to EMDB entry
    """

    def __init__(self, emdb_id, bmrb_dictionary, bmrb_logger):
        self.bmrb = BMRB(emdb_id)
        self.emdb_id = emdb_id
        self.bmrb_dictionary = bmrb_dictionary
        self.bmrb_logger = bmrb_logger

    def execute(self, proteins):
        pdb_id = []
        for protein in proteins:
            uniprot_id = protein.uniprot_id
            for pdb in protein.pdb:
                pdb = pdb.pdb_id
                if pdb != None:
                        pdb_id.append(pdb) if pdb not in pdb_id else pdb_id
            for key, value in self.bmrb_dictionary.items():
                if value == pdb_id:
                    bmrb_id = key
                    pdb_ids = (", ".join(pdb_id))
                    sample_id = protein.sample_id
                    self.bmrb_logger.info(f"{self.emdb_id}\t{sample_id}\t{bmrb_id}\t{pdb_ids}\t{uniprot_id}\tBMRB")
                    bmrb = BMRB()
                    bmrb.bmrb_id = bmrb_id
                    bmrb.pdb_id = pdb_id
                    bmrb.unp_id = uniprot_id
                    bmrb.provenance = "BMRB"
                    protein.bmrb = bmrb
        return proteins