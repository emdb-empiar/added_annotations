import os
import requests
import gzip
import csv

rfam_trans_table = {
    "23S LARGE SUBUNIT RIBOSOMAL RNA": "23S ribosomal RNA",
    "23S RIBOSOMAL RNA": "23S ribosomal RNA",
    "23S RIBOSOMAL RRNA": "23S ribosomal RNA",
    "23S RNA": "23S ribosomal RNA",
    "23S RRNA": "23S ribosomal RNA",
    "23S RibosomaL RNA from E. coli": "23S ribosomal RNA",
    "23S Ribosomal RNA": "23S ribosomal RNA",
    "23S Ribosomal RNA from E. coli": "23S ribosomal RNA",
    "23S rRNA": "23S ribosomal RNA",
    "23S rRNA (2899-MER)": "23S ribosomal RNA",
    "23S ribomosomal RNA": "23S ribosomal RNA",
    "23S ribosomal RNA": "23S ribosomal RNA",
    "23S ribosomal rna": "23S ribosomal RNA",
    "23S ribsomal RNA": "23S ribosomal RNA",
    "23s RNA": "23S ribosomal RNA",
    "50S 23S RIBOSOMAL RNA": "23S ribosomal RNA",
    "50S ribosomal RNA 23S": "23S ribosomal RNA",
    "RIBOSOMAL 23S RNA": "23S ribosomal RNA",
    "RRNA-23S RIBOSOMAL RNA": "23S ribosomal RNA",
    "rRNA-23S ribosomal RNA": "23S ribosomal RNA",
    "ribosomal 23S RNA": "23S ribosomal RNA",
    "ribosomal RNA 23S": "23S ribosomal RNA",
    "ribosome RNA 23S": "23S ribosomal RNA",
    "23S  ribosomal RNA": "23S ribosomal RNA",
    "23s ribosomal RNA": "23S ribosomal RNA",
    "16S RIBOSOMAL RNA": "16S ribosomal RNA",
    "16S RNA": "16S ribosomal RNA",
    "16S RRNA": "16S ribosomal RNA",
    "16S RRNA (E.COLI NUMBERING)": "16S ribosomal RNA",
    "16S Ribosomal RNA": "16S ribosomal RNA",
    "16S Ribosomal RNA from E. coli": "16S ribosomal RNA",
    "16S SMALL SUBUNIT RIBOSOMAL RNA": "16S ribosomal RNA",
    "16S rRNA": "16S ribosomal RNA",
    "16S rRNA (1504-MER)": "16S ribosomal RNA",
    "16S ribosomal RNA": "16S ribosomal RNA",
    "16s rRNA": "16S ribosomal RNA",
    "16s ribosomal RNA": "16S ribosomal RNA",
    "30S 16S RIBOSOMAL RNA": "16S ribosomal RNA",
    "30S 16S ribosomal RNA": "16S ribosomal RNA",
    "ribosomal RNA 16S": "16S ribosomal RNA",
    "ribosomal 16S RNA": "16S ribosomal RNA",
    "5s RNA": "5S ribosomal RNA",
    "5S ribosomal RNA (120-MER)": "5S ribosomal RNA",
    "5S ribomosomal RNA": "5S ribosomal RNA",
    "5S rRNA (119-MER)": "5S ribosomal RNA",
    "5S rRNA": "5S ribosomal RNA",
    "rRNA 5S": "5S ribosomal RNA",
    "5S Ribosomal RNA": "5S ribosomal RNA",
    "5S RRNA CHAIN OF THE LARGE RIBOSOMAL SUBUNIT": "5S ribosomal RNA",
    "5S RRNA": "5S ribosomal RNA",
    "5S RNA": "5S ribosomal RNA",
    "5S RIBOSOMAL RRNA": "5S ribosomal RNA",
    "5S RIBOSOMAL RNA": "5S ribosomal RNA",
    "5s ribosomal RNA": "5S ribosomal RNA",
    "5S LARGE SUBUNIT RIBOSOMAL RNA": "5S ribosomal RNA",
    "50S ribosomal RNA 5S": "5S ribosomal RNA",
    "50S 5S RIBOSOMAL RNA": "5S ribosomal RNA",
    "RIBOSOMAL 5S RNA": "5S ribosomal RNA",
    "RRNA-5S RIBOSOMAL RNA": "5S ribosomal RNA",
    "ribosome RNA 5S": "5S ribosomal RNA",
    "ribosomal RNA 5S": "5S ribosomal RNA",
    "ribosomal 5S RNA": "5S ribosomal RNA",
    "rRNA-5S ribosomal RNA": "5S ribosomal RNA",
    "Thermus thermophilus 5S rRNA": "5S ribosomal RNA",
    "5.8S RIBOSOMAL RNA": "5.8S ribosomal RNA",
    "5.8S RRNA": "5.8S ribosomal RNA",
    "5.8S RNA": "5.8S ribosomal RNA",
    "5.8S RRNA CHAIN OF THE LARGE RIBOSOMAL SUBUNIT": "5.8S ribosomal RNA",
    "5.8S Ribosomal RNA": "5.8S ribosomal RNA",
    "5.8S rRNA": "5.8S ribosomal RNA",
    "5.8s rRNA": "5.8S ribosomal RNA",
    "rRNA 5.8S": "5.8S ribosomal RNA",
    "18S RIBOSOMAL RNA": "18S Ribosomal RNA",
    "18S RRNA": "18S Ribosomal RNA",
    "18S RRNA 2": "18S Ribosomal RNA",
    "18S RRNA OF THE SMALL RIBOSOMAL SUBUNIT": "18S Ribosomal RNA",
    "18S rRNA": "18S Ribosomal RNA",
    "18S ribosomal RNA": "18S Ribosomal RNA",
    "18s ribosomal RNA": "18S Ribosomal RNA",
    "ribosomal RNA 18S": "18S Ribosomal RNA",
    "18S RNA": "18S Ribosomal RNA",
    "HUMAN 18S RIBOSOMAL RNA": "18S Ribosomal RNA",
    "Human 18S ribosomal RNA": "18S Ribosomal RNA",
    "28S RIBOSOMAL RNA": "28S ribosomal RNA",
    "28s ribosomal RNA": "28S ribosomal RNA",
    "28S RRNA": "28S ribosomal RNA",
    "28S Ribosomal RNA": "28S ribosomal RNA",
    "28S rRNA": "28S ribosomal RNA",
    "28S ribosomal RNA, mitochondial": "28S ribosomal RNA",
    "ALPHA CHAIN OF THE LARGE RIBOSOMAL SUBUNIT 28S RRNA": "28S ribosomal RNA",
    "BETA CHAIN OF THE LARGE RIBOSOMAL SUBUNIT 28S RRNA": "28S ribosomal RNA",
    "25S RIBOSOMAL RNA": "25S ribosomal RNA",
    "25s ribosomal RNA": "25S ribosomal RNA",
    "25S RRNA": "25S ribosomal RNA",
    "25S RNA": "25S ribosomal RNA",
    "25S rRNA": "25S ribosomal RNA",
    "25s rRNA": "25S ribosomal RNA",
    "25S Ribosomal RNA": "25S ribosomal RNA",
    "messenger RNA": "Messenger RNA",
    "MESSENGER RNA": "Messenger RNA",
    "mRNA": "Messenger RNA",
    "MRNA": "Messenger RNA",
    "mRNA fragment": "Messenger RNA",
    "Messenger RNA, mRNA": "Messenger RNA",
    "RNA message": "Messenger RNA",
    "FRAGMENT OF MESSENGER RNA": "Messenger RNA",
    "Fragment of messenger RNA": "Messenger RNA",
    "TRNA": "tRNA",
    "tRNA": "tRNA"
}

def generate_rfam_dictionary(workdir):
    ftpdir = os.path.join(workdir, "ftp_data")
    urls = [
        "https://ftp.ebi.ac.uk/pub/databases/Rfam/.preview/pdb_full_region.txt.gz",
        "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz"
    ]
    rfam_files = [
        os.path.join(ftpdir, "pdb_full_region.txt.gz"),
        os.path.join(ftpdir, "family.txt.gz")
    ]

    for url, filepath in zip(urls, rfam_files):
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(filepath, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)

    # Read the pdb_full_region file and build a dictionary
    rfam_dictionary = {}
    pdb_full_region_data = []
    family_data = {}

    with gzip.open(rfam_files[0], 'rt', encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            pdb_full_region_data.append(row)
            pdb_id = row[1]
            rfam_acc = row[0]
            if pdb_id not in rfam_dictionary:
                rfam_dictionary[pdb_id] = {}
            rfam_dictionary[pdb_id][rfam_acc] = None

    # Parse the family file and update the dictionary values
    with gzip.open(rfam_files[1], 'rt', encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            family_data[row[0]] = [row[3], row[1]]

    for pdb_id, acc_dict in rfam_dictionary.items():
        for rfam_acc in acc_dict:
            if rfam_acc in family_data:
                rfam_dictionary[pdb_id][rfam_acc] = family_data[rfam_acc]

    return rfam_dictionary


class RfamMapping:
    """
    Map Rfam IDs to EMDB ID and sample id along with the name
    """

    def __init__(self, rfam):
        self.rfam = rfam

    def execute(self, rfam_dictionary):
        for rf in self.rfam:
            rf = self.worker(rfam_dictionary, rf)
        return self.rfam

    def worker(self, rfam_dictionary, rf):
        if rf.pdb_id in rfam_dictionary:
            for rfam_acc, rfam_names in rfam_dictionary[rf.pdb_id].items():
                rfam_id = rfam_trans_table[rf.sample_name]
                rf.rfam_id = rfam_id
                if rfam_id == rfam_names[0]:
                    rf.rfam_acc = rfam_acc
                rf.provenance = "RFAM"

    def export_tsv(self, rfam_logger):
        for rf in self.rfam:
            rfam_logger.info(str(rf))