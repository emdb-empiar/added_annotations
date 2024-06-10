import argparse
import configparser
import os
from pathlib import Path
import csv
import requests

def rfam_mapping(ftp_file, rfam_log_file):
    with open(rfam_log_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter='\t')
        csv_writer.writerow(["EMDB_ID", "RFAM_ACCESSION", "RFAM_ID", "PROVENANCE"])

        # Open ftp_file and process lines
        with open(ftp_file, 'r') as f:
            next(f)  # Skip the header
            reader = csv.reader(f)
            for row in reader:
                emdb_ids = row[0].strip('"').split(',')
                pdb_id = row[1]
                rfam_acc = row[2]

                for emdb_id in emdb_ids:
                    # Use RFAM API to get rfam_id
                    rfam_api_url = f"https://rfam.org/family/{rfam_acc}/id"
                    response = requests.get(rfam_api_url)
                    if response.status_code == 200:
                        rfam_id = response.text.strip()
                    else:
                        rfam_id = "N/A"

                    csv_writer.writerow([emdb_id, rfam_acc, rfam_id, "PDBe"])

if __name__ == "__main__":
    prog = "RFAM (EMICSS)"
    usage = """
            EMICSS script for RFAM MAPPINGS
            Example:
            python fetch_rfam.py -w '[{"/path/to/working/folder"}]' -f '[{"/path/to/pdb_rfam/file"}]'
          """

    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path.")
    parser.add_argument('-f', '--ftpFile', type=Path, help="Path to the PDB-RFAM mapping file.")
    args = parser.parse_args()

    work_dir = args.workDir
    ftp_file = args.ftpFile
    rfam_log_file = os.path.join(work_dir, 'emdb_rfam.log')

    # Get config variables:
    config = configparser.ConfigParser()
    env_file = os.path.join(Path(__file__).parent.absolute(), "config.ini")
    config.read(env_file)

    rfam_mapping(ftp_file, rfam_log_file)