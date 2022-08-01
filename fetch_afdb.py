import argparse, configparser, os, json
from pathlib import Path
import logging
import requests
from glob import glob
import lxml.etree as ET

def get_afdb_ids(alphafold_ftp):
    alphafold_ids = set()
    with open(alphafold_ftp) as f:
        for line in f:
            id = line.split(',')[0]
            alphafold_ids.add(id)
    return alphafold_ids

if __name__ == "__main__":
    prog = "EMDBAddedAnnotations(EMICSS)"
    usage = """
            Collect citation information be read by added annotations
            Example:
            python fetch_pubmed.py -w '[{"/path/to/working/folder"}]'
            -f '[{"/path/to/EMDB/header/files/folder"}]'
            --citation
          """

    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path .")
    parser.add_argument('-a', '--afdbFile', type=Path, help="List of AlphaFold ids.")

    work_dir = args.workDir
    uniprot_file = os.path.join(work_dir, 'emdb_uniprot.log')
    afdb_map_file = os.path.join(work_dir, 'emdb_alphafold.log')
    afdb_file = args.afdbFile

    afdb = get_afdb_ids(afdb_file)

    uniprot_reader = open(uniprotFile, 'r')
    next(uniprot_reader) #Skip header
    afdb_writer = open(afdb_map_file, 'w')
EMDB_ID SAMPLE_ID   SAMPLE_NAME SAMPLE_COPIES   NCBI_ID UNIPROT_ID  PROVENANCE  SAMPLE_COMPLEX_IDS
    for line in uniprot_reader:
        line = line.strip()
        row = line.split("\t")
        if len(row) > 0:
            emdb_id = row[0]
            sample_id = row[1]





    

