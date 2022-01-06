import argparse, configparser, os
from pathlib import Path
from resources.OrcidMapping import OrcidMapping
from XMLParser import XMLParser
from glob import glob
import shutil
import logging
from joblib import Parallel, delayed
formatter = logging.Formatter('%(message)s')

def setup_logger(name, log_file, level=logging.INFO, mode='w'):
    """To setup as many loggers as you want"""

    handler = logging.FileHandler(log_file, mode=mode)        
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

def start_logger_if_necessary(log_name, log_file):
    logger = logging.getLogger(log_name)
    if len(logger.handlers) == 0:
        logger.setLevel(logging.INFO)
        fh = logging.FileHandler(log_file, mode='a')
        logger.addHandler(fh)
    return logger

def run(filename):
    id_num = filename.split('-')[1]
    print(f"Running EMD-{id_num}")
    xml_filepath = os.path.join(filename, f"header/emd-{id_num}-v30.xml")
    xml = XMLParser(xml_filepath)
    if orcid:
        orcid_log = start_logger_if_necessary("orcid_logger", orcid_log_file) if orcid else None
        orcid_mapping = OrcidMapping(xml.citations, pmc_api, emdb_pubmed)
        orcid_map = orcid_mapping.execute()
        orcid_mapping.export_tsv(orcid_log)

if __name__ == "__main__":
    ######### Command : python /Users/amudha/project/git_code/added_annotations/OrcidAnnotations.py
    # -w /Users/amudha/project/ -f /Users/amudha/project/EMD_XML/ --orcid

    prog = "EMDBAddedAnnotations"
    usage = """
            Mapping EMDB entries to ORCID ids and corresponding author name
            Example:
            python OrcidAnnotations.py -w '[{"/path/to/working/folder"}]'
            -f '[{"/path/to/EMDB/header/files/folder"}]'
            --orcid
          """

    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path .")
    parser.add_argument('-f', '--headerDir', type=Path, help="Directory path to the EMDB version 3.0 header files.")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of threads.")
    parser.add_argument("--orcid", type=bool, nargs='?', const=True, default=False, help="Mapping ORCID ID to Publications in entries ")

    args = parser.parse_args()
    orcid = args.orcid

    #Get config variables:
    config = configparser.ConfigParser()
    env_file = os.path.join(Path(__file__).parent.absolute(), "config.ini")
    config.read(env_file)
    pmc_api = config.get("api", "pmc")
    emdb_pubmed = config.get("file_paths", "emdb_pubmed")
    emdb_orcid = config.get("file_paths", "emdb_orcid")

    #Start loggers
    if orcid:
        orcid_log_file = os.path.join(args.workDir, 'emdb_orcid.log')
        orcid_log_file_backup = os.path.join(args.workDir, 'emdb_orcid_backup.log')
        orcid_log = setup_logger('orcid_logger', orcid_log_file)
        orcid_log.info("EMDB_ID\tAUTHOR_NAME\tORCID_ID\tPROVENANCE")

    Parallel(n_jobs=args.threads)(delayed(run)(file) for file in glob(os.path.join(args.headerDir, '*')))

    ######### Since emdb_orcid.log file is updated rarely and it consumes time, the file is backed up ####
    shutil.copyfile(orcid_log_file, orcid_log_file_backup)
