import argparse, configparser, os
from pathlib import Path
from resources.CitationMapping import CitationMapping
import logging
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

if __name__ == "__main__":
    ######### Command : python /Users/amudha/project/git_code/added_annotations/CitationAnnotations.py
    # -w /Users/amudha/project/ -f /Users/amudha/project/EMD_XML/ --citation

    prog = "EMDBAddedAnnotations(EMICSS)"
    usage = """
            Mapping EMDB entries to ORCID ids and corresponding author name
            Example:
            python CitationAnnotations.py -w '[{"/path/to/working/folder"}]'
            -f '[{"/path/to/EMDB/header/files/folder"}]'
            --citation
          """

    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path .")
    parser.add_argument('-f', '--headerDir', type=Path, help="Directory path to the EMDB version 3.0 header files.")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of threads.")
    parser.add_argument("--citation", type=bool, nargs='?', const=True, default=False, help="Mapping Pubmed,Pubmed Central, "
                                                                                            "DOI & ORCID ID to Publications in entries ")

    args = parser.parse_args()
    workDir = args.workDir
    citation = args.citation
    headerDir = args.headerDir

    #Get config variables:
    config = configparser.ConfigParser()
    env_file = os.path.join(Path(__file__).parent.absolute(), "config.ini")
    config.read(env_file)
    pmc_api = config.get("api", "pmc")

    config = configparser.ConfigParser()
    env_file = os.path.join(Path(__file__).parent.absolute(), "config.ini")
    config.read(env_file)
    pmc_api = config.get("api", "pmc")

    if citation:
        pubmed_log_file = os.path.join(args.workDir, 'EPMC_pubmed.log')
        pubmed_log = setup_logger('pubmed_logger', pubmed_log_file)
        pubmed_log.info("EMDB_ID\tPUBMED_ID\tPUBMEDCENTRAL_ID\tISSN\tDOI\tPROVENANCE")
        orcid_log_file = os.path.join(args.workDir, 'EPMC_orcid.log')
        orcid_log = setup_logger('orcid_logger', orcid_log_file)
        orcid_log.info("EMDB_ID\tAUTHOR_NAME\tORCID_ID\tAUTHOR_ORDER\tPROVENANCE")
        pubmed_log = start_logger_if_necessary("pubmed_logger", pubmed_log_file) if citation else None
        orcid_log = start_logger_if_necessary("orcid_logger", orcid_log_file) if citation else None
        citation_mapping = CitationMapping(headerDir, pmc_api, pubmed_log, orcid_log)
        citation_map = citation_mapping.execute()
