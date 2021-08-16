import argparse, configparser, os, sys, time
from pathlib import Path
import models
from resources.ComplexPortalMapping import CPMapping
from resources.ComponentsMapping import ComponentsMapping
from resources.UniprotMapping import UniprotMapping, generate_unp_dictionary, download_uniprot
from resources.StructureMapping import StructureMapping
from resources.SampleWeight import SampleWeight
from resources.EMPIARMapping import EMPIARMapping
from resources.PubmedMapping import PubmedMapping
from resources.GOMapping import GOMapping
from EMICSS.EmicssXML import EmicssXML
from XMLParser import XMLParser
from glob import glob
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
    if uniprot:
        uniprot_log = start_logger_if_necessary("uniprot_logger", uniprot_log_file)
        unp_mapping = UniprotMapping(args.workDir, xml.proteins, uniprot_tab, blast_db, blastp_bin)
        unip_map = unp_mapping.execute()
        unp_mapping.export_tsv(uniprot_log)
    if cpx:
        cpx_mapping = CPMapping(args.workDir, unp_mapping.proteins, xml.supras, CP_ftp)
        cpx_map = cpx_mapping.execute(args.threads)
        cpx_mapping.write_cpx_map()

"""
List of things to do:
  - Change the GO terms to be obtained from the UniProt instead of PMC
  - Adapt the unit tests to work with this version
"""

if __name__ == "__main__":
    ######### Command : python /Users/amudha/project/ComplexPortal/AddedAnnotations.py
    # -w /Users/amudha/project/ -f /Users/amudha/project/EMD_XML/ -p /Users/amudha/project/pdbeFiles/ --CPX --model
    # --component --uniprot --weight --empiar --pmc

    prog = "EMDBAddedAnnotations"
    usage = """
            Mapping EMDB entries to Complex portal, UNIPROT, chEMBL.
            Example:
            python AddedAnnotations.py -w '[{"/path/to/working/folder"}]'
            -f '[{"/path/to/EMDB/header/files/folder"}]'
            
            -p '[{"/path/to/PDBe/files/folder"}]'
            --download_uniprot --uniprot --CPX --component --model --weight --empiar --pmc --GO
          """
    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path .")
    parser.add_argument('-f', '--headerDir', type=Path, help="Directory path to the EMDB version 3.0 header files.")
    parser.add_argument('-p', '--PDBeDir', type=Path, help="Directory path to the PDBe Complex portal mapping files.")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of threads.")
    parser.add_argument("--all", type=bool, nargs='?', const=True, default=False, help="Fetch all external resources.")
    parser.add_argument("--download_uniprot", type=bool, nargs='?', const=True, default=False, help="Download uniprot tab file.")
    parser.add_argument("--uniprot", type=bool, nargs='?', const=True, default=False, help="Mapping to Complex Portal.")
    parser.add_argument("--CPX", type=bool, nargs='?', const=True, default=False, help="Mapping to Complex Portal.")
    parser.add_argument("--component", type=bool, nargs='?', const=True, default=False, help="Mapping to ChEMBL, ChEBI and DrugBank.")
    parser.add_argument("--model", type=bool, nargs='?', const=True, default=False, help="Collect MW from PDBe.")
    parser.add_argument("--weight", type=bool, nargs='?', const=True, default=False, help="Collect sample weight from header file.")
    parser.add_argument("--empiar", type=bool, nargs='?', const=True, default=False, help="Mapping EMPIAR ID to EMDB entries")
    parser.add_argument("--pmc", type=bool, nargs='?', const=True, default=False, help="Mapping publication ID to EMDB entries")
    parser.add_argument("--GO", type=bool, nargs='?', const=True, default=False, help="Mapping GO ids to EMDB entries")
    args = parser.parse_args()

    mapping_list = []

    uniprot = args.uniprot
    cpx = args.CPX
    component = args.component
    model = args.model
    weight = args.weight
    empiar = args.empiar
    pmc = args.pmc
    go = args.GO
    uniprot_dictionary = {}
    
    #CPX mapping requires Uniprot anotation
    if cpx:
        uniprot = True

    if args.all:
        uniprot = True
        cpx = True
        component = True
        model = True
        weight = True
        empiar = True
        pmc = True
        go = True

    #Get config variables:
    config = configparser.ConfigParser()
    env_file = os.path.join(Path(__file__).parent.absolute(), "config.ini")
    config.read(env_file)
    blast_db = config.get("file_paths", "BLAST_DB")
    blastp_bin = config.get("file_paths", "BLASTP_BIN")
    CP_ftp = config.get("file_paths", "CP_ftp")

    uniprot_tab = os.path.join(args.workDir, "uniprot.tsv")

    #Start loggers
    uniprot_log_file = os.path.join(args.workDir, 'emdb_uniprot.log')
    uniprot_log = setup_logger('uniprot_logger', uniprot_log_file)
    uniprot_log.info("EMDB_ID\tSAMPLE_ID\tSAMPLE_NAME\tSAMPLE_COPIES\tNCBI_ID\tUNIPROT_ID\tPROVENANCE\tSAMPLE_COMPLEX_IDS")

    if args.download_uniprot:
            download_uniprot(uniprot_tab)
    if uniprot:
        uniprot_dictionary = generate_unp_dictionary(uniprot_tab)

    Parallel(n_jobs=4)(delayed(run)(file) for file in glob(os.path.join(args.headerDir, '*')))

