import os, argparse
from pathlib import Path
import logging
import ComplexPortalMapping
import ChEMBLMapping

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(funcName)s:%(message)s')

file_handler = logging.FileHandler('logging_annotations.log', mode ='w')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

"""
List of things to do:
  - Add multi threading (Add try except to avoid unexopect closing)
  - Add config files to set up the paths (so we can run locally and in the cluster without having to change the code)
  - Generalize this code to work with any type of annotation instead of just complex portal and uniprot
  - Adapt the unit tests to work with this version
"""

if __name__ == "__main__":
    ######### Command : python /Users/amudha/project/ComplexPortal/AddedAnnotations.py
    # -w /Users/amudha/project/ -f /Users/amudha/project/EMD_XML/ -p /Users/amudha/project/pdbeFiles/ --CPX --chEMBL

    prog = "EMDBAddedAnnotations"
    usage = """
            Mapping EMDB entries to Complex portal, UNIPROT, chEMBL.
            Example:
            python AddedAnnotations.py -w '[{"/path/to/working/folder"}]'
            -f '[{"/path/to/EMDB/header/files/folder"}]'
            -p '[{"/path/to/PDBe/files/folder"}]'
            --CPX --chEMBL 
          """
    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path .")
    parser.add_argument('-f', '--headerDir', type=Path, help="Directory path to the EMDB version 3.0 header files.")
    parser.add_argument('-p', '--PDBeDir', type=Path, help="Directory path to the PDBe Complex portal mapping files.")
    parser.add_argument("--CPX", type=bool, nargs='?', const=True, default=False, help="Mapping to Complex Portal.")
    parser.add_argument("--chEMBL", type=bool, nargs='?', const=True, default=False, help="Mapping to chEMBL.")
    args = parser.parse_args()
    if args.CPX:
        cpx_mapping = ComplexPortalMapping.CPMapping(args.workDir, args.headerDir, args.PDBeDir)
        cpx_mapping.execute()
        cpx_mapping.write_cpx_map()
        cpx_mapping.write_uniprot_map()
    if args.chEMBL:
        che_mapping = ChEMBLMapping.ChEMBLMap(args.workDir, args.headerDir)
        che_mapping.execute_annotations()