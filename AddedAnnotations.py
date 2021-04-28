import os, sys, argparse
from pathlib import Path
from resources.ComplexPortalMapping import CPMapping
from resources.ComponentsMapping import ComponentsMap
from resources.UniprotMapping import UniprotMapping
from XMLParser import XMLParser

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
            --download_uniprot --CPX --components 
          """
    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path .")
    parser.add_argument('-f', '--headerDir', type=Path, help="Directory path to the EMDB version 3.0 header files.")
    parser.add_argument('-p', '--PDBeDir', type=Path, help="Directory path to the PDBe Complex portal mapping files.")
    parser.add_argument("--download_uniprot", type=bool, nargs='?', const=True, default=False, help="Download uniprot tab file.")
    parser.add_argument("--CPX", type=bool, nargs='?', const=True, default=False, help="Mapping to Complex Portal.")
    parser.add_argument("--components", type=bool, nargs='?', const=True, default=False, help="Mapping to chEMBL, "
                                                                                              "ChEBI and DrugBank.")
    args = parser.parse_args()

    xml = XMLParser(args.headerDir)
    xml.execute()

    #Uniprot is the first and mandatory since is the base for map many other resources
    unp_mapping = UniprotMapping(args.workDir, xml.proteins)
    if args.download_uniprot:
      unp_mapping.download_uniprot()
    unp_mapping.parseUniprot()
    unp_mapping.execute()
    unp_mapping.export_tsv()
    if args.CPX:
        cpx_mapping = CPMapping(args.workDir, xml.proteins)
        cpx_mapping.execute()
        cpx_mapping.write_cpx_map()
    if args.components:
        che_mapping = ComponentsMap(args.workDir, args.headerDir)
        che_mapping.execute_annotations()