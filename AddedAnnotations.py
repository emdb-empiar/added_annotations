import argparse
from pathlib import Path

import models
from resources.ComplexPortalMapping import CPMapping
from resources.ComponentsMapping import ComponentsMap
from resources.UniprotMapping import UniprotMapping
from resources.StructureMapping import StructureMapping
from resources.SampleWeight import SampleWeight
from resources.EMPIARMapping import EMPIARMapping
from resources.EmicssXML import EmicssXML
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
    # -w /Users/amudha/project/ -f /Users/amudha/project/EMD_XML/ -p /Users/amudha/project/pdbeFiles/ --CPX --model
    # --component --uniprot --weight --empiar

    prog = "EMDBAddedAnnotations"
    usage = """
            Mapping EMDB entries to Complex portal, UNIPROT, chEMBL.
            Example:
            python AddedAnnotations.py -w '[{"/path/to/working/folder"}]'
            -f '[{"/path/to/EMDB/header/files/folder"}]'
            -p '[{"/path/to/PDBe/files/folder"}]'
            --download_uniprot --uniprot --CPX --component --model --weight --empiar
          """
    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path .")
    parser.add_argument('-f', '--headerDir', type=Path, help="Directory path to the EMDB version 3.0 header files.")
    parser.add_argument('-p', '--PDBeDir', type=Path, help="Directory path to the PDBe Complex portal mapping files.")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of threads.")
    parser.add_argument("--download_uniprot", type=bool, nargs='?', const=True, default=False, help="Download uniprot tab file.")
    parser.add_argument("--uniprot", type=bool, nargs='?', const=True, default=False, help="Mapping to Complex Portal.")
    parser.add_argument("--CPX", type=bool, nargs='?', const=True, default=False, help="Mapping to Complex Portal.")
    parser.add_argument("--component", type=bool, nargs='?', const=True, default=False, help="Mapping to ChEMBL, "
                                                                                              "ChEBI and DrugBank.")
    parser.add_argument("--model", type=bool, nargs='?', const=True, default=False, help="Collect MW from PDBe.")
    parser.add_argument("--weight", type=bool, nargs='?', const=True, default=False, help="Collect sample weight from header file.")
    parser.add_argument("--empiar", type=bool, nargs='?', const=True, default=False, help="Mapping EMPIAR ID to EMDB entries")
    args = parser.parse_args()

    xml = XMLParser(args.headerDir)
    xml.execute()
    uniprot = False

    #CPX mapping requires Uniprot anotation
    if args.CPX or args.uniprot:
        uniprot = True

    if uniprot:
        unp_mapping = UniprotMapping(args.workDir, xml.proteins)
        unp_mapping.parseUniprot()
        unip_map = unp_mapping.execute(args.threads)
        unp_mapping.export_tsv()
        if args.download_uniprot:
            unp_mapping.download_uniprot()
    if args.CPX:
        cpx_mapping = CPMapping(args.workDir, unp_mapping.proteins, xml.supras)
        cpx_map = cpx_mapping.execute(args.threads)
        cpx_mapping.write_cpx_map()
    if args.component:
        che_mapping = ComponentsMap(args.workDir, xml.ligands)
        lig_map = che_mapping.execute(args.threads)
        che_mapping.write_ligands()
    if args.model:
        mw_mapping = StructureMapping(args.workDir, xml.models)
        mw_map = mw_mapping.execute(args.threads)
        mw_mapping.export_tsv()
    if args.weight:
        sw_mapping = SampleWeight(args.workDir, xml.weights)
        sw_map = sw_mapping.execute(args.threads)
    if args.empiar:
        empiar_mapping = EMPIARMapping(args.workDir, models.EMPIAR)
        empiar_map = empiar_mapping.execute()

    if args.uniprot and args.CPX and args.component and args.model and args.weight:
        write_annotation_xml = EmicssXML(args.workDir, unip_map, cpx_map, lig_map, mw_map, sw_map, empiar_map)
        write_annotation_xml.execute()
