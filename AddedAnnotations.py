import argparse, configparser, os
from pathlib import Path
import models
from resources.ComplexPortalMapping import CPMapping
from resources.ComponentsMapping import ComponentsMapping
from resources.UniprotMapping import UniprotMapping
from resources.StructureMapping import StructureMapping
from resources.SampleWeight import SampleWeight
from resources.EMPIARMapping import EMPIARMapping
from resources.PubmedMapping import PubmedMapping
from resources.GOMapping import GOMapping
from EMICSS.EmicssXML import EmicssXML
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

    xml = XMLParser(args.headerDir)
    xml.execute()
    mapping_list = []

    uniprot = args.uniprot
    cpx = args.CPX
    component = args.component
    model = args.model
    weight = args.weight
    empiar = args.empiar
    pmc = args.pmc
    go = args.GO
    
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
    config.read("config.ini")

    if uniprot:
        blast_db = config.get("file_paths", "BLAST_DB")
        blastp_bin = config.get("file_paths", "BLASTP_BIN")
        unp_mapping = UniprotMapping(args.workDir, xml.proteins, blast_db, blastp_bin)
        unp_mapping.parseUniprot()
        unip_map = unp_mapping.execute(args.threads)
        unp_mapping.export_tsv()
        if args.download_uniprot:
            unp_mapping.download_uniprot()
        mapping_list.extend(["UNIPROT", unip_map])
    if cpx:
        CP_ftp = config.get("file_paths", "CP_ftp")
        cpx_mapping = CPMapping(args.workDir, unp_mapping.proteins, xml.supras, CP_ftp)
        cpx_map = cpx_mapping.execute(args.threads)
        cpx_mapping.write_cpx_map()
        mapping_list.extend(["COMPLEX", cpx_map])
    if component:
        components_cif = config.get("file_paths", "components_cif")
        che_mapping = ComponentsMapping(args.workDir, xml.ligands, components_cif)
        lig_map = che_mapping.execute(args.threads)
        che_mapping.write_ligands()
        mapping_list.extend(["LIGANDS", lig_map])
    if model:
        assembly_ftp = config.get("file_paths", "assembly_ftp")
        mw_mapping = StructureMapping(args.workDir, xml.models, assembly_ftp)
        mw_map = mw_mapping.execute(args.threads)
        mw_mapping.export_tsv()
        mapping_list.extend(["MODEL", mw_map])
    if weight:
        sw_mapping = SampleWeight(args.workDir, xml.weights, xml.overall_mw)
        sw_map = sw_mapping.execute(args.threads)
        sw_mapping.export_overall_mw()
        mapping_list.extend(["WEIGHT", sw_map])
    if empiar:
        emdb_empiar_list = config.get("file_paths", "emdb_empiar_list")
        empiar_mapping = EMPIARMapping(args.workDir, models.EMPIAR)
        empiar_map = empiar_mapping.execute()
        mapping_list.extend(["EMPIAR", empiar_map])
    if pmc:
        pmc_ftp_gz = config.get("file_paths", "pmc_ftp_gz")
        pmc_ftp = config.get("file_paths", "pmc_ftp")
        pmc_mapping = PubmedMapping(args.workDir, xml.citations, pmc_ftp, pmc_ftp_gz)
        pmc_map = pmc_mapping.execute(args.threads)
        mapping_list.extend(["CITATION", pmc_map])
    if go:
        shifts_GO = config.get("file_paths", "sifts_GO")
        GO_obo = config.get("file_paths", "GO_obo")
        GO_mapping = GOMapping(args.workDir, xml.GOs, shifts_GO, GO_obo)
        GO_map = GO_mapping.execute(args.threads)
        mapping_list.extend(["GO", GO_map])

    write_annotation_xml = EmicssXML(args.workDir, mapping_list)
    write_annotation_xml.execute()
