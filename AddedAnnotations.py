import argparse, configparser, os, sys, time
from pathlib import Path
import models
from resources.ComplexPortalMapping import CPMapping
from resources.ComponentsMapping import ComponentsMapping, parseCCD
from resources.UniprotMapping import UniprotMapping, generate_unp_dictionary, download_uniprot
from resources.StructureMapping import StructureMapping
from resources.SampleWeight import SampleWeight
from resources.EMPIARMapping import EMPIARMapping, generate_emp_dictionary
from resources.PubmedMapping import PubmedMapping
from resources.ProteinTermsMapping import ProteinTermsMapping
from resources.PdbeKbMapping import PdbeKbMapping
from resources.AlphaFoldMapping import AlphaFoldMapping
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
        unp_mapping = UniprotMapping(args.workDir, xml.proteins, uniprot_dictionary, blast_db, blastp_bin)
        unip_map = unp_mapping.execute()
        unp_mapping.export_tsv(uniprot_log)
        mapping_list.extend(["UNIPROT", unip_map])
    if cpx:
        cpx_logger = start_logger_if_necessary("cpx_logger", cpx_log_file)
        cpx_mapping = CPMapping(unp_mapping.proteins, xml.supras, CP_ftp)
        cpx_map = cpx_mapping.execute()
        cpx_mapping.export_tsv(cpx_logger)
        mapping_list.extend(["COMPLEX", cpx_map])
    if component:
        chembl_log = start_logger_if_necessary("chembl_logger", chembl_log_file)
        chebi_log = start_logger_if_necessary("chebi_logger", chebi_log_file)
        drugbank_log = start_logger_if_necessary("drugbank_logger", drugbank_log_file)
        comp_mapping = ComponentsMapping(xml.ligands)
        comp_map = comp_mapping.execute(chembl_map, chebi_map, drugbank_map)
        comp_mapping.export_tsv(chembl_log, chebi_log, drugbank_log)
        mapping_list.extend(["LIGANDS", comp_map])
    if model:
        model_logger = start_logger_if_necessary("model_logger", model_log_file)
        mw_mapping = StructureMapping(xml.models, assembly_ftp)
        mw_map = mw_mapping.execute()
        mw_mapping.export_tsv(model_logger)
        mapping_list.extend(["MODEL", mw_map])
    if weight:
        weight_logger = start_logger_if_necessary("weight_logger", weight_log_file)
        weight_logger.info(f"{xml.emdb_id}\t{xml.overall_mw}")
        sw_mapping = SampleWeight(xml.weights)
        sw_map = sw_mapping.execute()
        mapping_list.extend(["WEIGHT", sw_map])
    if empiar:
        empiar_logger = start_logger_if_necessary("empiar_logger", empiar_log_file)
        empiar_mapping = EMPIARMapping(xml.emdb_id, empiar_dictionary, empiar_logger)
        empiar_map = empiar_mapping.execute()
        mapping_list.extend(["EMPIAR", empiar_map])
    if pmc:
        pmc_mapping = PubmedMapping(xml.citations, pmc_api)
        pmc_map = pmc_mapping.execute()
        mapping_list.extend(["CITATION", pmc_map])
    if go or interpro or pfam:
        go_log = start_logger_if_necessary("go_logger", go_log_file) if go else None
        interpro_log = start_logger_if_necessary("interpro_logger", interpro_log_file)  if interpro else None
        pfam_log = start_logger_if_necessary("pfam_logger", pfam_log_file)  if pfam else None
        PT_mapping = ProteinTermsMapping(unp_mapping.proteins, go, interpro, pfam)
        proteins_map = PT_mapping.execute()
        PT_mapping.export_tsv(go_log, interpro_log, pfam_log)
        mapping_list.extend(["PROTEIN-TERMS", proteins_map])
    if pdbekb:
        pdbekb_log = start_logger_if_necessary("pdbekb_logger", pdbekb_log_file)
        pdbekb_map = PdbeKbMapping()
        pdbekb_entries = pdbekb_map.execute(unp_mapping.proteins)
        pdbekb_map.export_tsv(pdbekb_log)
        mapping_list.extend(["PDBeKB", pdbekb_entries])
    if alphafold:
        alphafold_log = start_logger_if_necessary("alphafold_logger", alphafold_log_file)
        af_mapping = AlphaFoldMapping()
        af_entries = af_mapping.execute(unp_mapping.proteins)
        af_mapping.export_tsv(alphafold_log)
        mapping_list.extend(["ALPHAFOLD", af_entries])
    if emicss:
        # emicss_log = start_logger_if_necessary("emicss_logger", emicss_log_file)
        write_annotation_xml = EmicssXML(args.workDir, mapping_list)
        write_annotation_xml.execute()

"""
List of things to do:
  - Adapt the unit tests to work with this version
"""

if __name__ == "__main__":
    ######### Command : python /Users/amudha/project/git_code/added_annotations/AddedAnnotations.py
    # -w /Users/amudha/project/ -f /Users/amudha/project/EMD_XML/ -p /Users/amudha/project/pdbeFiles/ --CPX --model
    # --component --uniprot --weight --empiar --pmc

    prog = "EMDBAddedAnnotations"
    usage = """
            Mapping EMDB entries to Complex portal, UNIPROT, chEMBL.
            Example:
            python AddedAnnotations.py -w '[{"/path/to/working/folder"}]'
            -f '[{"/path/to/EMDB/header/files/folder"}]'
            
            -p '[{"/path/to/PDBe/files/folder"}]'
            --download_uniprot --uniprot --CPX --component --model --weight --empiar --pmc --GO --interpro --pfam --emicss
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
    parser.add_argument("--interpro", type=bool, nargs='?', const=True, default=False, help="Mapping InterPro ids to EMDB entries")
    parser.add_argument("--pfam", type=bool, nargs='?', const=True, default=False, help="Mapping pfam ids to EMDB entries")
    parser.add_argument("--pdbekb", type=bool, nargs='?', const=True, default=False, help="Mapping PDBeKB links to EMDB entries")
    parser.add_argument("--alphafold", type=bool, nargs='?', const=True, default=False, help="Mapping Alphafold links to EMDB entries")
    parser.add_argument("--emicss", type=bool, nargs='?', const=True, default=False, help="writting EMICSS XML file for each EMDB entry")
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
    interpro = args.interpro
    pfam = args.pfam
    pdbekb = args.pdbekb
    alphafold = args.alphafold
    emicss = args.emicss
    uniprot_dictionary = {}
    
    #CPX, GO, Interpro, Pfam mapping requires Uniprot anotation
    if cpx:
        uniprot = True
    if go:
        uniprot = True
    if interpro:
        uniprot = True
    if pfam:
        uniprot = True
    if pdbekb:
        uniprot = True
    if alphafold:
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
        interpro = True
        pfam = True
        pdbekb = True
        alphafold = True
        emicss = True

    #Get config variables:
    config = configparser.ConfigParser()
    env_file = os.path.join(Path(__file__).parent.absolute(), "config.ini")
    config.read(env_file)
    blast_db = config.get("file_paths", "BLAST_DB")
    blastp_bin = config.get("file_paths", "BLASTP_BIN")
    CP_ftp = config.get("file_paths", "CP_ftp")
    components_cif = config.get("file_paths", "components_cif")
    assembly_ftp = config.get("file_paths", "assembly_ftp")
    emdb_empiar_list = config.get("file_paths", "emdb_empiar_list")
    pmc_api = config.get("api", "pmc")
    uniprot_tab = os.path.join(args.workDir, "uniprot.tsv")
    GO_obo = config.get("file_paths", "GO_obo")

    #Start loggers
    if uniprot:
        uniprot_log_file = os.path.join(args.workDir, 'emdb_uniprot.log')
        uniprot_log = setup_logger('uniprot_logger', uniprot_log_file)
        uniprot_log.info("EMDB_ID\tSAMPLE_ID\tSAMPLE_NAME\tSAMPLE_COPIES\tNCBI_ID\tUNIPROT_ID\tPROVENANCE\tSAMPLE_COMPLEX_IDS")
    if cpx:
        cpx_log_file = os.path.join(args.workDir, 'emdb_cpx.log')
        cpx_log = setup_logger('cpx_logger', cpx_log_file)
        cpx_log.info("EMDB_ID\tEMDB_SAMPLE_ID\tSAMPLE_NAME\tSAMPLE_COPIES\tCPX_ID\tCPX_TITLE\tPROVENANCE\tSCORE")
    if component:
        chembl_log_file = os.path.join(args.workDir, 'emdb_chembl.log')
        chebi_log_file = os.path.join(args.workDir, 'emdb_chebi.log')
        drugbank_log_file = os.path.join(args.workDir, 'emdb_drugbank.log')
        chembl_log = setup_logger('chembl_logger', chembl_log_file)
        chembl_log.info("EMDB_ID\tSAMPLE_ID\tHET_CODE\tCOMP_NAME\tCOMP_COPIES\tChEMBL_ID\tPROVENANCE")
        chebi_log = setup_logger('chebi_logger', chebi_log_file)
        chebi_log.info("EMDB_ID\tSAMPLE_ID\tHET_CODE\tCOMP_NAME\tCOMP_COPIES\tChEBI_ID\tPROVENANCE")
        drugbank_log = setup_logger('drugbank_logger', drugbank_log_file)
        drugbank_log.info("EMDB_ID\tSAMPLE_ID\tHET_CODE\tCOMP_NAME\tCOMP_COPIES\tDRUGBANK_ID\tPROVENANCE")
    if model:
        model_log_file = os.path.join(args.workDir, 'emdb_model.log')
        model_log = setup_logger('model_logger', model_log_file)
        model_log.info("EMDB_ID\tPDB_ID\tASSEMBLY\tMOLECULAR_WEIGHT")
    if weight:
        weight_log_file = os.path.join(args.workDir, 'overall_mw.log')
        weight_log = setup_logger('weight_logger', weight_log_file)
        weight_log.info("EMDB_ID\tOVERALL_MW")
    if empiar:
        empiar_log_file = os.path.join(args.workDir, 'emdb_empiar.log')
        empiar_log = setup_logger('empiar_logger', empiar_log_file)
        empiar_log.info("EMDB_ID\tEMPIAR_ID\tPROVENANCE")
    if go:
        go_log_file = os.path.join(args.workDir, 'emdb_go.log')
        go_log = setup_logger('go_logger', go_log_file)
        go_log.info("EMDB_ID\tEMDB_SAMPLE_ID\tGO_ID\tGO_NAMESPACE\tGO_TYPE\tPROVENANCE")
    if interpro:
        interpro_log_file = os.path.join(args.workDir, 'emdb_interpro.log')
        interpro_log = setup_logger('interpro_logger', interpro_log_file)
        interpro_log.info("EMDB_ID\tEMDB_SAMPLE_ID\tINTERPRO_ID\tINTERPRO_NAMESPACE\tPROVENANCE")
    if pfam:
        pfam_log_file = os.path.join(args.workDir, 'emdb_pfam.log')
        pfam_log = setup_logger('pfam_logger', pfam_log_file)
        pfam_log.info("EMDB_ID\tEMDB_SAMPLE_ID\tPFAM_ID\tPFAM_NAMESPACE\tPROVENANCE")
    if pdbekb:
        pdbekb_log_file = os.path.join(args.workDir, 'emdb_pdbekb.log')
        pdbekb_log = setup_logger('pdbekb_logger', pdbekb_log_file)
        pdbekb_log.info("EMDB_ID\tEMDB_SAMPLE_ID\tPDBeKB_ID\tPROVENANCE")
    if alphafold:
        alphafold_log_file = os.path.join(args.workDir, 'emdb_alphafold.log')
        alphafold_log = setup_logger('alphafold_logger', alphafold_log_file)
        alphafold_log.info("EMDB_ID\tEMDB_SAMPLE_ID\tALPHAFOLDDB_ID\tPROVENANCE")
    if emicss:
        emicss_log_file = os.path.join(args.workDir, 'emdb_emicss.log')
        emicss_log = setup_logger('emicss_logger', emicss_log_file)

    if args.download_uniprot:
            download_uniprot(uniprot_tab)
    if uniprot:
        uniprot_dictionary = generate_unp_dictionary(uniprot_tab)
    if empiar:
        empiar_dictionary = generate_emp_dictionary(emdb_empiar_list)
    if component:
        chembl_map, chebi_map, drugbank_map = parseCCD(components_cif)

    Parallel(n_jobs=args.threads)(delayed(run)(file) for file in glob(os.path.join(args.headerDir, '*')))
