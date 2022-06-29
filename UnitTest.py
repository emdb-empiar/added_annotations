import argparse, configparser, os
from pathlib import Path
from unit_test.test_UniprotMapping import TestUniprotMapping
from unit_test.test_EMPIARMapping import TestEMPIARMapping
from unit_test.test_ComponentsMapping import TestComponentsMapping
from unit_test.test_StructureMapping import TestStructureMapping

if __name__ == "__main__":
    prog = "UnitTest(EMICSS)"
    usage = """
            UnitTest for all EMICSS
            Example:
            python UnitTest.py -w '[{"/path/to/working/folder"}]' --uniprot --empiar --component
          """

    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path .")
    parser.add_argument("--all", type=bool, nargs='?', const=True, default=False, help="Unit test for all annotations")
    parser.add_argument("--uniprot", type=bool, nargs='?', const=True, default=False, help="Unit test for UniProt annotations")
    parser.add_argument("--empiar", type=bool, nargs='?', const=True, default=False, help="Unit test for EMPIAR annotations")
    parser.add_argument("--component", type=bool, nargs='?', const=True, default=False, help="Unit test for ChEBML, ChEBI, DrugBank annotations")
    parser.add_argument("--model", type=bool, nargs='?', const=True, default=False, help="Unit test for MW from PDBe.")
    parser.add_argument('-N', default=500, help="Number of simultaneosly papers to be included in a query.")

    args = parser.parse_args()
    workDir = args.workDir
    uniprot = args.uniprot
    empiar = args.empiar
    component = args.component
    model = args.model
    if args.all:
        uniprot = True
        empiar = True
        component = True
        model = True
    N = int(args.N)

    # Get config variables:
    config = configparser.ConfigParser()
    env_file = os.path.join(Path(__file__).parent.absolute(), "config.ini")
    config.read(env_file)
    blast_db = config.get("file_paths", "BLAST_DB")
    blastp_bin = config.get("file_paths", "BLASTP_BIN")

    if uniprot:
        test_uniprot = TestUniprotMapping(workDir, blast_db, blastp_bin)
        test_uniprot.test_generate_unp_dictionary()
        test_uniprot.test_worker()
        test_uniprot.test_blastp()
        test_uniprot.test_extract_uniprot_from_blast()
    if empiar:
        test_empiar = TestEMPIARMapping()
        test_empiar.test_execute()
    if component:
        test_component = TestComponentsMapping()
        # test_component.test_parseCCD()
        test_component.test_worker()
    if model:
        test_model = TestStructureMapping()
        test_model.test_worker()
