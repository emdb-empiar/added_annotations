import argparse, configparser, os
from pathlib import Path
from unit_test.test_UniprotMapping import TestUniprotMapping
from unit_test.test_ComplexPortalMapping import TestComplexPortalMapping
from unit_test.test_EMPIARMapping import TestEMPIARMapping
from unit_test.test_ComponentsMapping import TestComponentsMapping
from unit_test.test_StructureMapping import TestStructureMapping
from unit_test.test_PublicationMapping import TestPublicationMapping
from unit_test.test_ProteinTermsMapping import TestProteinTermsMapping

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
    parser.add_argument("--cpx", type=bool, nargs='?', const=True, default=False, help="Unit test for ComplexPortal annotations")
    parser.add_argument("--empiar", type=bool, nargs='?', const=True, default=False, help="Unit test for EMPIAR annotations")
    parser.add_argument("--component", type=bool, nargs='?', const=True, default=False, help="Unit test for ChEBML, ChEBI, DrugBank annotations")
    parser.add_argument("--model", type=bool, nargs='?', const=True, default=False, help="Unit test for MW from PDBe.")
    parser.add_argument("--pmc", type=bool, nargs='?', const=True, default=False, help="Unit test for publication IDs")
    parser.add_argument("--orcid", type=bool, nargs='?', const=True, default=False, help="Unit test for ORCID IDs")
    parser.add_argument("--go", type=bool, nargs='?', const=True, default=False, help="Unit test for GO ids to EMDB entries")
    parser.add_argument("--interpro", type=bool, nargs='?', const=True, default=False, help="Unit test for InterPro annotations")
    parser.add_argument("--pfam", type=bool, nargs='?', const=True, default=False, help="Unit test for pfam annotations")
    parser.add_argument("--cath", type=bool, nargs='?', const=True, default=False, help="Unit test for Cath annotations")
    parser.add_argument("--scop", type=bool, nargs='?', const=True, default=False, help="Unit test for SCOP annotations")
    parser.add_argument("--scop2", type=bool, nargs='?', const=True, default=False, help="Unit test for SCOP2 annotations")
    parser.add_argument("--scop2B", type=bool, nargs='?', const=True, default=False, help="Unit test for SCOP2B annotations")
    parser.add_argument("--pdbekb", type=bool, nargs='?', const=True, default=False, help="Unit test for PDBeKB annotation")
    parser.add_argument("--alphafold", type=bool, nargs='?', const=True, default=False, help="Unit test for Alphafold annotation")
    parser.add_argument("--xml", type=bool, nargs='?', const=True, default=False, help="Unit test for EMICSS XML")
    parser.add_argument('-N', default=500, help="Number of simultaneosly papers to be included in a query.")

    args = parser.parse_args()
    workDir = args.workDir
    uniprot = args.uniprot
    cpx = args.cpx
    empiar = args.empiar
    component = args.component
    model = args.model
    pmc = args.pmc
    orcid = args.orcid
    go = args.go
    interpro = args.interpro
    pfam = args.pfam
    cath = args.cath
    scop = args.scop
    scop2 = args.scop2
    scop2B = args.scop2B
    pdbekb = args.pdbekb
    alphafold = args.alphafold
    xml = args.xml

    if args.all:
        uniprot = True
        cpx = True
        empiar = True
        component = True
        model = True
        pmc = True
        orcid = True
        go = True
        interpro = True
        pfam = True
        cath = True
        scop = True
        scop2 = True
        scop2B = True
        pdbekb = True
        alphafold = True
        xml = True
    N = int(args.N)

    # Get config variables:
    config = configparser.ConfigParser()
    env_file = os.path.join(Path(__file__).parent.absolute(), "config.ini")
    config.read(env_file)
    blast_db = config.get("file_paths", "BLAST_DB")
    blastp_bin = config.get("file_paths", "BLASTP_BIN")
    sifts_path = config.get("file_paths", "sifts")

    if uniprot:
        test_uniprot = TestUniprotMapping(workDir, blast_db, blastp_bin)
        test_uniprot.test_generate_unp_dictionary()
        test_uniprot.test_worker()
        test_uniprot.test_blastp()
        test_uniprot.test_extract_uniprot_from_blast()
    if cpx:
        test_cpx = TestComplexPortalMapping()
        test_cpx.test_worker()
    if empiar:
        test_empiar = TestEMPIARMapping()
        test_empiar.test_execute()
    if component:
        test_component = TestComponentsMapping()
        test_component.test_worker()
    if model:
        test_model = TestStructureMapping()
        test_model.test_worker()
    if pmc or orcid:
        test_publication = TestPublicationMapping()
        test_publication.test_generate_pubmed_dictionary()
        test_publication.test_worker()
    if go or interpro or pfam or cath or scop or scop2 or scop2B or pdbekb or alphafold:
        test_proteinterms = TestProteinTermsMapping(sifts_path, go, interpro, pfam, cath, scop, scop2, scop2B, pdbekb, alphafold)
        test_proteinterms.test_execute()