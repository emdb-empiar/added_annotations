import argparse
from pathlib import Path
from EMICSS.DBVersion import get_db_versions
from EMICSS.EmicssXML import EmicssXML
from EMICSS.EmicssModels import EmicssModels

if __name__ == "__main__":
    prog = "Write EMICSS annotations to XML files"
    usage = """
            EMICSS to EMICSS-XML 
            Example:
            python write_xml.py -w '[{"/path/to/working/folder"}]'
          """

    parser = argparse.ArgumentParser(prog=prog, usage=usage, add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit.")
    parser.add_argument('-w', '--workDir', type=Path, help="Main working directory path .")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of threads.")
    args = parser.parse_args()
    workDir = args.workDir

    db_list = ["pdbe", "empiar", "uniprot", "chembl", "chebi", "drugbank", "pubmed", "pubmedcentral", "issn",
               "orcid", "cpx", "go", "interpro", "pfam", "cath", "scop", "scop2", "scop2B", "pdbekb", "alphafold"]

    db_version = get_db_versions(db_list)
    emicss_models = EmicssModels(workDir)
    model_dict = emicss_models.worker()
    for emdb_id in list(model_dict.keys()):
        print(f"Writing {emdb_id} XML")
        packed_models = model_dict[emdb_id]
        write_annotation_xml = EmicssXML(emdb_id, workDir, db_version)
        write_annotation_xml.write(packed_models)






    

