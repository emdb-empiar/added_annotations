import argparse, os
from pathlib import Path
from glob import glob
from EMICSS.DBVersion import get_db_versions
from EMICSS.EmicssXML import EmicssXML
from EMICSS.EmicssModels import EmicssModels
from joblib import Parallel, delayed

def run(filename):
    id_num = filename.split('-')[1]
    emdb_id = f"EMD-{id_num}"
    print(f"Writing EMD-{id_num} XML")

    # log_files = glob(os.path.join(workDir, 'emdb_*.log'))
    log_files = ['/Users/amudha/project/emdb_overall_mw.log', '/Users/amudha/project/emdb_model.log']
    emicss_models = EmicssModels(emdb_id, log_files)
    packed_models = emicss_models.execute()
    print(packed_models)
    write_annotation_xml = EmicssXML(args.workDir, db_version)
    write_annotation_xml.write(packed_models)

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
    parser.add_argument('-f', '--headerDir', type=Path, help="Directory path to the EMDB version 3.0 header files.")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of threads.")
    args = parser.parse_args()
    workDir = args.workDir
    headerDir = args.headerDir

    db_list = ["pdbe", "empiar", "uniprot", "chembl", "chebi", "drugbank", "pubmed", "pubmedcentral", "issn",
               "orcid", "cpx", "go", "interpro", "pfam", "cath", "scop", "scop2", "scop2B", "pdbekb", "alphafold"]

    db_version = get_db_versions(db_list)
    files = glob(os.path.join(args.headerDir, '*'))

    Parallel(n_jobs=args.threads)(delayed(run)(file) for file in files)





    

