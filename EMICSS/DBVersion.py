import re
import json
from datetime import date
from datetime import timedelta
import requests
import lxml.etree as ET

def get_db_versions(db_list):
    """
    Database versions for EMDB_EMICSS.xml file
    """
    today = date.today()
    today_num = re.sub('-', '', str(today))
    offset = (today.weekday() - 2) % 7
    last_Wednesday = str(today - timedelta(days=offset))
    year = today.year
    month = '{:02d}'.format(today.month)
    year_month = f'{year}.{month}'
    week_num = today.isocalendar()[1]
    db_verison_list = {}

    db_verison_list['EMDB'] = None
    if "cpx" in db_list:
        url = "http://ftp.ebi.ac.uk/pub/databases/intact/complex/current/"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            html = response.content.decode('utf-8')
            cpx_ver = re.findall("\d*\-\w*\-\d*", html)[0]
            db_verison_list['Complex Portal'] = cpx_ver
    if "drugbank" in db_list:
        url = "https://go.drugbank.com/releases/latest"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            html = response.content.decode('utf-8')
            vers = re.findall("<td>\d*\.\d*\.\d*", html)[0]
            drugbank_ver = vers.split(">")[1]
            db_verison_list['DrugBank'] = drugbank_ver
    if "pfam" in db_list:
        url = "http://pfam.xfam.org/family/Piwi/acc?output=xml"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            root = ET.fromstring(response.content)
            for x in list(root.iter('pfam')):
                pfam_ver = x.attrib['release']
            db_verison_list['Pfam'] = pfam_ver
    if "interpro" in db_list:
        url = f"https://www.ebi.ac.uk/interpro/api/"
        response = requests.get(url)
        if response.status_code == 200:
            res_text = response.text
            data = json.loads(res_text)
            if 'databases' in data:
                ipr_ver = data['databases']['interpro']['version']
            db_verison_list['InterPro'] = ipr_ver
    if "cath" in db_list:
        url = "https://www.cathdb.info/"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            html = response.content.decode('utf-8')
            vers = re.findall("<h1>CATH / Gene3D <small>\w*\d*\.\d*", html)[0]
            cath_ver = vers.split("small>")[1]
            db_verison_list['CATH'] = cath_ver
    if "scop" in db_list:
        db_verison_list['SCOP'] = "1.75" # TODO: Replace hardcoded version
    if "scop2" in db_list:
        db_verison_list['SCOP2'] = "2.0" # TODO: Replace hardcoded version
    if "chembl" in db_list:
        url = f"https://www.ebi.ac.uk/chembl/api/data/status/"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            root = ET.fromstring(response.content)
            for x in list(root.iter('response')):
                chembl_ver = x.find('chembl_db_version').text
            db_verison_list['ChEMBL'] = chembl_ver
    if "chebi" in db_list:
        db_verison_list['ChEBI'] = None # year_month
    if "go" in db_list:
        url = "http://current.geneontology.org/release_stats/go-stats-summary.json"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            data = json.loads(response.content)
            db_verison_list['GO'] = data['release_date']
    if "uniprot" in db_list:
        url = "https://ftp.uniprot.org/pub/databases/uniprot/relnotes.txt"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            db_verison_list['UniProt'] = response.text.split("\n")[0]
    if "alphafold" in db_list:
        db_verison_list['AlphaFold DB'] = "2.0" #TODO: Replace hardcoded version
    if "empiar" in db_list:
        db_verison_list['EMPIAR'] = None
    if "pdbekb" in db_list:
        db_verison_list["PDBe-KB"] = None
    if "pmc" in db_list:
        db_verison_list["PubMed"] = None
        db_verison_list["PubMed Central"] = None
    return db_verison_list