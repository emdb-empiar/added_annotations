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

    if "cpx" in db_list:
        url = "http://ftp.ebi.ac.uk/pub/databases/intact/complex/current/"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            html = response.content.decode('utf-8')
            cpx_ver = re.findall("\d*\-\w*\-\d*", html)[0]
            db_verison_list['cpx'] =  cpx_ver
    if "drugbank" in db_list:
        url = "https://go.drugbank.com/releases/latest"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            html = response.content.decode('utf-8')
            vers = re.findall("<td>\d*\.\d*\.\d*", html)[0]
            drugbank_ver = vers.split(">")[1]
            db_verison_list['drugbank'] = drugbank_ver
    if "pfam" in db_list:
        url = "http://pfam.xfam.org/family/Piwi/acc?output=xml"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            root = ET.fromstring(response.content)
            for x in list(root.iter('pfam')):
                pfam_ver = x.attrib['release']
            db_verison_list['pfam'] = pfam_ver
    if "interpro" in db_list:
        url = f"https://www.ebi.ac.uk/interpro/api/"
        response = requests.get(url)
        if response.status_code == 200:
            res_text = response.text
            data = json.loads(res_text)
            if 'databases' in data:
                ipr_ver = data['databases']['interpro']['version']
            db_verison_list['interpro'] = ipr_ver
    if "cath" in db_list:
        url = "https://www.cathdb.info/"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            html = response.content.decode('utf-8')
            vers = re.findall("<h1>CATH / Gene3D <small>\w*\d*\.\d*", html)[0]
            cath_ver = vers.split("small>")[1]
            db_verison_list['cath'] = cath_ver
    if "scop" in db_list:
        db_verison_list['scop'] = "1.75"
    if "scop2" in db_list:
        db_verison_list['scop2'] = "2.0"
    if "chembl" in db_list:
        url = f"https://www.ebi.ac.uk/chembl/api/data/status/"
        response = requests.get(url)
        if response.status_code == 200 and response.content:
            root = ET.fromstring(response.content)
            for x in list(root.iter('response')):
                chembl_ver = x.find('chembl_db_version').text
            db_verison_list['chembl'] = chembl_ver
    if "chebi" in db_list:
        db_verison_list['chebi'] = year_month
    if "go" in db_list:
        go_ver = re.sub('-', '', str(last_Wednesday))
        db_verison_list['go'] = go_ver
    if "uniprot" in db_list:
        db_verison_list['uniprot'] = f'{year}_{month}'
    if "alphafold" in db_list:
        db_verison_list['alphafold'] = "2.0" #TODO: Get a dynamic version
    if "empiar" in db_list:
        db_verison_list['empiar'] = today_num
    return db_verison_list