import re
import json
from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager
from datetime import date
from datetime import timedelta
import requests
import lxml.etree as ET

class DBVersion:
    """
    Database versions for EMDB_EMICSS.xml file
    """

    def __init__(self, db_list):
        self.db_list = db_list

    def execute(self):
        db_ver_list = self.db_versions(self.db_list)
        return db_ver_list

    def db_versions(self, db_list):
        today = date.today()
        offset = (today.weekday() - 2) % 7
        last_Wednesday = str(today - timedelta(days=offset))
        year = today.year
        month = '{:02d}'.format(today.month)
        year_month = f'{year}.{month}'
        week_num = today.isocalendar()[1]
        db_verison_list = []

        if "cpx" or "drugbank" in db_list:
            options = webdriver.ChromeOptions()
            options.headless = True
            driver = webdriver.Chrome(ChromeDriverManager().install(), options=options)
            if "drugbank" in db_list:
                driver.get("https://go.drugbank.com/releases/latest")
                drugbank_ver = driver.find_element_by_xpath("//*[@class='table table-bordered']/tbody//tr//td[3]").text
                db_verison_list.extend(["drugbank", drugbank_ver])
            if "cpx" in db_list:
                driver.get("http://ftp.ebi.ac.uk/pub/databases/intact/complex/current/")
                ftp_ver = driver.find_element_by_xpath("//html/body/pre").text
                cpxv = ftp_ver.split("complextab/")[1]
                cpx_ver = cpxv.split()[0]
                driver.close()
                driver.quit()
                db_verison_list.extend(["cpx", cpx_ver])
        if "pfam" in db_list:
            url = "https://pfam.xfam.org/family/Piwi/acc?output=xml"
            response = requests.get(url)
            if response.status_code == 200 and response.content:
                root = ET.fromstring(response.content)
                for x in list(root.iter('pfam')):
                    pfam_ver = x.attrib['release']
            db_verison_list.extend(["pfam", pfam_ver])
        if "interpro" in db_list:
            url = f"https://www.ebi.ac.uk/interpro/api/"
            response = requests.get(url)
            if response.status_code == 200:
                res_text = response.text
                data = json.loads(res_text)
                if 'databases' in data:
                    ipr_ver = data['databases']['interpro']['version']
            db_verison_list.extend(["interpro", ipr_ver])
        if "chembl" in db_list:
            url = f"https://www.ebi.ac.uk/chembl/api/data/status/"
            response = requests.get(url)
            if response.status_code == 200 and response.content:
                root = ET.fromstring(response.content)
                for x in list(root.iter('response')):
                    chembl_ver = x.find('chembl_db_version').text
            db_verison_list.extend(["chembl", chembl_ver])
        if "chebi" in db_list:
            db_verison_list.extend(["chebi", year_month])
        if "go" in db_list:
            go_ver = re.sub('-', '', str(last_Wednesday))
            db_verison_list.extend(["go", go_ver])
        if "uniprot" in db_list:
            db_verison_list.extend(["uniprot", year_month])
        if "pdbe" in db_list:
            pdbe_ver = f'{week_num}.{str(year)[-2:]}'
            db_verison_list.extend(["pdbe", pdbe_ver])
        if "pdbekb" in db_list:
            db_verison_list.extend(["pdbekb", year_month])
        if "alphafold" in db_list:
            db_verison_list.extend(["alphafold", "1.0"])
        if "empiar" in db_list:
            empiar_ver = re.sub('-', '', str(today))
            db_verison_list.extend(["empiar", empiar_ver])

        return db_verison_list
