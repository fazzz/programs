import argparse as arg
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
import chromedriver_binary

parser = arg.ArgumentParser(description='Download PDBefold msa file.')
parser.add_argument('pdbid', metavar='str', type=str, help='pdbid')
parser.add_argument('cid', metavar='str', type=str, help='chainid')
args = parser.parse_args()

pdbid = args.pdbid
cid = args.cid

driver = webdriver.Chrome()
url = "http://www.ebi.ac.uk/msd-srv/ssm/cgi-bin/ssmserver?q=" + pdbid + ",c=" + cid
driver.get(url)
WebDriverWait(driver, 60).until(EC.presence_of_element_located((By.CLASS_NAME,"title-line")))
time.sleep(10)
driver.find_element_by_name('download_algns').click()
html = driver.page_source
msa = pdbid + "_" + cid + ".msa"
with open(msa, 'w', encoding='utf-8') as f:
    f.write(html)
driver.quit()
