from glob import glob
import lxml.etree as ET
import os, re, subprocess
import urllib.parse
import urllib.request
import requests
import logging
from fuzzywuzzy import fuzz

BLAST_DB = "/nfs/public/rw/pdbe/httpd-em/software/ncbi-blast-2.11.0+/database/uniprot_sprot"  #  uniprotkb_swissprot
BLASTP_BIN = "/nfs/public/rw/pdbe/httpd-em/software/ncbi-blast-2.11.0+/bin/blastp"

uniprot_api = "www.uniprot.org/uniprot/?query=\"%s\" AND database:(type:pdb %s)&format=tab&limit=100&columns=id,organism-id&sort=score"

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(funcName)s:%(message)s')

file_handler = logging.FileHandler('logging_uniprot.log', mode ='w')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

class Protein:
	"""
	Defines the attributes of a protein sample in a EMDB entry
	"""
	def __init__(self, emdb_id, sample_id):
		self.emdb_id = emdb_id
		self.sample_id = sample_id
		self.sample_name = ""
		self.sample_organism = None
		self.pdb_ids = []
		self.sample_complexes = []
		self.uniprot_id = None
		self.method = None
		self.sequence = ""

	def __str__(self):
		return "%s (%s)\n%s (%s) - %s [%s]\nComplexes: %s\nPDB: %s\n%s" % (self.sample_name, self.sample_organism, self.emdb_id, self.sample_id, self.uniprot_id, self.method, str(self.sample_complexes), str(self.pdb_ids), self.sequence)

	def get_tsv(self):
		complex_str = ';'.join([str(elem) for elem in self.sample_complexes])
		if self.method:
			return ("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.emdb_id, self.sample_id, self.sample_name, self.sample_organism, self.uniprot_id, self.method, complex_str))
		else:
			return ""

class UniprotMapping:
	"""
	Map EMDB protein samples to Uniprot IDs
	"""
	def __init__(self, headerDir, workDir):
		self.header_dir = headerDir
		self.output_dir = workDir
		self.uniprot_tab = os.path.join(self.output_dir, "uniprot.tsv")
		self.uniprot = {}
		self.proteins = []

	def parseUniprot(self):
		"""
		Parse the Uniprot tab file containg all entries with models and resulting in
		 a dictionary of pdb_id -> [(Uniprot_id, protein_names)]
		"""
		with open(self.uniprot_tab, 'r') as fr:
			next(fr) #Skip header
			for line in fr:
				line = line.strip()
				uniprot_id, pdb_list, protein_names = line.split('\t')
				pdb_list = pdb_list.split(';')[:-1]
				for pdb_id in pdb_list:
					pdb_id = pdb_id.lower()
					if pdb_id in self.uniprot:
						self.uniprot[pdb_id].append((uniprot_id, protein_names))
					else:
						self.uniprot[pdb_id] = [(uniprot_id, protein_names)]

	def export_tsv(self):
		filepath = os.path.join(self.output_dir, "emdb_uniprot.tsv")
		with open(filepath, 'w') as fw:
			fw.write("#EMDB_ID\tSAMPLE_ID\tSAMPLE_NAME\tNCBI_ID\tUNIPROT_ID\tMETHOD\tSAMPLE_COMPLEX_IDS\n")
			for protein in self.proteins:
				fw.write(protein.get_tsv())

	def execute(self):
		for fn in glob(os.path.join(str(self.header_dir), '*')):
			id_num = fn.split('-')[1]
			xml_filename = "emd-" + id_num + "-v30.xml"
			xml_dirpath = os.path.join(str(self.header_dir), fn, "header")
			xml_filepath = os.path.join(xml_dirpath, xml_filename)
			proteins = self.read_xml(xml_filepath)

			for protein in proteins:
				if not protein.method:
					#If contain a model: Uniprot search
					found = False
					if len(protein.pdb_ids) > 0:
						found = self.query_uniprot(protein)
					
					if not found:
						if protein.sequence:
							self.blastp(protein)
			self.proteins += proteins

	def blastp(self, protein):
		#Create fasta file
		fasta_file = os.path.join(self.output_dir, protein.emdb_id + ".fasta")
		with open(fasta_file, "w") as f:
			f.write(">seq\n%s" % protein.sequence)
		qout = os.path.join(self.output_dir, protein.emdb_id + ".xml")
		#### Using biopython for Blastp ##
		#blastp_command = NcbiblastpCommandline(query=fasta_file, db=db_path, out=qout, evalue='1e-40')
		#blastp_command()
		blastp_command = [BLASTP_BIN, "-query", fasta_file, "-db", BLAST_DB, "-out", qout, "-outfmt", "5",
						  "-evalue", "1e-40"]
		subprocess.call(blastp_command)
		if os.path.isfile(qout):
			uniprot_id = self.extract_uniprot_from_blast(qout, protein.sample_organism)
			if uniprot_id:
				protein.uniprot_id = uniprot_id
				protein.method = "BLASTP"

	def extract_uniprot_from_blast(self, out_file, ncbi_id):
		"""
		Extracting the UNIPROT ID from BLASTP XML output
		"""
		with open(out_file, 'r') as outFile:
			tree = ET.parse(outFile)
			root = tree.getroot()
			seq_len = root.find('BlastOutput_query-len').text
			for x in list(root.iter('Hit')):
				Hit_def = x.find('Hit_def').text
				Hit_split = Hit_def.split("OS")
				uniprot_id = Hit_split[0]
				tax_id = re.findall('OX=(.*)=', Hit_def, re.S)[0].split("GN")[0].strip()
				ident_num = x.find('Hit_len').text
				if tax_id == ncbi_id:
					coverage = int(ident_num) / int(seq_len)
					if coverage < 1.0:
						continue
					return uniprot_id
		return None

	def query_uniprot(self, protein):
		pdb_found = False
		for pdb_id in protein.pdb_ids:
			if pdb_id in self.uniprot:
				uniprot_list = self.uniprot[pdb_id]
				pdb_found = True
				break

		if not pdb_found:
			return False

		best_score = 0
		best_match = ""

		for uniprot_id, uniprot_names in uniprot_list:
			score = fuzz.token_set_ratio(protein.sample_name.lower(), uniprot_names.lower())
			if score > best_score:
				best_score = score
				best_match = uniprot_id

		if best_match and best_score > 80:
			protein.uniprot_id = best_match
			protein.method = "PDB+UNIPROT"
			return True
		return False

	def read_xml(self, xml_file):
		proteins = []
		pdb_ids = set()

		with open(xml_file, 'r') as filexml:
			tree = ET.parse(filexml)
			root = tree.getroot()
			a = root.attrib
			emd_id = a.get('emdb_id')
			prt_cpx = {} #Macromolecule -> Supramolecule
			for x in list(root.iter('pdb_reference')):
				pdb_id = x.find('pdb_id').text.lower()
				pdb_ids.add(pdb_id)

			if list(root.iter('complex_supramolecule')):			
				for x in list(root.iter('complex_supramolecule')):
					complex_id = x.attrib['supramolecule_id']
					for y in list(x.iter('macromolecule_id')):
						protein_id = y.text
						if protein_id in prt_cpx:
							prt_cpx[protein_id].add(complex_id)
						else:
							prt_cpx[protein_id] = set(complex_id)

			if list(root.iter('protein_or_peptide')):
				for x in list(root.iter('protein_or_peptide')):
					sample_id = x.attrib['macromolecule_id']
					protein = Protein(emd_id,sample_id)
					protein.pdb_ids = list(pdb_ids)
					protein.sample_name = x.find('name').text
					if sample_id in prt_cpx:
						protein.sample_complexes = list(prt_cpx[sample_id])

					if x.find('natural_source') is not None:
						qs = x.find('natural_source')
						if qs.find('organism') is not None:
							if 'ncbi' in qs.find('organism').attrib:
								ncbi_id = qs.find('organism').attrib['ncbi']
								protein.sample_organism = ncbi_id
					
					qs = x.find('sequence')
					if qs.find('external_references') is not None:
						if qs.find('external_references').attrib['type'] == 'UNIPROTKB':
							uniprot_id = qs.find('external_references').text
							protein.uniprot_id = uniprot_id
							protein.method = "AUTHOR"
					if qs.find('string') is not None:
						seq = qs.find('string').text
						seq = re.sub(r'\(\s*UNK\s*\)', 'X', seq)
						seq = seq.replace("\n", "")
						protein.sequence = seq
					proteins.append(protein)
		return proteins

	def download_uniprot(self):
		os.system('wget "https://www.uniprot.org/uniprot/?query=database:(type:pdb)&format=tab&limit=100000&columns=id,database(PDB),protein names&sort=score" -O %s' % self.uniprot_tab)


