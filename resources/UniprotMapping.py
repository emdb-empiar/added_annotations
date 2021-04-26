from glob import glob
import lxml.etree as ET
import os, re, subprocess
import urllib.parse
import urllib.request
import requests

sifts_file = "/Users/neli/EBI/annotations/uniprot_pdb.csv"
#sifts_file = r'/nfs/ftp/pub/databases/msd/sifts/csv/uniprot_pdb.csv'
BLAST_DB = "uniprotkb_swissprot"  #  uniprot_sprot

uniprot_api = "www.uniprot.org/uniprot/?query=\"%s\" AND database:(type:pdb %s)&format=tab&limit=100&columns=id,organism-id&sort=score"

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
		if self.method:
			return ("%s\t%s\t%s\t%s\t%s\t%s\n" % (self.emdb_id, self.sample_id, self.sample_name, self.sample_organism, self.uniprot_id, self.method))
		else:
			return ""
			
class UniprotMapping:
	"""
	Map EMDB protein samples to Uniprot IDs
	"""
	def __init__(self, headerDir, workDir):
		self.header_dir = headerDir
		self.output_dir = workDir
		self.sifts = {}
		self.proteins = []

		#parse Sifts file
		self.parseSifts()


	def parseSifts(self):
		"""
		Parse the PDB-Uniprot sifts file resulting in a dictionary of pdb_id -> [Uniprot_id]
		"""
		with open(sifts_file, 'r') as fr:
			for line in fr:
				line = line.strip()
				if line[0] != '#':
					uniprot_id, pdb_list = line.split(',')
					pdb_list = pdb_list.split(';')
					for pdb_id in pdb_list:
						if pdb_id in self.sifts:
							self.sifts[pdb_id].add(uniprot_id)
						else:
							self.sifts[pdb_id] = set([uniprot_id])

	def export_tsv(self):
		filepath = os.path.join(self.output_dir, "emdb_uniprot.tsv")
		with open(filepath, 'w') as fw:
			fw.write("#EMDB_ID\tSAMPLE_ID\tSAMPLE_NAME\tNCBI_ID\tUNIPROT_ID\tMETHOD\n")
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
					#If contains model: SIFTS + Uniprot search API
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
		db_path = os.path.join(self.output_dir, BLAST_DB)  #
		qout = os.path.join(self.output_dir, protein.emdb_id + ".xml")
		#### Using biopython for Blastp ##
		#blastp_command = NcbiblastpCommandline(query=fasta_file, db=db_path, out=qout, evalue='1e-40')
		#blastp_command()
		blastp_command = ["blastp", "-query", fasta_file, "-db", db_path, "-out", qout, "-outfmt", "5",
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
		params = {
			'query': '%s AND database:(type:pdb %s)' % (protein.sample_name, protein.pdb_ids[0]),
			'limit': 100,
			'format': 'tab',
			'columns': 'id,organism-id,database(PDB)',
			'sort': 'score'
			}
		data = urllib.parse.urlencode(params)
		data = data.encode('utf-8')

		#url = uniprot_api % (protein.sample_name, protein.pdb_ids[0])
		#url = url.replace(" ","%20")
		#url = "https://" + url
		url = "https://www.uniprot.org/uniprot/"
		req = urllib.request.Request(url,data)
		with urllib.request.urlopen(req) as f:
			for line in f:
				line = line.decode('utf-8')
				line = line.strip()
				
				unp_id, ncbi_id, pdbs = line.split('\t')

				if protein.sample_organism:
					if ncbi_id == protein.sample_organism:
						protein.uniprot_id = unp_id
						protein.method = "PDB+UNIPROT"
						return True
				else:
					protein.uniprot_id = unp_id
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


