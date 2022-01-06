import requests, json
from models import Alphafold

class AlphaFoldMapping:
	"""
	Collect AlphaFoldDB entries according to the UniProtKb ID
	"""
	def __init__(self, alphafold_ftp):
		self.proteins = []
		self.alphafold_ftp = alphafold_ftp

	def execute(self, proteins):
		for protein in proteins:
			if protein.uniprot_id and protein.sample_id:
				uid = protein.uniprot_id
				############# if Alphafold API is used ##########
				# url = f"https://alphafold.ebi.ac.uk/api/prediction/{uid}"
				# try:
				# 	response = requests.get(url)
				# 	if response.status_code == 200 and response.content:
				# 		if json.loads(response.content):
				# 			alphafold = Alphafold()
				# 			alphafold.unip_id = uid
				# 			alphafold.link = f"https://alphafold.ebi.ac.uk/entry/{uid}"
				# 			alphafold.provenance = "AlphaFoldDB"
				# 			protein.alphafold.append(alphafold)
				# 			self.proteins.append(protein)
				# except:
				# 	print(f"WARNING: AlphaFold {uid} failed!")
				#########################
				with open(self.alphafold_ftp) as f:
					if uid in f.read():
						alphafold = Alphafold()
						alphafold.unip_id = uid
						alphafold.link = f"https://alphafold.ebi.ac.uk/entry/{uid}"
						alphafold.provenance = "AlphaFoldDB"
						protein.alphafold.append(alphafold)
						self.proteins.append(protein)
		return self.proteins

	def export_tsv(self, logger):
		for protein in self.proteins:
			row = f"{protein.emdb_id}\t{protein.sample_id}\t{protein.uniprot_id}\tAlphaFold DB"
			logger.info(row)