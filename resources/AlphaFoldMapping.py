import requests, json

class AlphaFoldMapping:
	"""
	Collect AlphaFoldDB entries according to the UniProtKb ID
	"""
	def __init__(self):
		self.proteins = []

	def execute(self, proteins):
		for protein in proteins:
			if protein.uniprot_id:
				uid = protein.uniprot_id
				url = f"https://alphafold.ebi.ac.uk/api/prediction/{uid}"
				response = requests.get(url)
				if response.status_code == 200 and response.content:
					if json.loads(response.content):
						self.proteins.append(protein)

		return self.proteins

	def export_tsv(self, logger):
		for protein in self.proteins:
			row = f"{protein.emdb_id}\t{protein.sample_id}\t{protein.uniprot_id}\tAlphaFoldDB"
			logger.info(row)