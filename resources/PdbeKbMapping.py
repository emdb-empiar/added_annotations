import requests

class PdbeKbMapping:
	"""
	Collect PDBe-KB entries according to the UniProtKb ID
	"""
	def __init__(self):
		self.proteins = []

	def execute(self, proteins):
		for protein in proteins:
			if protein.uniprot_id:
				uid = protein.uniprot_id
				url = f"https://www.uniprot.org/uniprot/?query=id:{uid}%20database:(type:pdb)&sort=score&columns=id&format=tab"
				response = requests.get(url)
				if response.status_code == 200 and response.content:
					self.proteins.append(protein)

		return self.proteins

	def export_tsv(self, logger):
		for protein in self.proteins:
			row = f"{protein.emdb_id}\t{protein.sample_id}\t{protein.uniprot_id}\tUniProtKb"
			logger.info(row)