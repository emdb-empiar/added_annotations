import lxml.etree as ET
from models import Protein, Supra, Ligand, Model, Weight, Citation, GO, Sample, Interpro, Pfam, Author
import re

class XMLParser:
	"""
	Parse the XML files and store the objects
	"""

	def __init__(self, file):
		self.xml_file = file
		self.emdb_id = ""
		self.proteins = []
		self.supras = []
		self.ligands = []
		self.models = []
		self.weights = []
		self.citation = None
		self.overall_mw = 0.0
		self.read_xml()

	def get_mw(self, sample):
		if sample.xpath('.//molecular_weight/experimental'):
			return float(sample.find('molecular_weight/experimental').text)
		elif sample.xpath('.//molecular_weight/theoretical'):
			return float(sample.find('molecular_weight/theoretical').text)
		return None

	def get_n_copies(self, sample):
		try:
			number_copies = int(sample.find('number_of_copies').text)
		except AttributeError:
			number_copies = 1

		return number_copies

	def get_multiplier(self, node):
		multiplier = node.copies
		stack = node.parent

		while(len(stack) > 0):
			current_node = stack.pop()
			multiplier *= current_node.copies if current_node.copies else 1
			for parent in current_node.parent:
				stack.append(parent)
		return multiplier

	def sum_mw(self, samples, start_nodes):
		stack = []
		single_mw_list = []
		overall_mw = 0.0
		nodes_counted = set()

		for start_node_id in start_nodes:
			start_node = samples[start_node_id]
			if start_node.mw:
				if (start_node_id not in nodes_counted):
					multiplier = self.get_multiplier(start_node)
					try:
						overall_mw += (start_node.mw*multiplier)
					except TypeError:
						overall_mw += 0.0
					nodes_counted.add(start_node_id)
			else:
				stack.append(start_node)

		while(len(stack) > 0):
			current_node = stack.pop()

			if current_node.mw:
				if (current_node.id not in nodes_counted):
					multiplier = self.get_multiplier(current_node)
					try:
						overall_mw += (current_node.mw*multiplier)
					except TypeError:
						overall_mw += 0.0
					nodes_counted.add(current_node.id)
			else:
				for child in current_node.children:
					child_obj = samples[child.id]
					stack.append(child)
		return overall_mw

	def read_xml(self):
		with open(self.xml_file, 'r') as filexml:
			tree = ET.parse(filexml)
			root = tree.getroot()
			a = root.attrib
			self.emdb_id = a.get('emdb_id')
			prt_cpx = {} #Macromolecule -> Supramolecule
			for x in list(root.iter('pdb_reference')):
				pdb_id = x.find('pdb_id').text.lower()
				model = Model(self.emdb_id, pdb_id)
				self.models.append(model)

			if list(root.iter('complex_supramolecule')):			
				for x in list(root.iter('complex_supramolecule')):
					complex_id = x.attrib['supramolecule_id']
					supra = Supra(self.emdb_id, complex_id)
					supra.supra_id = "supra_" + complex_id
					supra.kind = "supra"
					if x.find('name') is not None:
						sup_name = x.find('name').text
						supra.supra_name = sup_name + "_" + complex_id
					self.supras.append(supra)

					for y in list(x.iter('macromolecule_id')):
						protein_id = y.text
						if protein_id in prt_cpx:
							prt_cpx[protein_id].add(complex_id)
						else:
							prt_cpx[protein_id] = set(complex_id)

			if list(root.iter('protein_or_peptide')):
				for x in list(root.iter('protein_or_peptide')):
					sample_id = x.attrib['macromolecule_id']
					protein = Protein(self.emdb_id,sample_id)
					protein.pdb = self.models
					protein.sample_name = x.find('name').text
					if sample_id in prt_cpx:
						protein.sample_complexes = list(prt_cpx[sample_id])
					if x.find('number_of_copies') is not None:
						protein.sample_copies = x.find('number_of_copies').text
					else:
						protein.sample_copies = "1"

					if x.find('natural_source') is not None:
						nat_sor = x.find('natural_source')
						if nat_sor.find('organism') is not None:
							if 'ncbi' in nat_sor.find('organism').attrib:
								ncbi_id = nat_sor.find('organism').attrib['ncbi']
								protein.sample_organism = ncbi_id

					qs = x.find('sequence')
					if qs.find('external_references') is not None:
						if qs.find('external_references').attrib['type'] == 'UNIPROTKB':
							uniprot_id = qs.find('external_references').text
							protein.uniprot_id = uniprot_id
							protein.provenance = "EMDB"
							for t in list(qs.iter('external_references')):
								if t.attrib['type'] == 'GO':
									go = GO()
									go.add_from_author(t.text, uniprot_id)
									if go.id and go.namespace and go.type:
										protein.go.add(go)
								elif t.attrib['type'] == 'INTERPRO':
									ipr = Interpro()
									ipr.add_from_author(t.text, uniprot_id)
									if ipr.id and ipr.namespace:
										protein.interpro.add(ipr)
								elif t.attrib['type'] == 'PFAM':
									pfam = Pfam()
									pfam.add_from_author(t.text, uniprot_id)
									if pfam.id:
										protein.pfam.add(pfam)
					if qs.find('string') is not None:
						seq = qs.find('string').text
						#seq = re.sub(r'\(\s*UNK\s*\)', 'X', seq)
						seq = re.sub(r'\(.*?\)', 'X', seq)
						seq = seq.replace("\n", "")
						protein.sequence = seq
					self.proteins.append(protein)

			supramolecule_list = ["cell_supramolecule", "complex_supramolecule", "organelle_or_cellular_component_supramolecule",
								  "sample_supramolecule", "virus_supramolecule"]
			for element in supramolecule_list:
				if list(root.iter(element)):
					for x in list(root.iter(element)):
						weight = Weight(self.emdb_id)
						weight.provenance = "EMDB"
						# t.macro_exp_unit = exp_weight_unit
						self.weights.append(weight)

			if list(root.iter('ligand')):
				for x in list(root.iter('ligand')):
					if x is not None:
						ligand_id = x.attrib['macromolecule_id']
						ligand = Ligand(self.emdb_id, ligand_id)
						HET = x.find('formula')
						if HET is not None:
							ligand.HET = HET.text
							lig_name = x.find('name').text
							if lig_name:
								ligand.lig_name = lig_name
							if x.find('number_of_copies') is not None:
								lig_copies = x.find('number_of_copies').text
								if lig_copies:
									ligand.lig_copies = lig_copies
								else:
									ligand.lig_copies = "1"

							for ref in x.iter('external_references'):
								if ref.attrib['type'] == 'CHEMBL':
									chembl_id = ref.text
									ligand.chembl_id = chembl_id
									ligand.provenance_chembl = "EMDB"
								if ref.attrib['type'] == 'CHEBI':
									chebi_id = ref.text
									ligand.chebi_id = chebi_id
									ligand.provenance_chebi = "EMDB"
								if ref.attrib['type'] == 'DRUGBANK':
									drugbank_id = ref.text
									ligand.drugbank_id = drugbank_id
									ligand.provenance_drugbank = "EMDB"
							self.ligands.append(ligand)

			if list(root.iter('primary_citation')):
				for y in list(root.iter('primary_citation')):
					citation = Citation(self.emdb_id)
					pub = y.find('journal_citation')
					for auth in y.iter('author'):
						author = Author(auth.text, int(auth.attrib['order']))
						if 'ORCID' in auth.attrib:
							author.orcid = auth.attrib['ORCID']
						citation.authors.append(author)
					nas = pub.find('title').text
					title = nas.split('\n\n', 1)[0]
					citation.title = title
					pubStatus = pub.attrib['published']
					if pubStatus == 'true':
						citation.status = "published"
					if pubStatus is None:
						continue
					for child in pub:
						pmedid = child.text
						pmedty = (child.attrib).get('type')
						if pmedty is not None:
							if pmedty == 'PUBMED':
								citation.pmedid = pmedid
								citation.provenance_pm = "EMDB"
							if pmedty == 'DOI':
								doi = pmedid.split(":")[1]
								citation.doi = doi
								citation.provenance_doi = "EMDB"
							if pmedty == 'ISSN':
								citation.issn = pmedid
								citation.provenance_issn = "EMDB"
					self.citation = citation

			#MW calculation
			sample_dic = {}
			start_nodes = set()
			macromolecules = set()
			sample = root.find('sample')
			supramolecule_list = root.xpath('.//sample/supramolecule_list/*')
			macromolecule_list = root.xpath('.//sample/macromolecule_list/*')

			for sample in macromolecule_list:
				sample_id = 'm' + sample.attrib['macromolecule_id']
				
				number_copies = self.get_n_copies(sample)
				mw = self.get_mw(sample)

				sample_obj = Sample(sample_id,mw,number_copies)
				sample_dic[sample_id] = sample_obj

				if mw:
					macromolecules.add(sample_id)

			for sample in supramolecule_list:
				sample_id = 's' + sample.attrib['supramolecule_id']
				try:
					parent = sample.find('parent').text
				except AttributeError:
					parent = '0'

				number_copies = self.get_n_copies(sample)
				mw = self.get_mw(sample)

				sample_obj = Sample(sample_id,mw,number_copies)
				
				if parent == '0':
					start_nodes.add(sample_id)
				else:
					parent_id = 's' + parent
					if parent_id in sample_dic:
						parent_obj = sample_dic[parent_id]
						sample_obj.add_parent(parent_obj)
						parent_obj.add_child(sample_obj)

				if sample.xpath('.//macromolecule_list'):
					child_list = sample.xpath('.//macromolecule_list/macromolecule/macromolecule_id/text()')
					for child in child_list:
						child_id = "m" + child
						if child_id in sample_dic:
							child_obj = sample_dic[child_id]
							child_obj.add_parent(sample_obj)
							sample_obj.add_child(child_obj)
			
				sample_dic[sample_id] = sample_obj

			for molecule_id in macromolecules:
				molecule = sample_dic[molecule_id]
				if len(molecule.parent) == 0:
					start_nodes.add(molecule_id)

			self.overall_mw = self.sum_mw(sample_dic, start_nodes)
