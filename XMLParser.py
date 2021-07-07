import lxml.etree as ET
from glob import glob
from models import Protein, Supra, Ligand, Model, Weight
import os, re

class XMLParser:
	"""
	Parse the XML files and store the objects
	"""

	def __init__(self, header_dir):
		self.header_dir = header_dir
		self.proteins = []
		self.supras = []
		self.ligands = []
		self.models = []
		self.weights = []

	def execute(self):
		for fn in glob(os.path.join(str(self.header_dir), '*')):
			id_num = fn.split('-')[1]
			xml_filename = "emd-" + id_num + "-v30.xml"
			xml_dirpath = os.path.join(str(self.header_dir), fn, "header")
			xml_filepath = os.path.join(xml_dirpath, xml_filename)
			proteins, supras, ligands, models , weight = self.read_xml(xml_filepath)
			self.proteins += proteins
			self.supras += supras
			self.ligands += ligands
			self.models += models
			self.weights += weight

	def read_xml(self, xml_file):
		proteins = []
		supras = []
		pdb_ids = []
		ligands = []
		weights = []

		with open(xml_file, 'r') as filexml:
			tree = ET.parse(filexml)
			root = tree.getroot()
			a = root.attrib
			emd_id = a.get('emdb_id')
			prt_cpx = {} #Macromolecule -> Supramolecule
			for x in list(root.iter('pdb_reference')):
				pdb_id = x.find('pdb_id').text.lower()
				model = Model(emd_id, pdb_id)
				pdb_ids.append(model)

			if list(root.iter('complex_supramolecule')):			
				for x in list(root.iter('complex_supramolecule')):
					complex_id = x.attrib['supramolecule_id']
					supra = Supra(emd_id, complex_id)
					supra.supra_id = "supra_" + complex_id
					supra.kind = "supra"
					if x.find('name') is not None:
						supra.supra_name = x.find('name').text
					supras.append(supra)

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
					protein.pdb = pdb_ids
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
							protein.provenance = "AUTHOR"
					if qs.find('string') is not None:
						seq = qs.find('string').text
						#seq = re.sub(r'\(\s*UNK\s*\)', 'X', seq)
						seq = re.sub(r'\(.*?\)', 'X', seq)
						seq = seq.replace("\n", "")
						protein.sequence = seq
					proteins.append(protein)

			supramolecule_list = ["cell_supramolecule", "complex_supramolecule", "organelle_or_cellular_component_supramolecule",
								  "sample_supramolecule", "virus_supramolecule"]
			for element in supramolecule_list:
				if list(root.iter(element)):
					for x in list(root.iter(element)):
						weight = Weight(emd_id)
						weight.provenance = "AUTHOR"
						if x.find('parent') is not None:
							par_child = x.find('parent').text
							weight.type = "parent"
						else:
							par_child = None
						if par_child == "0" or par_child is None:
							if x.find('molecular_weight') is not None:
								if x.find('number_of_copies') is not None:
									num_copies = x.find('number_of_copies').text
								else:
									num_copies = 1
								sup_wei = x.find('molecular_weight')
								if sup_wei.find('theoretical') is not None:
									sup_th_wei = float(sup_wei.find('theoretical').text)*float(num_copies)
									(weight.sup_th_weight).append(sup_th_wei)
									if 'units' in sup_wei.find('theoretical').attrib:
										sup_th_unit = sup_wei.find('theoretical').attrib['units']
										weight.sup_th_unit = sup_th_unit
								if sup_wei.find('experimental') is not None:
									sup_exp_wei = float(sup_wei.find('experimental').text)*float(num_copies)
									(weight.sup_exp_weight).append(sup_exp_wei)
									if 'units' in sup_wei.find('experimental').attrib:
										sup_exp_unit = sup_wei.find('experimental').attrib['units']
										weight.sup_exp_unit = sup_exp_unit
						if par_child != 0:
							if not weight.sup_th_weight and not weight.sup_exp_weight:
								macromolecule_list = ["protein_or_peptide", "ligand"]
								for item in macromolecule_list:
									if list(root.iter(item)):
										for x in list(root.iter(item)):
											if x.find('molecular_weight') is not None:
												weight.type = "child"
												if x.find('number_of_copies') is not None:
													num_copies = x.find('number_of_copies').text
												else:
													num_copies = 1
												mol_wei = x.find('molecular_weight')
												if mol_wei.find('theoretical') is not None:
													th_wei = float(mol_wei.find('theoretical').text)*float(num_copies)
													(weight.macro_th_weight).append(th_wei)
													if 'units' in mol_wei.find('theoretical').attrib:
														th_weight_unit = mol_wei.find('theoretical').attrib['units']
														weight.macro_th_unit = th_weight_unit
												if mol_wei.find('experimental') is not None:
													exp_wei = float(mol_wei.find('experimental').text)*float(num_copies)
													(weight.macro_exp_weight).append(exp_wei)
													if 'units' in mol_wei.find('experimental').attrib:
														exp_weight_unit = mol_wei.find('experimental').attrib['units']
														weight.macro_exp_unit = exp_weight_unit
						weights.append(weight)

			if list(root.iter('ligand')):
				for x in list(root.iter('ligand')):
					if x is not None:
						ligand_id = x.attrib['macromolecule_id']
						ligand = Ligand(emd_id, ligand_id)
						HET = x.find('formula')
						if HET is not None:
							ligand.HET = HET.text
						lig_name = x.find('name').text
						if lig_name:
							ligand.lig_name = lig_name
						lig_copies = x.find('number_of_copies').text
						if lig_copies:
							ligand.lig_copies = lig_copies
						else:
							ligand.lig_copies = "1"

						for ref in x.iter('external_references'):
							if ref.attrib['type'] == 'CHEMBL':
								chembl_id = ref.text
								ligand.chembl_id = chembl_id
								ligand.provenance = "AUTHOR"
							if ref.attrib['type'] == 'CHEBI':
								chebi_id = ref.text
								ligand.chebi_id = chebi_id
								ligand.provenance = "AUTHOR"
							if ref.attrib['type'] == 'DRUGBANK':
								drugbank_id = ref.text
								ligand.drugbank_id = drugbank_id
								ligand.provenance = "AUTHOR"
					ligands.append(ligand)
		return proteins, supras, ligands, pdb_ids, weights