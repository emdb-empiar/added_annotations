import lxml.etree as ET
from glob import glob
from models import Protein, Supra, Ligand, Model
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

	def execute(self):
		for fn in glob(os.path.join(str(self.header_dir), '*')):
			id_num = fn.split('-')[1]
			xml_filename = "emd-" + id_num + "-v30.xml"
			xml_dirpath = os.path.join(str(self.header_dir), fn, "header")
			xml_filepath = os.path.join(xml_dirpath, xml_filename)
			proteins, supras, ligands, models = self.read_xml(xml_filepath)
			self.proteins += proteins
			self.supras += supras
			self.ligands += ligands
			self.models += models

	def read_xml(self, xml_file):
		proteins = []
		supras = []
		pdb_ids = []
		ligands = []

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
					supra.kind = "supra"
					if x.find('name') is not None:
						supra.supra_name = x.find('name').text
					if x.find('parent') is not None:
						par_child = x.find('parent').text
						if par_child == "0":
							supra.type = "parent"
						else:
							supra.type = "child"
					if x.find('molecular_weight') is not None:
						mol_wei = x.find('molecular_weight')
						if mol_wei.find('theoretical') is not None:
							th_wei = mol_wei.find('theoretical').text
							supra.sup_th_weight = th_wei
							if 'units' in mol_wei.find('theoretical').attrib:
								th_weight_unit = mol_wei.find('theoretical').attrib['units']
								supra.sup_th_unit = th_weight_unit
						if mol_wei.find('experimental') is not None:
							exp_wei = mol_wei.find('experimental').text
							supra.sup_exp_weight = exp_wei
							if 'units' in mol_wei.find('experimental').attrib:
								exp_weight_unit = mol_wei.find('experimental').attrib['units']
								supra.sup_exp_unit = exp_weight_unit
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

					if x.find('molecular_weight') is not None:
						mol_wei = x.find('molecular_weight')
						if mol_wei.find('theoretical') is not None:
							th_wei = mol_wei.find('theoretical').text
							protein.mol_th_weight = th_wei
							if 'units' in mol_wei.find('theoretical').attrib:
								th_weight_unit = mol_wei.find('theoretical').attrib['units']
								protein.mol_th_unit = th_weight_unit
						if mol_wei.find('experimental') is not None:
							exp_wei = mol_wei.find('experimental').text
							protein.mol_exp_weight = exp_wei
							if 'units' in mol_wei.find('experimental').attrib:
								exp_weight_unit = mol_wei.find('experimental').attrib['units']
								protein.mol_exp_unit = exp_weight_unit
					
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
		return proteins, supras, ligands, pdb_ids