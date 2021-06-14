from multiprocessing import Pool
import os
import lxml.etree as ET

#assembly_ftp = '/Users/neli/EBI/annotations/data/pdbe/assembly/'
assembly_ftp = '/Users/amudha/project/ftp_data/pdbe/assembly/'
#assembly_ftp = '/nfs/ftp/pub/databases/msd/assemblies/split/'

class StructureMapping:
	"""
	Map PDB ids and they overall molecular weight
	"""
	def __init__(self, workDir, models):
		self.output_dir = workDir
		self.models = models

	def execute(self, threads):
		with Pool(processes=threads) as pool:
			self.models = pool.map(self.worker, self.models)
		return self.models

	def worker(self, model):
		assembly_file = os.path.join(assembly_ftp, "%s/%s/%s-assembly.xml" % (model.pdb_id[1:3],model.pdb_id,model.pdb_id))
		mw, order = self.parse_assembly(assembly_file)
		if mw:
			model.molecular_weight = mw
			model.assembly = order
			return model
		return None

	def parse_assembly(self, file):
		with open(file, 'r') as fr:
			tree = ET.parse(fr)
			root = tree.getroot()
			assemblies = tree.xpath("//assembly")
			for assembly in assemblies:
				if assembly.attrib['prefered'] == 'True':
					mw = float(assembly.attrib['molecular_weight'])
					order = int(assembly.attrib['order'])
					return mw,order
			return None, None

	def export_tsv(self):
		filepath = os.path.join(self.output_dir, "emdb_model.tsv")
		with open(filepath, 'w') as fw:
			fw.write("EMDB_ID\tPDB_ID\tASSEMBLY\tMOLECULAR_WEIGHT\n")
			for model in self.models:
				fw.write(str(model))


