import os
import models

class EmicssModels:
    """
     Convert the annotations to models to write the XML files
    """

    def __init__(self, emdb_id, log_files):
        self.emdb_id = emdb_id
        self.log_files = log_files

    def execute(self):
        packed_models = self.worker()
        return packed_models

    def worker(self):
        packed_models = {}
        packed_models['EMDB_ID'] = self.emdb_id
        for file in self.log_files:
            filename = os.path.basename(file)

            if filename == "emdb_overall_mw.log":
                row = self.read_filecont(file)
                wgt = models.Weight(emdb_id=self.emdb_id, overall_mw=float(row[1]), units="MDa", provenance="EMDB")
                packed_models["WEIGHT"] = wgt

            if filename == "emdb_model.log":
                row = self.read_filecont(file)
                model = models.Model(emdb_id=self.emdb_id, pdb_id=row[1], assembly=row[2], molecular_weight=row[3])
                packed_models["MODEL"] = [model]
        return packed_models

    def read_filecont(self, file):
        with open(file, 'r') as txtin:
            next(txtin)
            for line in txtin:
                if line.startswith(self.emdb_id):
                    row = line.strip('\n').split('\t')
            return row