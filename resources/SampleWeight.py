from multiprocessing import Pool
import os

class SampleWeight:
    """
    Calculating the total weight of the sample provided by the author
    """

    def __init__(self, workDir, weights, overall_mw):
        self.workDir = workDir
        self.weights = weights
        self.entries = overall_mw

    def execute(self, threads):
        with Pool(processes=threads) as pool:
            self.weights = pool.map(self.worker, self.weights)
        return self.weights

    def worker(self, weight):
        if weight.sup_th_weight:
            weight.kind = "supra"
            weight.method = "theoretical"
            weight.sample_th_weight = sum(weight.sup_th_weight)
            weight.th_unit = weight.sup_th_unit
        if weight.sup_exp_weight:
            weight.kind = "supra"
            weight.method = "experimental"
            weight.sample_exp_weight = sum(weight.sup_exp_weight)
            weight.exp_unit = weight.sup_exp_unit
        if weight.macro_th_weight:
            weight.kind = "macro"
            weight.method = "theoretical"
            weight.sample_th_weight = sum(weight.macro_th_weight)
            weight.th_unit = weight.macro_th_unit
        if weight.macro_exp_weight:
            weight.kind = "macro"
            weight.method = "experimental"
            weight.sample_exp_weight = sum(weight.macro_exp_weight)
            weight.exp_unit = weight.macro_exp_unit
        # print(weight.__dict__)
        return weight

    def export_overall_mw(self):
        filepath = os.path.join(self.workDir, "overall_mw.tsv")
        with open(filepath, 'w') as fw:
            fw.write("EMDB_ID\tOVERALL_MW\n")
            for entry in self.entries:
                fw.write(f"{str(entry)}\n")



