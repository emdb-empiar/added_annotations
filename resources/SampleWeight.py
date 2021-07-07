import logging
from multiprocessing import Pool

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(funcName)s:%(message)s')
file_handler = logging.FileHandler('logging_molecularweight.log')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

class SampleWeight:
    """
    Calculating the author provided total sample weight
    """

    def __init__(self, workDir, weights):
        self.workDir = workDir
        self.weights = weights

    def execute(self, threads):
        with Pool(processes=threads) as pool:
            self.weights = pool.map(self.worker, self.weights)
        return self.weights

    def worker(self, weight):
        if weight.sup_th_weight:
            weight.kind = "supra"
            weight.method = "theoretical"
            weight.sample_weight = sum(weight.sup_th_weight)
            weight.weight_unit = weight.sup_th_unit
        elif weight.sup_exp_weight:
            weight.kind = "supra"
            weight.method = "experimental"
            weight.sample_weight = sum(weight.sup_exp_weight)
            weight.weight_unit = weight.sup_exp_unit
        elif weight.macro_th_weight:
            weight.kind = "macro"
            weight.method = "theoretical"
            weight.sample_weight = sum(weight.macro_th_weight)
            weight.weight_unit = weight.macro_th_unit
        elif weight.macro_exp_weight:
            weight.kind = "macro"
            weight.method = "experimental"
            weight.sample_weight = sum(weight.macro_exp_weight)
            weight.weight_unit = weight.macro_exp_unit
        return weight



