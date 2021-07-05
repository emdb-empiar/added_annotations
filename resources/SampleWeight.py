import os
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

    def __init__(self, workDir, samples):
        self.workDir = workDir
        self.samples = samples
        self.given_weight = {}
        #
        # for sample in samples:
        #     if



    def execute(self, threads):
        self.given_weight = self.calculate_total_weight

        with Pool(processes=threads) as pool:
            self.given_weight = pool.map(self.worker, self.samples)
        return self.samples

    def worker(self, proteins):

        if protein.mol_th_weight is not None:


        return total_weight