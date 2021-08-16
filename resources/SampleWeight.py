import os

class SampleWeight:
    """
    Calculating the total weight of the sample provided by the author
    """

    def __init__(self, weights):
        self.weights = weights

    def execute(self):
        for weight in self.weights:
            weight = self.worker(weight)
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
        return weight



