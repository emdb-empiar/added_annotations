class SampleWeight:
    """
    Calculating the total weight of the sample provided by the author
    """

    def __init__(self, weights, overall_weight):
        self.weights = weights
        self.overall_weight = overall_weight

    def execute(self):
        for weight in self.weights:
            weight = self.worker(weight)
        return self.weights

    def worker(self, weight):
        weight.overall_mw = self.overall_weight
        weight.units = "MDa"
        weight.provenance = "EMDB"

        return weight



