import os
import json

class EMPIARMapping:
    """
    Mapping EMPIAR ID to EMDB entry
    """

    def __init__(self, workDir, EMPIAR, emdb_empiar_list):
        self.output_dir = workDir
        self.EMPIAR = EMPIAR
        self.emdb_empiar_list = emdb_empiar_list

    def execute(self):
        empiars = []
        with open(self.emdb_empiar_list, "r") as file:
            data = json.load(file)
            for key, value in data.items():
                for item in value:
                    empiar = self.EMPIAR(key)
                    empiar.emdb_id = key
                    empiar.empiar_id = item
                    empiar.provenance = "AUTHOR"
                    empiars.append(empiar)
        return empiars