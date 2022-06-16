import json
from models import Empiar

def generate_emp_dictionary(emdb_empiar_list):
    empiar_dictionary = {}
    with open(emdb_empiar_list, "r") as file:
        data = json.load(file)
        for emdb, empiar_list in data.items():
            empiar_dictionary[emdb] = empiar_list
    return empiar_dictionary

class EMPIARMapping:
    """
    Mapping EMPIAR ID to EMDB entry
    """

    def __init__(self, emdb_id, empiar_dictionary, empiar_logger):
        self.empiar = Empiar(emdb_id)
        self.emdb_id = emdb_id
        self.empiar_dictionary = empiar_dictionary
        self.empiar_logger = empiar_logger

    def execute(self):
        empiars = []
        if self.emdb_id in self.empiar_dictionary:
            for empiar_id in self.empiar_dictionary[self.emdb_id]:
                self.empiar_logger.info(f"{self.emdb_id}\t{empiar_id}\tEMDB")
                empiar = Empiar(self.emdb_id)
                empiar.emdb_id = self.emdb_id
                empiar.empiar_id = empiar_id
                empiars.append(empiar)
        return empiars