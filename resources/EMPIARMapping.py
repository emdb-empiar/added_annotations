import os
import json

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
        if emdb_id in empiar_dictionary:
            for empiar_id in empiar_dictionary[emdb_id]:
                empiar_logger.info(f"{emdb_id}\t{empiar_id}\tAUTHOR")