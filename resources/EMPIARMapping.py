import os
import json

# empiar_ftp = r'/nfs/ftp/pub/databases/emtest/'
empiar_ftp = r'/users/amudha/project/ftp_data/EMPIAR/'

class EMPIARMapping:
    """
    Mapping EMPIAR ID to EMDB entry
    """

    def __init__(self, workDir, EMPIAR):
        self.output_dir = workDir
        self.EMPIAR = EMPIAR

    def execute(self):
        empiars = []
        json_file = os.path.join(str(empiar_ftp), "emdb_empiar_list.json")
        with open(json_file, "r") as file:
            data = json.load(file)
            for key, value in data.items():
                for item in value:
                    empiar = self.EMPIAR(key)
                    empiar.emdb_id = key
                    empiar.empiar_id = item
                    empiar.provenance = "AUTHOR"
                    empiars.append(empiar)
        return empiars