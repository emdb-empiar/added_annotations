import unittest
import resources.EMPIARMapping
from unittest import mock

class TestEMPIARMapping(unittest.TestCase):
    """
    UnitTest for EMPIAR annotations
    """
    def test_generate_emp_dictionary(self):
        emp_list = {"EMD-7456": ["EMPIAR-10176", "EMPIAR-10177", "EMPIAR-10178"], "EMD-8959": ["EMPIAR-10240"]}
        emp_dict = {'EMD-7456': ['EMPIAR-10176', 'EMPIAR-10177', 'EMPIAR-10178'],'EMD-8959': ['EMPIAR-10240']}
        self.assertEqual(resources.EMPIARMapping.generate_emp_dictionary(emp_list), emp_dict)

    def test_execute(self):
        mock_log = mock.MagicMock()
        emdb_id = ["EMD-7456", "EMD-8959"]
        emp_id = [['EMPIAR-10176', 'EMPIAR-10177', 'EMPIAR-10178'], ['EMPIAR-10240']]
        empiar_log = "EMD-7456	EMPIAR-10176	EMDB \n EMD-7456	EMPIAR-10177	EMDB \n EMD-7456	EMPIAR-10178	EMDB \n EMD-8959	EMPIAR-10240	EMDB"
        mock_log.info.return_value.fetchall.return_value = empiar_log
        for n in range(len(emdb_id)):
            EMPIARMap = resources.EMPIARMapping.EMPIARMapping(emdb_id[n], {emdb_id[n]: emp_id[n]}, mock_log)
            for x in range(len(emp_id[n])):
                self.assertEqual(EMPIARMap.execute()[x].emdb_id, emdb_id[n])
                self.assertEqual(EMPIARMap.execute()[x].empiar_id, emp_id[n][x])