import unittest
import resources.EMPIARMapping
from unittest import mock

class TestEMPIARMapping(unittest.TestCase):
    """
    UnitTest for EMPIAR annotations
    """
    def setUp(self):
        super(TestEMPIARMapping, self).setUp()
        self.emdb_id = ["EMD-7456", "EMD-8959"]
        self.emp_id = [['EMPIAR-10176', 'EMPIAR-10177', 'EMPIAR-10178'], ['EMPIAR-10240']]
        self.empiar_log = "EMD-7456 EMPIAR-10176    EMDB \n EMD-7456    EMPIAR-10177    EMDB \n EMD-7456    EMPIAR-10178    EMDB \n EMD-8959    EMPIAR-10240    EMDB"


    def test_execute(self):
        mock_log = mock.MagicMock()
        mock_log.info.return_value.fetchall.return_value = self.empiar_log
        for n in range(len(self.emdb_id)):
            EMPIARMap = resources.EMPIARMapping.EMPIARMapping(self.emdb_id[n], {self.emdb_id[n]: self.emp_id[n]}, mock_log)
            for x in range(len(self.emp_id[n])):
                self.assertEqual(EMPIARMap.execute()[x].emdb_id, self.emdb_id[n])
                self.assertEqual(EMPIARMap.execute()[x].empiar_id, self.emp_id[n][x])