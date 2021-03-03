#!/usr/bin/env/python

import unittest
import ComplexPortalMapping as CPM
import os

base_path = r'/Users/amudha/project/'

class TestComplexPortalMapping(unittest.TestCase):
    def setUp(self):
        print("setUp")
        self.filename = os.path.join(base_path + "EMD_XML/EMD-1831/header/emd-1831-v30.xml")
        self.pdbefile = os.path.join(base_path + "pdbeFiles/complex_portal_output_complete_complexes.csv")

    def test_extracting_IDs(self):
        extracted_IDs = ['EMD-1831', 'EMDB_ID', '2xzb', 'PDB_ID', 'complex_id:CPX-2195', 'PDBe',
                         'P19156', 'UNIPROT', 'P05027', 'UNIPROT']
        self.assertListEqual(CPM.CPMapping.extracting_IDs(self, self.filename, self.pdbefile), extracted_IDs)

    def test_extract_uniprot_from_blast(self):
        self.fileOut_all = os.path.join(base_path + "extract_uniprot_all.txt")
        self.fs = os.path.join(base_path + "out/EMD-7981/2_EMD-7981.txt")
        extracted_uniprot = ['Q14573', 'BLAST']
        self.assertListEqual(CPM.CPMapping.extract_uniprot_from_blast(self, self.fileOut_all,
                                                                      self.fs), extracted_uniprot)

    def test_quering_IDs(self):
        self.mapFile = os.path.join(base_path + "EMDB_CPX_map.tsv")
        self.query = ['EMD-1831', 'EMDB_ID', '2xzb', 'PDB_ID', 'complex_id:CPX-2195', 'PDBe',
                         'P19156', 'UNIPROT', 'P05027', 'UNIPROT']
        complex_IDs = ['CPX-2195', 'CPX-57']
        complex_name = ['Hydrogen:potassium-exchanging ATPase complex', 'Sodium:potassium-exchanging ATPase complex']
        self.assertEqual(CPM.CPMapping.quering_IDs(self, self.query, self.mapFile), (complex_IDs, complex_name))

if __name__ == '__main__':
    unittest.main()
