import unittest
import resources.ComponentsMapping
from unit_test.setters import set_ligand

class TestComponentsMapping(unittest.TestCase):
    """
    UnitTest for Ligand annotations
    """
    def setUp(self):
        super(TestComponentsMapping, self).setUp()
        self.ligands = [
            set_ligand("EMD-8959", "2", "CA", "CALCIUM ION", "4", None, "29108", "DB14577", None, "PDBe-CCD", "PDBe-CCD"),
            set_ligand("EMD-8959", "3", "D12", "DODECANE", "8", "CHEMBL30959", "28817", None, "PDBe-CCD", "PDBe-CCD", None),
            set_ligand("EMD-8959", "4", "D10", "DECANE", "2", "CHEMBL134537", None, None, "PDBe-CCD", None, "PDBe-CCD")]
        self.chembl_map = {'D12': 'CHEMBL30959', 'D10': 'CHEMBL134537'}
        self.chebi_map = {'D12': '28817', 'CA': '29108'}
        self.drugbank_map = {'CA': 'DB14577'}


    def test_worker(self):
        LigandMap = resources.ComponentsMapping.ComponentsMapping(self.ligands)
        for n in range(len(self.ligands)):
            self.assertEqual(LigandMap.worker(self.ligands[n], self.chembl_map, self.chebi_map, self.drugbank_map).chembl_id, self.chembl_map.get(self.ligands[n].HET))
            self.assertEqual(LigandMap.worker(self.ligands[n], self.chembl_map, self.chebi_map, self.drugbank_map).chebi_id, self.chebi_map.get(self.ligands[n].HET))
            self.assertEqual(LigandMap.worker(self.ligands[n], self.chembl_map, self.chebi_map, self.drugbank_map).drugbank_id, self.drugbank_map.get(self.ligands[n].HET))