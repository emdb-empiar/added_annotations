import unittest
import resources.ComponentsMapping
from models import Ligand

class TestComponentsMapping(unittest.TestCase):
    """
    UnitTest for Ligand annotations
    """
    def setUp(self):
        super(TestComponentsMapping, self).setUp()
        self.ligands = [
            Ligand("EMD-8959", "2", HET="CA", lig_name="CALCIUM ION", lig_copies="4"),
            Ligand("EMD-8959", "3", HET="D12", lig_name="DODECANE", lig_copies="8"),
            Ligand("EMD-8959", "4", HET="D10", lig_name="DECANE", lig_copies="2")]
        self.chembl_map = {'D12': 'CHEMBL30959', 'D10': 'CHEMBL134537'}
        self.chebi_map = {'D12': '28817', 'CA': '29108'}
        self.drugbank_map = {'CA': 'DB14577'}


    def test_worker(self):
        LigandMap = resources.ComponentsMapping.ComponentsMapping(self.ligands)
        for n in range(len(self.ligands)):
            self.assertEqual(LigandMap.worker(self.ligands[n], self.chembl_map, self.chebi_map, self.drugbank_map).chembl_id, self.chembl_map.get(self.ligands[n].HET))
            self.assertEqual(LigandMap.worker(self.ligands[n], self.chembl_map, self.chebi_map, self.drugbank_map).chebi_id, self.chebi_map.get(self.ligands[n].HET))
            self.assertEqual(LigandMap.worker(self.ligands[n], self.chembl_map, self.chebi_map, self.drugbank_map).drugbank_id, self.drugbank_map.get(self.ligands[n].HET))