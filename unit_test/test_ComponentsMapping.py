import unittest
import io
from unittest import mock
import resources.ComponentsMapping
from models import Ligand

class TestComponentsMapping(unittest.TestCase):
    """
    UnitTest for Ligand annotations
    """

    @mock.patch("builtins.open")
    def test_parseCCD(self, open_mock):
        open_mock.return_value = io.StringIO("""\
        data_000
        #
        loop_ 
        _pdbe_chem_comp_external_mappings.resource"
        _pdbe_chem_comp_external_mappings.resource_id"
        D12 UniChem        ChEMBL      CHEMBL30959"
        D12 UniChem         ChEBI            28817"
        #
        loop_
        _pdbe_chem_comp_external_mappings.resource
        _pdbe_chem_comp_external_mappings.resource_id
        CA UniChem      DrugBank    DB14577
        CA UniChem         ChEBI      29108
        #
        loop_
        _pdbe_chem_comp_external_mappings.resource
        _pdbe_chem_comp_external_mappings.resource_id
        D10 UniChem     ChEMBL     CHEMBL134537
        D10 UniChem      ChEBI            41808""")
        chembl_map = {'D12': 'CHEMBL30959', 'D10': 'CHEMBL134537'}
        chebi_map = {'D12': '28817', 'CA': '29108'}
        drugbank_map = {'CA': 'DB14577'}
        self.assertEqual(resources.ComponentsMapping.parseCCD("ccd_file"), (chembl_map, chebi_map, drugbank_map))

    def ligand(self, emdb_id, sample_id, HET, lig_name, lig_copies, chembl_id, chebi_id, drugbank_id, provenance_chembl, provenance_chebi, provenance_drugbank):
        lig = Ligand(emdb_id, sample_id)
        lig.HET = HET
        lig.lig_name = lig_name
        lig.lig_copies = lig_copies
        lig.chembl_id = chembl_id
        lig.chebi_id = chebi_id
        lig.drugbank_id = drugbank_id
        lig.provenance_chembl = provenance_chembl
        lig.provenance_chebi = provenance_chebi
        lig.provenance_drugbank = provenance_drugbank
        return lig

    def test_worker(self):
        self.ligands = [
            self.ligand("EMD-8959", "2", "CA", "CALCIUM ION", "4", None, "29108", "DB14577", None, "PDBe-CCD", "PDBe-CCD"),
            self.ligand("EMD-8959", "3", "D12", "DODECANE", "8", "CHEMBL30959", "28817", None, "PDBe-CCD", "PDBe-CCD", None),
            self.ligand("EMD-8959", "4", "D10", "DECANE", "2", "CHEMBL134537", None, None, "PDBe-CCD", None, "PDBe-CCD")]
        chembl_map = {'D12': 'CHEMBL30959', 'D10': 'CHEMBL134537'}
        chebi_map = {'D12': '28817', 'CA': '29108'}
        drugbank_map = {'CA': 'DB14577'}
        LigandMap = resources.ComponentsMapping.ComponentsMapping(self.ligands)
        for n in range(len(self.ligands)):
            self.assertEqual(LigandMap.worker(self.ligands[n], chembl_map, chebi_map, drugbank_map).chembl_id, chembl_map.get(self.ligands[n].HET))
            self.assertEqual(LigandMap.worker(self.ligands[n], chembl_map, chebi_map, drugbank_map).chebi_id, chebi_map.get(self.ligands[n].HET))
            self.assertEqual(LigandMap.worker(self.ligands[n], chembl_map, chebi_map, drugbank_map).drugbank_id, drugbank_map.get(self.ligands[n].HET))