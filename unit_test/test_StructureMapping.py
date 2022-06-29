import unittest
import io
from unittest import mock
import resources.StructureMapping
from models import Model

class TestStructureMapping(unittest.TestCase):
    """
    UnitTest for PDB ids and their overall molecular weight annotations
    """

    def model(self, emdb_id, pdb_id, assembly, molecular_weight):
        mod = Model(emdb_id, pdb_id)
        mod.emdb_id = emdb_id
        mod.pdb_id = pdb_id
        mod.assembly = assembly
        mod.molecular_weight = molecular_weight
        return mod

    @mock.patch("builtins.open")
    def test_worker(self, open_mock):
        self.models = [self.model("EMD-3061", "5a63", "3",	"38843.82"),
                       self.model("EMD-4678", "6gh5", "8", "483169.49")]
        open_mock.side_effect = [io.StringIO("""\
                <assembly composition="protein structure" id="1" molecular_weight="38843.82" name="monomer" order="3" prefered="True" type="homo">
                    <entity chain_ids="A" class="protein" count="1" entity_id="1" source_entry="2FVO" type="polymer"/>
                </assembly>
                """),
                                 io.StringIO("""\
                         <assembly composition="DNA/protein complex" id="1" molecular_weight="483169.4909999999" name="octamer" order="8" prefered="True" type="hetero">
                            <entity chain_ids="A,B" class="protein" count="2" entity_id="1" source_entry="6GH5" type="polymer"/>
                        </assembly>
                        """)
                                 ]
        test_mw = [38843.82, 483169.491]
        test_assembly = [3, 8]
        ModelMap = resources.StructureMapping.StructureMapping(self.models, "assembly_ftp")
        for n in range(len(self.models)):
            Model = ModelMap.worker(self.models[n])
            self.assertEqual(Model.molecular_weight, test_mw[n])
            self.assertEqual(Model.assembly, test_assembly[n])