import unittest
import io
import resources.UniprotMapping
from unit_test.setters import set_protein
import mock
from models import Model

class TestUniprotMapping(unittest.TestCase):
    """
    UnitTest for Uniprot annotations
    """

    def setUp(self):
        super(TestUniprotMapping, self).setUp()
        self.workDir = "/dummy"
        self.blast_db = "/blast_db"
        self.blastp_bin = "/blastp_bin"
        self.uniprot_dict = {'6GFW': [('P0A7Z4', 'DNA-directed RNA polymerase subunit alpha'), ('P0A8T7', 'DNA-directed RNA polymerase subunit beta')]}
        self.proteins = [
            set_protein("EMD-0001", "1", "DNA-directed RNA polymerase subunit alpha", "83333", [Model('EMD-0001', '6GFW')], ['1', '2'],
                         None, None, "", "2", "", "", "", "", "", "", "", "", ""),
            set_protein("EMD-0001", "2", "DNA-directed RNA polymerase subunit beta", "83333", [Model('EMD-0001', '6GFW')], ['1', '2'],
                         None, None, "", "1", "", "", "", "", "", "", "", "", ""),
            set_protein("EMD-0001", "4", "DNA-directed RNA polymerase subunit omega", "83333", "", "", None,
                         None, "", "1", "", "", "", "", "", "", "", "", "")]

    @mock.patch("builtins.open")
    def test_generate_unp_dictionary(self, file_mock):
        file_mock.return_value = io.StringIO("""\
        uni_id	pdb_id	name
        A0A0J9X294	5AFT;	Dynactin subunit 2
        A0A5S8WF48	6EE9;	Stress-response Peptide-1
        R4GRT5	3ZKT;	Tau-cnva""")
        uniprot = {'5aft': [('A0A0J9X294', 'Dynactin subunit 2')], '6ee9': [('A0A5S8WF48', 'Stress-response Peptide-1')], '3zkt': [('R4GRT5','Tau-cnva')]}
        uniprot_with_models = {'A0A0J9X294', 'A0A5S8WF48', 'R4GRT5'}
        self.assertEqual(resources.UniprotMapping.generate_unp_dictionary("filename"), (uniprot, uniprot_with_models))

    def test_worker(self):
        uniprot_model = {'6GFW': [('P0A7Z4', 'DNA-directed RNA polymerase subunit alpha'), ('P0A8T7', 'DNA-directed RNA polymerase subunit beta')]}
        uniprot_seq = ['P0A800', 'DNA-directed RNA polymerase subunit omega']
        ProteinMap = resources.UniprotMapping.UniprotMapping(self.workDir, self.proteins, self.uniprot_dict, self.blast_db, self.blastp_bin)

        self.assertEqual(ProteinMap.worker(self.proteins[0]).uniprot_id, "P0A7Z4")
        self.assertEqual(ProteinMap.worker(self.proteins[1]).uniprot_id, "P0A8T7")
        self.assertIsNone(ProteinMap.worker(self.proteins[2]).uniprot_id)

    @mock.patch("builtins.open")
    def test_extract_uniprot_from_blast(self, open_mock):
        open_mock.side_effect = [io.StringIO("""\
        <BlastOutput>
        <BlastOutput_query-len>651</BlastOutput_query-len>
        <Hit>
        <Hit_def>sp|Q12931|TRAP1_HUMAN Heat shock protein 75 kDa, mitochondrial OS=Homo sapiens OX=9606 GN=TRAP1 PE=1 SV=3</Hit_def>
        <Hit_accession>490430</Hit_accession>
        <Hit_len>704</Hit_len>
        </Hit>
        </BlastOutput>""")]
        ncbi_id = "9606"
        extracted_uniprot = "Q12931"
        self.assertEqual(resources.UniprotMapping.UniprotMapping.extract_uniprot_from_blast(self, "fastafile", ncbi_id), extracted_uniprot)