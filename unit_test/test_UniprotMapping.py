import unittest
import io, os
import resources.UniprotMapping
from models import Protein
import mock

class TestUniprotMapping(unittest.TestCase):
    """
    UnitTest for Uniprot annotations
    """
    def __init__(self, workDir, blast_db, blastp_bin):
        super().__init__()
        self.workDir = workDir
        self.blast_db = blast_db
        self.blastp_bin = blastp_bin
        self.uniprot_dict = {'6GFW': [('P0A7Z4', 'DNA-directed RNA polymerase subunit alpha'), ('P0A8T7', 'DNA-directed RNA polymerase subunit beta')]}

    def protein(self, emdb_id, sample_id, sample_name, sample_organism, pdb, sample_complexes, uniprot_id, provenance,
                sequence, sample_copies, go, interpro, pfam, cath, scop, scop2, scop2B, pdbekb, alphafold):
        prot = Protein(emdb_id, sample_id)
        prot.emdb_id = emdb_id
        prot.sample_id = sample_id
        prot.sample_name = sample_name
        prot.sample_organism = sample_organism
        prot.pdb = pdb
        prot.sample_complexes = sample_complexes
        prot.uniprot_id = uniprot_id
        prot.provenance = provenance
        prot.sequence = sequence
        prot.sample_copies = sample_copies
        prot.go = go
        prot.interpro = interpro
        prot.pfam = pfam
        prot.cath = cath
        prot.scop = scop
        prot.scop2 = scop2
        prot.scop2B = scop2B
        prot.pdbekb = pdbekb
        prot.alphafold = alphafold
        return prot

    @mock.patch("builtins.open")
    def test_generate_unp_dictionary(self, open_mock):
        open_mock.return_value = io.StringIO("""\
        uni_id	pdb_id	name
        A0A0J9X294	5AFT;	Dynactin subunit 2
        A0A5S8WF48	6EE9;	Stress-response Peptide-1
        R4GRT5	3ZKT;	Tau-cnva""")
        uniprot = {'5aft': [('A0A0J9X294', 'Dynactin subunit 2')], '6ee9': [('A0A5S8WF48', 'Stress-response Peptide-1')], '3zkt': [('R4GRT5','Tau-cnva')]}
        uniprot_with_models = {'A0A0J9X294', 'A0A5S8WF48', 'R4GRT5'}
        self.assertEqual(resources.UniprotMapping.generate_unp_dictionary("filename"), (uniprot, uniprot_with_models))

    def test_worker(self):
        self.proteins = [
            self.protein("EMD-0001", "1", "DNA-directed RNA polymerase subunit alpha", "83333", ['6GFW'], ['1', '2'],
                         "P0A7Z4", "UNIPROT",
                         "MQGSVTEFLKPRLVDIEQVSSTHAKVTLEPLERGFGHTLGNALRRILLSSMPGCAVTEVEIDGVLHEYSTKEGVQEDILEILLNLKG"
                         "LAVRVQGKDEVILTLNKSGIGPVTAADITHDGDVEIVKPQHVICHLTDENASISMRIKVQRGRGYVPASTRIHSEEDERPIGRLLVDA"
                         "CYSPVERIAYNVEAARVEQRTDLDKLVIEMETNGTIDPEEAIRRAATILAEQLEAFVDLRDVRQPEVKEEKPEFDPILLRPVDDLELT"
                         "VRSANCLKAEAIHYIGDLVQRTEVELLKTPNLGKKSLTEIKDVLASRGLSLGMRLENWPPASIADE",
                         "2", "", "", "", "", "", "", "", "", ""),
            self.protein("EMD-0001", "2", "DNA-directed RNA polymerase subunit beta", "83333", ['6GFW'], ['1', '2'],
                         "P0A8T7", "UNIPROT",
                         "MVYSYTEKKRIRKDFGKRPQVLDVPYLLSIQLDSFQKFIEQDPEGQYGLEAAFRSVFPIQSYSGNSELQYVSYRLGEPVF DVQECQIRGVTYS"
                         "APLRVKLRLVIYEREAPEGTVKDIKEQEVYMGEIPLMTDNGTFVINGTERVIVSQLHRSPGVFFDSD KGKTHSSGKVLYNARIIPYRGSWLDFE"
                         "FDPKDNLFVRIDRRRKLPATIILRALNYTTEQILDLFFEKVIFEIRDNKLQME LVPERLRGETASFDIEANGKVYVEKGRRITARHIRQLEKDDV"
                         "KLIEVPVEYIAGKVVAKDYIDESTGELICAANMELSLD LLAKLSQSGHKRIETLFTNDLDHGPYISETLRVDPTNDRLSALVEIYRMMRPGEPPT",
                         "1", "", "", "", "", "", "", "", "", ""),
            self.protein("EMD-0001", "4", "DNA-directed RNA polymerase subunit omega", "83333", "", "", "P0A800", "UNIPROT",
                         "MARVTVQDAVEKIGNRFDLVLVAARRARQMQVGGKDPLVPEENDKTTVIALREIEEGLINNQILDVRERQEQQEQEAAEL QAVTAIAEGRR",
                         "1", "", "", "", "", "", "", "", "", "")]

        uniprot_model = {'6GFW': [('P0A7Z4', 'DNA-directed RNA polymerase subunit alpha'), ('P0A8T7', 'DNA-directed RNA polymerase subunit beta')]}
        uniprot_seq = ['P0A800', 'DNA-directed RNA polymerase subunit omega']
        self.ProteinMap = resources.UniprotMapping.UniprotMapping(self.workDir, self.proteins, self.uniprot_dict, self.blast_db, self.blastp_bin)
        for n in range(len(self.proteins)):
            if self.proteins[n].pdb:
                for x in range(len(self.proteins[n].pdb)):
                    uni_list = uniprot_model.get(self.proteins[n].pdb[x])
                    self.assertEqual(self.ProteinMap.worker(self.proteins[n]).uniprot_id, uni_list[n][x])
            if not self.proteins[n].pdb:
                self.assertEqual(self.ProteinMap.worker(self.proteins[n]).uniprot_id, uniprot_seq[0])

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

    def test_blastp(self):
        uni_map = ["P0A7Z4", "UniProt"]
        self.ProteinMap = resources.UniprotMapping.UniprotMapping(self.workDir, self.proteins, self.uniprot_dict, self.blast_db, self.blastp_bin)
        self.assertEqual(self.ProteinMap.blastp(self.proteins[0]).uniprot_id, uni_map[0])
        self.assertEqual(self.ProteinMap.blastp(self.proteins[0]).provenance, uni_map[1])
