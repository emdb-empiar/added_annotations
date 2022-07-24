import unittest
import io
import resources.ComplexPortalMapping
from models import Protein, CPX, EMDB_complex, Supra
import mock
from mock import patch

class TestComplexPortalMapping(unittest.TestCase):
    """
    UnitTest for Complex Portal annotations
    """
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

    def supra(self, emdb_id, supra_id, supra_name, kind):
        supra = Supra(emdb_id, supra_id)
        supra.emdb_id = emdb_id
        supra.supra_id = supra_id
        supra.supra_name = supra_name
        supra.kind = kind
        return supra

    def EMDB_comp(self, emdb_id, sample_id, supra_name, sample_copies, complex_sample_id, cpx_list, proteins, provenance, score):
        em_cpx = EMDB_complex(emdb_id, sample_id, supra_name, sample_copies, complex_sample_id)
        em_cpx.emdb_id = emdb_id
        em_cpx.sample_id = sample_id
        em_cpx.supra_name = supra_name
        em_cpx.sample_copies = sample_copies
        em_cpx.complex_sample_id = complex_sample_id
        em_cpx.cpx_list = cpx_list
        em_cpx.proteins = proteins
        em_cpx.provenance = provenance
        em_cpx.score = score
        return em_cpx

    # @patch.object(resources.ComplexPortalMapping.CPX_database, 'cpx_database', return_value='cpx_db')
    @mock.patch("resources.ComplexPortalMapping.CPX_database.get_from_cpx")
    @mock.patch("resources.ComplexPortalMapping.CPX_database.get_from_uniprot")
    def test_worker(self, uniprot_mock, cpx_mock):
        uniprot_mock.side_effect = [['CPX-4881', 'CPX-4883', 'CPX-4885'], ['CPX-4881', 'CPX-4883', 'CPX-4885'], [],
                        ['CPX-4881', 'CPX-4883', 'CPX-4885'], ['CPX-4881', 'CPX-4883', 'CPX-4885'], [],
                        []]
        row = [['CPX-4881', 'DNA-directed RNA polymerase holoenzyme complex, Sigma70 variant', 'DNA', '83333', 'CHEBI:18420(1)|CHEBI:29105(2)|P00579(1)|P0A7Z4(2)|P0A800(1)|P0A8T7(1)|P0A8V2(1)', 'ECO:0000353(evidence used in manual assertion)', 'intact:EBI-22113905', 'GO:0003677(DNA binding)|GO:0003899(DNA-directed RNA polymerase activity)', 'pubmed:28666008(see-also)|pubmed:30109846(see-also)|wwpdb:3lu0(subset)|complex portal:CPX-4883(complex-primary)|complex portal:CPX-4881(inferred-from)|pubmed:21398637(see-also)|rhea:RHEA:21248(identity)|intenz:2.7.7.6(identity)'],
               ['CPX-4883', 'DNA-directed RNA polymerase holoenzyme complex, SigmaS variant', 'DNA', '83333', 'CHEBI:18420(1)|CHEBI:29105(2)|P0A7Z4(2)|P0A800(1)|P0A8T7(1)|P0A8V2(1)|P13445(1)', 'ECO:0005546(biological system)', 'intact:EBI-22113905', 'GO:2000142(regulation)', 'pubmed:28666008(see-also)|pubmed:30109846(see-also)|wwpdb:3lu0(subset)|complex portal:CPX-4883(complex-primary)|complex portal:CPX-4881(inferred-from)|pubmed:21398637(see-also)|rhea:RHEA:21248(identity)|intenz:2.7.7.6(identity)'],
               ['CPX-4885', 'DNA-directed RNA polymerase holoenzyme complex, SigmaE variant', 'DNA', '83333', 'CHEBI:18420(1)|CHEBI:29105(2)|P0A7Z4(2)|P0A800(1)|P0A8T7(1)|P0A8V2(1)|P0AGB6(1)', 'ECO:0000353(physical interaction evidence)', 'intact:EBI-22113905', 'GO:0090605(submerged biofilm formation)|GO:0036460(cellular response to cell envelope stress)', 'pubmed:28666008(see-also)|pubmed:30109846(see-also)|wwpdb:3lu0(subset)|complex portal:CPX-4883(complex-primary)|complex portal:CPX-4881(inferred-from)|pubmed:21398637(see-also)|rhea:RHEA:21248(identity)|intenz:2.7.7.6(identity)']]
        cpx_mock.side_effect = [CPX(row[0]), CPX(row[1]), CPX(row[2]),
                                CPX(row[0]), CPX(row[1]), CPX(row[2]),
                                CPX(row[0]), CPX(row[1]), CPX(row[2])]
        proteins = [
            self.protein("EMD-0001", "1", "RNA polymerase-sigma54 holoenzyme", "83333", ['6GFW'], ['1', '2'],
                         "P0A7Z4", "UNIPROT",
                         "MQGSVTEFLKPRLVDIEQVSSTHAKVTLEPLERGFGHTLGNALRRILLSSMPGCAVTEVEIDGVLHEYSTKEGVQEDILEILLNLKG"
                         "LAVRVQGKDEVILTLNKSGIGPVTAADITHDGDVEIVKPQHVICHLTDENASISMRIKVQRGRGYVPASTRIHSEEDERPIGRLLVDA"
                         "CYSPVERIAYNVEAARVEQRTDLDKLVIEMETNGTIDPEEAIRRAATILAEQLEAFVDLRDVRQPEVKEEKPEFDPILLRPVDDLELT"
                         "VRSANCLKAEAIHYIGDLVQRTEVELLKTPNLGKKSLTEIKDVLASRGLSLGMRLENWPPASIADE",
                         "2", "", "", "", "", "", "", "", "", ""),
            self.protein("EMD-0001", "2", "DNA-directed RNA polymerase", "83333", ['6GFW'], ['1', '2'],
                         "P0A8T7", "UNIPROT",
                         "MVYSYTEKKRIRKDFGKRPQVLDVPYLLSIQLDSFQKFIEQDPEGQYGLEAAFRSVFPIQSYSGNSELQYVSYRLGEPVF DVQECQIRGVTYS"
                         "APLRVKLRLVIYEREAPEGTVKDIKEQEVYMGEIPLMTDNGTFVINGTERVIVSQLHRSPGVFFDSD KGKTHSSGKVLYNARIIPYRGSWLDFE"
                         "FDPKDNLFVRIDRRRKLPATIILRALNYTTEQILDLFFEKVIFEIRDNKLQME LVPERLRGETASFDIEANGKVYVEKGRRITARHIRQLEKDDV"
                         "KLIEVPVEYIAGKVVAKDYIDESTGELICAANMELSLD LLAKLSQSGHKRIETLFTNDLDHGPYISETLRVDPTNDRLSALVEIYRMMRPGEPPT",
                         "1", "", "", "", "", "", "", "", "", ""),
            self.protein("EMD-0001", "4", "DNA", "83333", "", "" , "P0A7Z4", "UNIPROT",
                         "MARVTVQDAVEKIGNRFDLVLVAARRARQMQVGGKDPLVPEENDKTTVIALREIEEGLINNQILDVRERQEQQEQEAAEL QAVTAIAEGRR",
                         "1", "", "", "", "", "", "", "", "", "")]
        supras = [self.supra('EMD-0001', 'supra_1', 'RNA polymerase-sigma54 holoenzyme', 'supra'),
                  self.supra('EMD-0001', 'supra_2', 'DNA-directed RNA polymerase', 'supra'),
                  self.supra('EMD-0001', 'supra_4', 'DNA', 'supra')]
        EMDB_complex = [self.EMDB_comp('EMD-0001', 'EMD-0001_1', 'RNA polymerase-sigma54 holoenzyme', '1', '', ['CPX-4881', 'CPX-4883', 'CPX-4885'],
                                       {'P0A800', 'P0A8T7', 'P0A7Z4'}, 'Complex Portal', '0.5'),
                        self.EMDB_comp('EMD-0001', 'EMD-0001_2', 'DNA-directed RNA polymerase', '1', '', ['CPX-4881', 'CPX-4883', 'CPX-4885'],
                                       {'P0A800', 'P0A8T7', 'P0A7Z4'}, 'Complex Portal', '0.6'),
                        self.EMDB_comp('EMD-0001', 'EMD-0001_4', 'DNA', '1', '',[], {'P0987'}, '', '')]

        cpx_list = [('CPX-4881', 'CPX-4883', 'CPX-4885'), ('CPX-4881', 'CPX-4883', 'CPX-4885'), '']
        score = (0.6, 0.6, 0.0)
        provenance = ("Complex Portal", "Complex Portal", "")
        ComplexMap = resources.ComplexPortalMapping.CPMapping(proteins, supras, "CP_FTP")
        for n in range(len(proteins)):
            EM_comp = ComplexMap.worker(EMDB_complex[n])
            cpx_id_list = ()
            if EM_comp != None:
                for x in range(len(EM_comp.cpx_list)):
                    cpx_id_list = cpx_id_list + (EM_comp.cpx_list[x].cpx_id,)
                self.assertEqual(cpx_id_list, cpx_list[n])
                self.assertEqual(EM_comp.score, score[n])
                self.assertEqual(EM_comp.provenance, provenance[n])
            else:
                if n == 2:
                    self.assertTrue(True)
                else:
                    self.assertTrue(False)