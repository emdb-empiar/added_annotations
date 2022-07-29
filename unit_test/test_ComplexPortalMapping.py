import unittest
import resources.ComplexPortalMapping
from models import CPX
import mock
from unit_test.setters import set_protein, set_supra, set_EMDB_comp


class TestComplexPortalMapping(unittest.TestCase):
    """
    UnitTest for Complex Portal annotations
    """

    def setUp(self):
        super(TestComplexPortalMapping, self).setUp()
        self.uniprot_side_effect = [['CPX-4881', 'CPX-4883', 'CPX-4885'], ['CPX-4881', 'CPX-4883', 'CPX-4885'], [],
                        ['CPX-4881', 'CPX-4883', 'CPX-4885'], ['CPX-4881', 'CPX-4883', 'CPX-4885'], [],
                        []]
        self.row = [['CPX-4881', 'DNA-directed RNA polymerase holoenzyme complex, Sigma70 variant', 'DNA', '83333', 'CHEBI:18420(1)|CHEBI:29105(2)|P00579(1)|P0A7Z4(2)|P0A800(1)|P0A8T7(1)|P0A8V2(1)', 'ECO:0000353(evidence used in manual assertion)', 'intact:EBI-22113905', 'GO:0003677(DNA binding)|GO:0003899(DNA-directed RNA polymerase activity)', 'pubmed:28666008(see-also)|pubmed:30109846(see-also)|wwpdb:3lu0(subset)|complex portal:CPX-4883(complex-primary)|complex portal:CPX-4881(inferred-from)|pubmed:21398637(see-also)|rhea:RHEA:21248(identity)|intenz:2.7.7.6(identity)'],
               ['CPX-4883', 'DNA-directed RNA polymerase holoenzyme complex, SigmaS variant', 'DNA', '83333', 'CHEBI:18420(1)|CHEBI:29105(2)|P0A7Z4(2)|P0A800(1)|P0A8T7(1)|P0A8V2(1)|P13445(1)', 'ECO:0005546(biological system)', 'intact:EBI-22113905', 'GO:2000142(regulation)', 'pubmed:28666008(see-also)|pubmed:30109846(see-also)|wwpdb:3lu0(subset)|complex portal:CPX-4883(complex-primary)|complex portal:CPX-4881(inferred-from)|pubmed:21398637(see-also)|rhea:RHEA:21248(identity)|intenz:2.7.7.6(identity)'],
               ['CPX-4885', 'DNA-directed RNA polymerase holoenzyme complex, SigmaE variant', 'DNA', '83333', 'CHEBI:18420(1)|CHEBI:29105(2)|P0A7Z4(2)|P0A800(1)|P0A8T7(1)|P0A8V2(1)|P0AGB6(1)', 'ECO:0000353(physical interaction evidence)', 'intact:EBI-22113905', 'GO:0090605(submerged biofilm formation)|GO:0036460(cellular response to cell envelope stress)', 'pubmed:28666008(see-also)|pubmed:30109846(see-also)|wwpdb:3lu0(subset)|complex portal:CPX-4883(complex-primary)|complex portal:CPX-4881(inferred-from)|pubmed:21398637(see-also)|rhea:RHEA:21248(identity)|intenz:2.7.7.6(identity)']]

        self.cpx_side_effect = [CPX(self.row[0]), CPX(self.row[1]), CPX(self.row[2]),
                                CPX(self.row[0]), CPX(self.row[1]), CPX(self.row[2]),
                                CPX(self.row[0]), CPX(self.row[1]), CPX(self.row[2])]
        self.proteins = [
            set_protein("EMD-0001", "1", "RNA polymerase-sigma54 holoenzyme", "83333", ['6GFW'], ['1', '2'],
                         "P0A7Z4", "UNIPROT",
                         "MQGSVTEFLKPRLVDIEQVSSTHAKVTLEPLERGFGHTLGNALRRILLSSMPGCAVTEVEIDGVLHEYSTKEGVQEDILEILLNLKG"
                         "LAVRVQGKDEVILTLNKSGIGPVTAADITHDGDVEIVKPQHVICHLTDENASISMRIKVQRGRGYVPASTRIHSEEDERPIGRLLVDA"
                         "CYSPVERIAYNVEAARVEQRTDLDKLVIEMETNGTIDPEEAIRRAATILAEQLEAFVDLRDVRQPEVKEEKPEFDPILLRPVDDLELT"
                         "VRSANCLKAEAIHYIGDLVQRTEVELLKTPNLGKKSLTEIKDVLASRGLSLGMRLENWPPASIADE",
                         "2", "", "", "", "", "", "", "", "", ""),
            set_protein("EMD-0001", "2", "DNA-directed RNA polymerase", "83333", ['6GFW'], ['1', '2'],
                         "P0A8T7", "UNIPROT",
                         "MVYSYTEKKRIRKDFGKRPQVLDVPYLLSIQLDSFQKFIEQDPEGQYGLEAAFRSVFPIQSYSGNSELQYVSYRLGEPVF DVQECQIRGVTYS"
                         "APLRVKLRLVIYEREAPEGTVKDIKEQEVYMGEIPLMTDNGTFVINGTERVIVSQLHRSPGVFFDSD KGKTHSSGKVLYNARIIPYRGSWLDFE"
                         "FDPKDNLFVRIDRRRKLPATIILRALNYTTEQILDLFFEKVIFEIRDNKLQME LVPERLRGETASFDIEANGKVYVEKGRRITARHIRQLEKDDV"
                         "KLIEVPVEYIAGKVVAKDYIDESTGELICAANMELSLD LLAKLSQSGHKRIETLFTNDLDHGPYISETLRVDPTNDRLSALVEIYRMMRPGEPPT",
                         "1", "", "", "", "", "", "", "", "", ""),
            set_protein("EMD-0001", "4", "DNA", "83333", "", "" , "P0A7Z4", "UNIPROT",
                         "MARVTVQDAVEKIGNRFDLVLVAARRARQMQVGGKDPLVPEENDKTTVIALREIEEGLINNQILDVRERQEQQEQEAAEL QAVTAIAEGRR",
                         "1", "", "", "", "", "", "", "", "", "")]
        self.supras = [set_supra('EMD-0001', 'supra_1', 'RNA polymerase-sigma54 holoenzyme', 'supra'),
                  set_supra('EMD-0001', 'supra_2', 'DNA-directed RNA polymerase', 'supra'),
                  set_supra('EMD-0001', 'supra_4', 'DNA', 'supra')]
        self.EMDB_complex = [set_EMDB_comp('EMD-0001', 'EMD-0001_1', 'RNA polymerase-sigma54 holoenzyme', '1', '', ['CPX-4881', 'CPX-4883', 'CPX-4885'],
                                       {'P0A800', 'P0A8T7', 'P0A7Z4'}, 'Complex Portal', ''),
                        set_EMDB_comp('EMD-0001', 'EMD-0001_2', 'DNA-directed RNA polymerase', '1', '', ['CPX-4881', 'CPX-4883', 'CPX-4885'],
                                       {'P0A800', 'P0A8T7', 'P0A7Z4'}, 'Complex Portal', ''),
                        set_EMDB_comp('EMD-0001', 'EMD-0001_4', 'DNA', '1', '',[], {'P0987'}, '', '')]



    # @patch.object(resources.ComplexPortalMapping.CPX_database, 'cpx_database', return_value='cpx_db')
    @mock.patch("resources.ComplexPortalMapping.CPX_database.get_from_cpx")
    @mock.patch("resources.ComplexPortalMapping.CPX_database.get_from_uniprot")
    def test_worker(self, uniprot_mock, cpx_mock):
        uniprot_mock.side_effect = self.uniprot_side_effect
        cpx_mock.side_effect = self.cpx_side_effect
        cpx_list = [('CPX-4881', 'CPX-4883', 'CPX-4885'), ('CPX-4881', 'CPX-4883', 'CPX-4885'), '']
        score = (0.6, 0.6, 0.0)
        provenance = ("Complex Portal", "Complex Portal", "")
        ComplexMap = resources.ComplexPortalMapping.CPMapping(self.proteins, self.supras, "CP_FTP")

        self.assertEqual(len(self.proteins), 3)
        for n in range(len(self.proteins)-1):
            EM_comp = ComplexMap.worker(self.EMDB_complex[n])
            cpx_id_list = ()
            for x in range(len(EM_comp.cpx_list)):
                cpx_id_list = cpx_id_list + (EM_comp.cpx_list[x].cpx_id,)
            self.assertEqual(cpx_id_list, cpx_list[n])
            self.assertEqual(EM_comp.score, score[n])
            self.assertEqual(EM_comp.provenance, provenance[n])

        self.assertIsNone(ComplexMap.worker(self.EMDB_complex[-1]))