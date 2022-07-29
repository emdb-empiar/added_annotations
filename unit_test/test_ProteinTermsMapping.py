import unittest
import resources.ProteinTermsMapping
from models import Model
from unit_test.setters import set_protein, set_go, set_ipro, set_csss, set_kbalpha


class TestProteinTermsMapping(unittest.TestCase):
    """
    UnitTest for GO, InterPro, Pfam, CATH, SCOP, SCOP2, SCOP2B, pdbekb and alphafold DB
    """
    def setUp(self):
        super(TestProteinTermsMapping, self).setUp()
        self.sifts_prefix = "/"
        self.is_go = True
        self.is_interpro = True
        self.is_pfam = True
        self.is_cath = True
        self.is_scop = True
        self.is_scop2 = True
        self.is_scop2B = True
        self.is_pdbekb = True
        self.is_AFDB = True


    def test_execute(self):
        model = Model("EMD-0001", '6GFW')
        proteins = [
            set_protein("EMD-0001", "1", "DNA-directed RNA polymerase subunit alpha", "83333", [model], ['1', '2'],
                         "P0A7Z4", "UNIPROT",
                         "MQGSVTEFLKPRLVDIEQVSSTHAKVTLEPLERGFGHTLGNALRRILLSSMPGCAVTEVEIDGVLHEYSTKEGVQEDILEILLNLKG"
                         "LAVRVQGKDEVILTLNKSGIGPVTAADITHDGDVEIVKPQHVICHLTDENASISMRIKVQRGRGYVPASTRIHSEEDERPIGRLLVDA"
                         "CYSPVERIAYNVEAARVEQRTDLDKLVIEMETNGTIDPEEAIRRAATILAEQLEAFVDLRDVRQPEVKEEKPEFDPILLRPVDDLELT"
                         "VRSANCLKAEAIHYIGDLVQRTEVELLKTPNLGKKSLTEIKDVLASRGLSLGMRLENWPPASIADE",
                         "2", [set_go('GO:0003677', 'DNA binding', 'molecular function', 'P0A7Z4', 'PDBe'),
                               set_go('GO:0005737', 'cytopasm', 'cellular component', 'P0A7Z4', 'PDBe')],
                         [set_ipro('IPR011263', 'DNA-dir_RNA_pol_RpoA/D/Rpb3', 'P0A7Z4', 'PDBe', '23', '234', '23', '234'),
                          set_ipro('IPR011260', 'RNAP_asu_C', 'P0A7Z4', 'PDBe', '245', '308', '245', '308')],
                         [set_ipro('PF01193','RNA_pol_L', 'P0A7Z4', 'PDBe', '29', '228', '29', '228'),
                          set_ipro('PF03118', 'RNA_pol_A_CTD', 'P0A7Z4', 'PDBe', '243', '309', '243', '309')],
                         [set_csss('3.30.1360.10', 'P0A7Z4', 'PDBe', '4', '238', '4', '238'),
                          set_csss('2.170.120.12', 'P0A7Z4', 'PDBe', '53', '177', '53', '177')],
                         [set_csss('SF-DOMID:80359', 'P0A7Z4', 'PDBe', '1', '230', '1', '230')],
                         [set_csss('SF-DOMID:8035907', 'P0A7Z4', 'PDBe', '4', '232', '4', '232')],
                         [set_csss('SF-DOMID:8042964', 'P0A7Z4', 'PDBe', '250', '321', '250', '321')],
                         set_kbalpha('P0A7Z4', 'UniProt'), set_kbalpha('P0A7Z4', 'AlphaFold DB')),
            set_protein("EMD-0001", "4", "DNA-directed RNA polymerase subunit omega", "83333", [model], "", "P0A800",
                         "UNIPROT",
                         "MARVTVQDAVEKIGNRFDLVLVAARRARQMQVGGKDPLVPEENDKTTVIALREIEEGLINNQILDVRERQEQQEQEAAEL QAVTAIAEGRR",
                         "1", [set_go('GO:0030880', 'RNA polymerase complex', 'cellular component', 'P0A800', 'PDBe'),
                               set_go('GO:0003899', 'DNA-directed 5-3 RNA polymerase activity', 'molecular function', 'P0A800', 'PDBe')],
                         [set_ipro('IPR036161', 'RPB6/omega-like_sf', 'P0A800', 'PDBe', '1', '79', '1', '79')],
                         [set_ipro('PF01192', 'RNA_pol_Rpb6', 'P0A800', 'PDBe', '8', '60', '8', '60')],
                         [set_csss('3.90.940.10', 'P0A800', 'PDBe', '1', '75', '1', '75')], "", "", "",
                         set_kbalpha('P0A800', 'UniProt'), set_kbalpha('P0A800', 'AlphaFold DB'))]
        uniprot_model = {'6GFW': [('P0A7Z4', 'DNA-directed RNA polymerase subunit alpha'),
                                  ('P0A800', 'DNA-directed RNA polymerase subunit omega')]}
        alphafold_ids = {'Q54TQ2', 'O07175', 'K7KVD5', 'Q9W1U1', 'P0A7Z4', 'P0A800'}

        ProteinMap = resources.ProteinTermsMapping.ProteinTermsMapping(proteins, self.sifts_prefix, alphafold_ids,
                                                                            self.is_go, self.is_interpro, self.is_pfam, self.is_cath, self.is_scop,
                                                                            self.is_scop2, self.is_scop2B, self.is_pdbekb, self.is_AFDB)

        self.assertEqual(len(proteins), 2)
        for n in range(len(proteins)):
            PT = ProteinMap.execute(uniprot_model)
            lg = len(proteins[n].go)
            self.assertGreater(lg, 0)
            for x in range(lg):
                self.assertEqual(PT[n].go[x], proteins[n].go[x])
            li = len(proteins[n].interpro)
            self.assertGreater(li, 0)
            for x in range(li):
                self.assertEqual(PT[n].interpro[x], proteins[n].interpro[x])
            lp = len(proteins[n].pfam)
            self.assertGreater(lp, 0)
            for x in range(lp):
                self.assertEqual(PT[n].pfam[x], proteins[n].pfam[x])
            lc = len(proteins[n].cath)
            self.assertGreater(lc, 0)
            for x in range(lc):
                self.assertEqual(PT[n].cath[x], proteins[n].cath[x])
            ls = len(proteins[n].scop)

            for x in range(ls):
                self.assertEqual(PT[n].scop[x], proteins[n].scop[x])
            ls2 = len(proteins[n].scop2)

            for x in range(ls2):
                self.assertEqual(PT[n].scop2[x], proteins[n].scop2[x])
            ls2b = len(proteins[n].scop2B)

            for x in range(ls2b):
                self.assertEqual(PT[n].scop2B[x], proteins[n].scop2B[x])
            self.assertEqual(PT[n].pdbekb, proteins[n].pdbekb)
            self.assertEqual(PT[n].alphafold, proteins[n].alphafold)

