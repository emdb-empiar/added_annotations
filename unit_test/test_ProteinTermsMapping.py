import unittest
import resources.ProteinTermsMapping
from models import Protein, GO, Interpro, Cath, Pdbekb, Model

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

    def go(self, id, namespace, type, unip_id, provenance):
        go = GO()
        go.id = id
        go.namespace = namespace
        go.type = type
        go.unip_id = unip_id
        go.provenance = provenance
        return go
    
    def ipro(self, id, namespace, unip_id, provenance, start, end, unp_start, unp_end):
        intpfam = Interpro()
        intpfam.id = id
        intpfam.namespace = namespace
        intpfam.unip_id = unip_id
        intpfam.provenance = provenance
        intpfam.start = start
        intpfam.end = end
        intpfam.unp_start = unp_start
        intpfam.unp_end = unp_end
        return intpfam

    def csss(self, id, unip_id, provenance, start, end, unp_start, unp_end):
        csss = Cath()
        csss.id = id
        csss.unip_id = unip_id
        csss.provenance = provenance
        csss.start = start
        csss.end = end
        csss.unp_start = unp_start
        csss.unp_end = unp_end
        return csss

    def kbalpha(self, unip_id, provenance):
        kbalpha = Pdbekb(unip_id, provenance)
        kbalpha.unip_id = unip_id
        kbalpha.provenance = provenance
        return kbalpha

    def test_execute(self):
        model = Model("EMD-0001", '6GFW')
        proteins = [
            self.protein("EMD-0001", "1", "DNA-directed RNA polymerase subunit alpha", "83333", [model], ['1', '2'],
                         "P0A7Z4", "UNIPROT",
                         "MQGSVTEFLKPRLVDIEQVSSTHAKVTLEPLERGFGHTLGNALRRILLSSMPGCAVTEVEIDGVLHEYSTKEGVQEDILEILLNLKG"
                         "LAVRVQGKDEVILTLNKSGIGPVTAADITHDGDVEIVKPQHVICHLTDENASISMRIKVQRGRGYVPASTRIHSEEDERPIGRLLVDA"
                         "CYSPVERIAYNVEAARVEQRTDLDKLVIEMETNGTIDPEEAIRRAATILAEQLEAFVDLRDVRQPEVKEEKPEFDPILLRPVDDLELT"
                         "VRSANCLKAEAIHYIGDLVQRTEVELLKTPNLGKKSLTEIKDVLASRGLSLGMRLENWPPASIADE",
                         "2", [self.go('GO:0003677', 'DNA binding', 'molecular function', 'P0A7Z4', 'PDBe'),
                               self.go('GO:0005737', 'cytopasm', 'cellular component', 'P0A7Z4', 'PDBe')],
                         [self.ipro('IPR011263', 'DNA-dir_RNA_pol_RpoA/D/Rpb3', 'P0A7Z4', 'PDBe', '23', '234', '23', '234'),
                          self.ipro('IPR011260', 'RNAP_asu_C', 'P0A7Z4', 'PDBe', '245', '308', '245', '308')],
                         [self.ipro('PF01193','RNA_pol_L', 'P0A7Z4', 'PDBe', '29', '228', '29', '228'),
                          self.ipro('PF03118', 'RNA_pol_A_CTD', 'P0A7Z4', 'PDBe', '243', '309', '243', '309')],
                         [self.csss('3.30.1360.10', 'P0A7Z4', 'PDBe', '4', '238', '4', '238'),
                          self.csss('2.170.120.12', 'P0A7Z4', 'PDBe', '53', '177', '53', '177')],
                         [self.csss('SF-DOMID:80359', 'P0A7Z4', 'PDBe', '1', '230', '1', '230')],
                         [self.csss('SF-DOMID:8035907', 'P0A7Z4', 'PDBe', '4', '232', '4', '232')],
                         [self.csss('SF-DOMID:8042964', 'P0A7Z4', 'PDBe', '250', '321', '250', '321')],
                         self.kbalpha('P0A7Z4', 'UniProt'), self.kbalpha('P0A7Z4', 'AlphaFold DB')),
            self.protein("EMD-0001", "4", "DNA-directed RNA polymerase subunit omega", "83333", [model], "", "P0A800",
                         "UNIPROT",
                         "MARVTVQDAVEKIGNRFDLVLVAARRARQMQVGGKDPLVPEENDKTTVIALREIEEGLINNQILDVRERQEQQEQEAAEL QAVTAIAEGRR",
                         "1", [self.go('GO:0030880', 'RNA polymerase complex', 'cellular component', 'P0A800', 'PDBe'),
                               self.go('GO:0003899', 'DNA-directed 5-3 RNA polymerase activity', 'molecular function', 'P0A800', 'PDBe')],
                         [self.ipro('IPR036161', 'RPB6/omega-like_sf', 'P0A800', 'PDBe', '1', '79', '1', '79')],
                         [self.ipro('PF01192', 'RNA_pol_Rpb6', 'P0A800', 'PDBe', '8', '60', '8', '60')],
                         [self.csss('3.90.940.10', 'P0A800', 'PDBe', '1', '75', '1', '75')], "", "", "",
                         self.kbalpha('P0A800', 'UniProt'), self.kbalpha('P0A800', 'AlphaFold DB'))]
        uniprot_model = {'6GFW': [('P0A7Z4', 'DNA-directed RNA polymerase subunit alpha'),
                                  ('P0A800', 'DNA-directed RNA polymerase subunit omega')]}
        alphafold_ids = {'Q54TQ2', 'O07175', 'K7KVD5', 'Q9W1U1', 'P0A7Z4', 'P0A800'}

        ProteinMap = resources.ProteinTermsMapping.ProteinTermsMapping(proteins, self.sifts_prefix, alphafold_ids,
                                                                            self.is_go, self.is_interpro, self.is_pfam, self.is_cath, self.is_scop,
                                                                            self.is_scop2, self.is_scop2B, self.is_pdbekb, self.is_AFDB)
        for n in range(len(proteins)):
            PT = ProteinMap.execute(uniprot_model)
            lg = len(proteins[n].go)
            for x in range(lg):
                self.assertEqual(PT[n].go[x], proteins[n].go[x])
            li = len(proteins[n].interpro)
            for x in range(li):
                self.assertEqual(PT[n].interpro[x], proteins[n].interpro[x])
            lp = len(proteins[n].pfam)
            for x in range(lp):
                self.assertEqual(PT[n].pfam[x], proteins[n].pfam[x])
            lc = len(proteins[n].cath)
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

