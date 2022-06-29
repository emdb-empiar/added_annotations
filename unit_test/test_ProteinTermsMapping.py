import os
import unittest
from resources.ComplexPortalMapping import CPMapping
from resources.UniprotMapping import UniprotMapping

base_path = r'/Users/amudha/project/'

class TestUniprotMapping(unittest.TestCase):
    """
    UnitTest for UNIPROT annotation
    """
    def __init__(self):
        self.proteins = [
            Protein(emdb_id="EMD-0001", sample_id="1", sample_name ="DNA-directed RNA polymerase subunit alpha"),
            Protein(sample_organism="83333", pdb="[6GFW]"),
            Protein(sample_complexes="['1', '2']", uniprot_id ="P0A7Z4", provenance="UNIPROT"),
            Protein(sequence = "MQGSVTEFLKPRLVDIEQVSSTHAKVTLEPLERGFGHTLGNALRRILLSSMPGCAVTEVEIDGVLHEYSTKEGVQEDILEILLNLKG"
                               "LAVRVQGKDEVILTLNKSGIGPVTAADITHDGDVEIVKPQHVICHLTDENASISMRIKVQRGRGYVPASTRIHSEEDERPIGRLLVDA"
                               "CYSPVERIAYNVEAARVEQRTDLDKLVIEMETNGTIDPEEAIRRAATILAEQLEAFVDLRDVRQPEVKEEKPEFDPILLRPVDDLELT"
                               "VRSANCLKAEAIHYIGDLVQRTEVELLKTPNLGKKSLTEIKDVLASRGLSLGMRLENWPPASIADE"),
            Protein(sample_copies="2", go="set()", interpro="set()", pfam="set()", cath="set()", scop="set()", scop2="set()"),
            Protein(pdbekb="[]", alphafold="[]")
        ]

    def setUp(self):
        print("setUp")
        self.filename = os.path.join(base_path + "ONE_XML/EMD-3061/header/emd-3061-v30.xml")

    # def test_extracting_IDs(self):
    #     extracted_IDs = ['EMD-1831', 'EMDB_ID', '2xzb', 'PDB_ID', 'complex_id:CPX-2195', 'PDBe',
    #                      'P19156', 'UNIPROT', 'P05027', 'UNIPROT']
    #     self.assertListEqual(CPMapping.CPMapping.extracting_IDs(self, self.filename, self.pdbefile), extracted_IDs)

    def test_blastp(self):
        proteins = [
            Protein(emdb_id="EMD-0001", sample_id="1", sample_name="DNA-directed RNA polymerase subunit alpha"),
            Protein(sample_organism="83333", pdb="[6GFW]"),
            Protein(sample_complexes="['1', '2']", uniprot_id="P0A7Z4", provenance="UNIPROT"),
            Protein(sequence="MQGSVTEFLKPRLVDIEQVSSTHAKVTLEPLERGFGHTLGNALRRILLSSMPGCAVTEVEIDGVLHEYSTKEGVQEDILEILLNLKG"
                             "LAVRVQGKDEVILTLNKSGIGPVTAADITHDGDVEIVKPQHVICHLTDENASISMRIKVQRGRGYVPASTRIHSEEDERPIGRLLVDA"
                             "CYSPVERIAYNVEAARVEQRTDLDKLVIEMETNGTIDPEEAIRRAATILAEQLEAFVDLRDVRQPEVKEEKPEFDPILLRPVDDLELT"
                             "VRSANCLKAEAIHYIGDLVQRTEVELLKTPNLGKKSLTEIKDVLASRGLSLGMRLENWPPASIADE"),
            Protein(sample_copies="2", go="set()", interpro="set()", pfam="set()", cath="set()", scop="set()",
                    scop2="set()"),
            Protein(pdbekb="[]", alphafold="[]")
        ]
        self.fastafile = os.path.join(base_path + "fasta/EMD-10085.xml")
        self.ncbi_id = "559292"
        extracted_uniprot = ['Q14573', 'BLAST']
        self.assertListEqual(CPMapping.CPMapping.extract_uniprot_from_blast(self, self.fastafile,
                                                                        self.ncbi_id), extract_uniprot_from_blast)
    def test_extract_uniprot_from_blast(self):
        self.fastafile = os.path.join(base_path + "fasta/EMD-10085.xml")
        self.ncbi_id = "559292"
        extracted_uniprot = ['Q14573', 'BLAST']
        self.assertListEqual(CPMapping.CPMapping.extract_uniprot_from_blast(self, self.fastafile,
                                                                      self.ncbi_id), extract_uniprot_from_blast)

    def test_quering_IDs(self):
        self.mapFile = os.path.join(base_path + "EMDB_CPX_map.tsv")
        self.query = ['EMD-1831', 'EMDB_ID', '2xzb', 'PDB_ID', 'complex_id:CPX-2195', 'PDBe',
                         'P19156', 'UNIPROT', 'P05027', 'UNIPROT']
        complex_IDs = ['CPX-2195', 'CPX-57']
        complex_name = ['Hydrogen:potassium-exchanging ATPase complex', 'Sodium:potassium-exchanging ATPase complex']
        self.assertEqual(CPMapping.CPMapping.quering_IDs(self, self.query, self.mapFile), (complex_IDs, complex_name))

if __name__ == '__main__':
    unittest.main()
