import unittest
import io
import resources.PublicationMapping
from models import Citation, Author
import mock

class TestPublicationMapping(unittest.TestCase):
    """
    UnitTest for Citation annotations
    """
    def citation(self, emdb_id, pmedid, pmcid, doi, issn, authors, status, title, provenance_pm, provenance_pmc, provenance_issn, provenance_doi, provenance_orcid, url):
        cite = Citation(emdb_id)
        cite.emdb_id = emdb_id
        cite.pmedid = pmedid
        cite.pmcid = pmcid
        cite.doi = doi
        cite.issn = issn
        cite.authors = authors
        cite.status = status
        cite.title = title
        cite.provenance_pm = provenance_pm
        cite.provenance_pmc = provenance_pmc
        cite.provenance_issn = provenance_issn
        cite.provenance_doi = provenance_doi
        cite.provenance_orcid = provenance_orcid
        cite.url = url
        return cite

    @mock.patch("builtins.open")
    def test_generate_pubmed_dictionary(self, open_mock):
        open_mock.return_value = io.StringIO(
        "PMID\tPMC\tDOI\tISSN\t[AUTHORS ORCID]\n"
        "34812732\tPMC8735999\t10.7554/elife.73724\t2050-084X\t0000-0002-5119-3039	0000-0002-6290-8853\n"
        "21224846\tPMC3105306\t10.1038/ncomms1154\t2041-1723\t0000-0003-2681-5921	0000-0003-4835-154X	0000-0002-8070-1493\n")
        self.pubmed_dict = {'34812732': {'pmid': '34812732', 'pmcid': 'PMC8735999', 'doi': '10.7554/elife.73724', 'issn': '2050-084X', 'authors': ['0000-0002-5119-3039', '0000-0002-6290-8853']},
                       '21224846': {'pmid': '21224846', 'pmcid': 'PMC3105306', 'doi': '10.1038/ncomms1154', 'issn': '2041-1723', 'authors': ['0000-0003-2681-5921', '0000-0003-4835-154X', '0000-0002-8070-1493']}
        }
        self.assertEqual.__self__.maxDiff = None
        self.assertEqual(resources.PublicationMapping.generate_pubmed_dictionary("filename"), self.pubmed_dict)

    def test_worker(self):
        pubmed_dict = {
            '34812732': {'pmid': '34812732', 'pmcid': 'PMC8735999', 'doi': '10.7554/elife.73724', 'issn': '2050-084X',
                         'authors': ['0000-0002-5119-3039', '0000-0002-6290-8853']},
            '21224846': {'pmid': '21224846', 'pmcid': 'PMC3105306', 'doi': '10.1038/ncomms1154', 'issn': '2041-1723',
                         'authors': ['0000-0003-2681-5921', '0000-0003-4835-154X', '0000-0002-8070-1493']}
            }
        pmedid = '34812732'
        pmcid = 'PMC8735999'
        doi = '10.7554/elife.73724'
        issn = '2050-084X'
        authors = [('0000-0002-5119-3039', '1', 'EuropePMC'), ('0000-0002-6290-8853', '2', 'EuropePMC')]
        provenance = "EuropePMC"
        cite = self.citation("EMD-13501", "34812732", "PMC8735999", "10.7554/elife.73724", "2050-084X", [Author('', '1', '0000-0002-5119-3039', 'EuropePMC'), Author('', '2', '0000-0002-6290-8853', 'EuropePMC')],
                          "", "High-resolution structures of the actomyosin-V complex in three nucleotide states provide insights into the force generation mechanism.",
                        "EMDB", "EuropePMC", "EuropePMC", "EuropePMC", "EuropePMC", "")
        CitationMap = resources.PublicationMapping.PublicationMapping(cite)
        pub = CitationMap.worker(pubmed_dict)
        self.assertEqual(pub.pmedid, pmedid)
        self.assertEqual(pub.pmcid, pmcid)
        self.assertEqual(pub.provenance_pmc, provenance)
        self.assertEqual(pub.doi, doi)
        self.assertEqual(pub.provenance_doi, provenance)
        self.assertEqual(pub.issn, issn)
        self.assertEqual(pub.provenance_issn, provenance)
        for n in range(len(pub.authors)):
            self.assertEqual(pub.authors[n].orcid, authors[n][0])
            self.assertEqual(pub.authors[n].order, authors[n][1])
            self.assertEqual(pub.authors[n].provenance, authors[n][2])