import os
import itertools

def generate_pubmed_dictionary(workDir):
    pubmed_dict = {}
    epmc_pubmed = os.path.join(workDir, "EPMC_pubmed.tsv")
    with open(epmc_pubmed, 'r') as f:
        for line in f.readlines()[1:]:
            row = line.strip('\n').split('\t')
            pubmed_dict[row[0]] = {'pmid': row[0], 'pmcid': row[1], 'doi': row[2], 'issn': row[3], 'authors': row[4:]}
    return pubmed_dict

class PublicationMapping:
    """
     If pubmed id available then annotations collected from EuropePMC API (PubMed, PubMED central, DOI and ORCID IDs),
     if no pubmed id then author provided annotations for the publication.
    """

    def __init__(self, citations):
        self.citations = citations

    def execute(self, pubmed_dict):
        for citation in self.citations:
            citation = self.worker(citation, pubmed_dict)
        return self.citations

    def worker(self, citation, pubmed_dict):
        if citation.pmedid:
            pmid = citation.pmedid
            if pmid in pubmed_dict:
                pm = pubmed_dict[pmid]
                
                # ORCID
                authors = pm['authors']
                for i, orcid in enumerate(authors):
                    if orcid:
                        citation.addExternalOrcid(orcid, i+1, "EuropePMC")

                # PMC
                if pm['pmcid']:
                    citation.pmcid = pm['pmcid']
                    citation.provenance_pmc = "EuropePMC"
                if pm['doi']:
                    citation.doi = pm['doi']
                    citation.provenance_doi = "EuropePMC"
                if pm['issn']:
                    citation.issn = pm['issn']
                    citation.provenance_issn = "EuropePMC"
            
        return citation

    def export_tsv(self, pubmed_logger, orcid_logger):
        for citation in self.citations:
            if citation.pmedid or citation.doi:
                pubmed_logger.info(str(citation))
            for author in citation.authors:
                orcid_logger.info(f"{citation.emdb_id}\t{str(author)}")
