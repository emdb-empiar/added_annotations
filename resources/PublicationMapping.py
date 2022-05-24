import os
import itertools

def generate_orcid_dictionary(workDir):
    orcid_dict = {}
    epmc_orcid = os.path.join(workDir, "EPMC_orcid.log")
    with open(epmc_orcid, 'r') as f:
        for line in f.readlines()[1:]:
            line = line.strip('\n')
            emdb_id = line.split('\t')[0]
            name = line.split('\t')[1]
            id = line.split('\t')[2]
            order = line.split('\t')[3]
            pvn = line.split('\t')[4]
            if emdb_id not in orcid_dict:
                ind = -1
                orcid_dict[emdb_id] = {}
            if emdb_id in orcid_dict:
                ind = ind + 1
            orc_list = ["name_"+str(ind), name, "id_"+str(ind), id, "order_"+str(ind), order, "provenance_orcid_"+str(ind), pvn]
            list_dict = dict(itertools.zip_longest(*[iter(orc_list)] * 2, fillvalue=""))
            for k in list_dict.keys():
                orcid_dict[emdb_id][k] = list_dict[k]
            orcid_dict[emdb_id]["ind"] = ind + 1
    return orcid_dict

def generate_pubmed_dictionary(workDir):
    pubmed_dict = {}
    epmc_pubmed = os.path.join(workDir, "EPMC_pubmed.log")
    with open(epmc_pubmed, 'r') as f:
        for line in f.readlines()[1:]:
            line = line.strip('\n')
            emdb_id = line.split('\t')[0]
            pub_id = line.split('\t')[1]
            pmc_id = line.split('\t')[2]
            issn = line.split('\t')[3]
            doi = line.split('\t')[4]
            pvn = line.split('\t')[5]
            if emdb_id not in pubmed_dict:
                pubmed_dict[emdb_id] = {}
            pub_list = ["pubmed_id", pub_id, "pubmed_central", pmc_id, "issn", issn, "doi", doi, "provenance", pvn]
            list_dict = dict(itertools.zip_longest(*[iter(pub_list)] * 2, fillvalue=""))
            for k in list_dict.keys():
                pubmed_dict[emdb_id][k] = list_dict[k]
    return pubmed_dict

class PublicationMapping:
    """
     If pubmed id available then annotations collected from EuropePMC API (PubMed, PubMED central, DOI and ORCID IDs),
     if no pubmed id then author provided annotations for the publication.
    """

    def __init__(self, citations, is_pmc=True, pubmed_dict={}, is_orcid=True, orcid_dict={}):
        self.citations = citations
        self.is_pmc = is_pmc
        self.pubmed_dict = pubmed_dict
        self.is_orcid = is_orcid
        self.orcid_dict = orcid_dict

    def execute(self):
        for citation in self.citations:
            citation = self.worker(citation)
        return self.citations

    def worker(self, citation):
        if self.is_orcid:
            orc = self.orcid_dict.get(citation.emdb_id)
            if orc:
                citation.orcid_ids = orc
            if not orc:
                orcid_ids = {}
                ind = -1
                for ord_name, id in (citation.name_order).items():
                    ind = ind + 1
                    name = ord_name.split('_')[1]
                    order = ord_name.split('_')[0]
                    orc_list = ["name_" + str(ind), name, "id_" + str(ind), id, "order_" + str(ind), order,
                                "provenance_orcid_" + str(ind), citation.provenance_orcid]
                    list_dict = dict(itertools.zip_longest(*[iter(orc_list)] * 2, fillvalue=""))
                    for k in list_dict.keys():
                        orcid_ids[k] = list_dict[k]
                    orcid_ids["ind"] = ind + 1
                    citation.orcid_ids = orcid_ids

        if self.is_pmc:
            pm = self.pubmed_dict.get(citation.emdb_id)
            if pm:
                citation.pmedid = pm.get('pubmed_id')
                citation.pmcid = pm.get('pubmed_central')
                citation.issn = pm.get('issn')
                citation.doi = pm.get('doi')
                citation.provenance_pm = pm.get('provenance')
                citation.provenance_pmc = pm.get('provenance')
                citation.provenance_doi = pm.get('provenance')
                citation.provenance_issn = pm.get('provenance')
            if not pm:
                if citation.pmcid:
                    citation.provenance_pmc = "EMDB"
                if citation.doi:
                    citation.provenance_doi = "EMDB"
                if citation.pmedid:
                    citation.provenance_pm = "EMDB"
                if citation.issn:
                    citation.provenance_issn = "EMDB"
        return citation

    def export_tsv(self, pubmed_logger, orcid_logger):
        for citation in self.citations:
            row = f"{citation.emdb_id}\t{citation.pmedid}\t{citation.pmcid}\t{citation.issn}\t{citation.doi}"
            pubmed_logger.info(row)
            ind = (citation.orcid_ids).get("ind")
            if ind != None:
                for i in range(int(ind)):
                    name = (citation.orcid_ids).get('name_' + str(i))
                    id = (citation.orcid_ids).get('id_' + str(i))
                    order = (citation.orcid_ids).get('order_' + str(i))
                    pvn = (citation.orcid_ids).get('provenance_orcid_' + str(i))
                    line = f"{citation.emdb_id}\t{name}\t{id}\t{order}\t{pvn}"
                    orcid_logger.info(line)