# Added annotations: EMICSS (**E**MDB **In**tegration with **C**omplexes, **S**tructures and **S**equences)

This repository provides tools and scripts for extracting and adding annotations to EMDB entries, which are used to enhance the metadata associated with EM datasets.

### Table of Contents

* Installation
* Configuration
* Usage
* Contributing
* License

### Installation

To install the necessary dependencies, run: 
pip install -r requirements.txt

### Configuration

The repository uses a config.ini file for configuration, which is not included in the repository. This file should be created in the root directory of the project with the following structure:

```
[file_paths]
uniprot_tab: <path_to_file>/uniprot.tsv
CP_ftp: <path_to_file>/complextab
components_cif: <path_to_file>/components.cif
chem_comp_list: <path_to_file>/chem_comp_list.xml
pmc_ftp_gz: <path_to_file>/PMID_PMCID_DOI.csv.gz
pmc_ftp: <path_to_file>/PMID_PMCID_DOI.csv
emdb_pubmed: <path_to_file>/emdb_pubmed.log
emdb_orcid: <path_to_file>/emdb_orcid.log
assembly_ftp: <path_to_file>/assembly/
BLAST_DB: <path_to_file>/ncbi-blast-2.13.0+/database/uniprot_sprot
BLASTP_BIN: blastp
sifts_GO: <path_to_file>/pdb_chain_go.csv
GO_obo: <path_to_file>/go.obo
GO_interpro: /nfs/ftp/pub/databases/GO/goa/external2go/interpro2go
sifts: <path_to_file>/split_xml/
alphafold_ftp: <path_to_file>/accession_ids.txt
rfam_ftp: <path_to_file>/rfam_files_combined.txt

[api]
pmc: https://www.ebi.ac.uk/europepmc/webservices/rest/searchPOST
```

#### File Sources and Download Links
| File        | 	Descritption         | 	Download Link                                                                                                                                     |	
|-------------|-----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------|
| uniprot.tsv | 	UniProt annpotations | 	 https://rest.uniprot.org/uniprotkb/stream?fields=accession,xref_pdb,protein_name&query=((database:pdb))&format=tsv&compressed=false              |
| complextab | 	Complex Portal data | 	 https://ftp.ebi.ac.uk/pub/databases/complexportal/complexes.tab.gz                                                                               |
| components.cif | 	Chemical components data | 	 https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/ccd/components.cif                                                                           |
| chem_comp_list.xml | 	Chemical component list | 	 https://ftp.ebi.ac.uk/pub/databases/msd/pdbechem_v2/ccd/chem_comp_list.xml                                                                       |
| PMID_PMCID_DOI.csv.gz | 	Europe PMC dataset (compressed) | 	 https://europepmc.org/pub/databases/pmc/DOI/PMID_PMCID_DOI.csv.gz                                                                                |
| PMID_PMCID_DOI.csv | 	Unzipped version of the Europe PMC dataset | 	 https://ftp.ebi.ac.uk/pub/databases/pmc/DOI/PMID_PMCID_DOI.csv                                                                                   |
| assembly_ftp | 	PDB assemblies | 	 https://ftp.ebi.ac.uk/pub/databases/msd/assemblies/split/                                                                                        |
| BLAST_DB | 	UniProt BLAST database | 	 https://ftp.uniprot.org/pub/databases/uniprot/uniprot_sprot/uniprot_sprot.fasta.gz                                                               |    
| sifts_GO | 	PDB chain Gene Ontology mapping | 	 https://ftp.ebi.ac.uk/pub/databases/msd/sifts/pdb_chain_go.csv                                                                                   |    
| GO_obo | 	Gene Ontology definitions | 	 https://current.geneontology.org/ontology/go.obo                                                                                                 |    
| GO_interpro | 	InterPro to GO mapping | 	 https://ftp.ebi.ac.uk/pub/databases/GO/goa/external2go/interpro2go                                                                               |    
| sifts | 	SIFTS data | 	 https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/                                                                                         |    
| alphafold_ftp | 	AlphaFold DB accession IDs | 	 https://ftp.ebi.ac.uk/pub/databases/alphafold/accession_ids.csv                                                                                  |    
| rfam_ftp | 	RFAM files | 	  https://www.ebi.ac.uk/pdbe/search/pdb/select?q=emdb_id:*%20AND%20rfam:%5B*%20TO%20*%5D&wt=csv&fl=emdb_id,pdb_id,rfam,rfam_id,entity_id&rows=9999999 |
| emd-xxxx-v30.xml | EMDB metadata | https://ftp.ebi.ac.uk/pub/databases/emdb/ |
| xxxxx.xml | EMPIAR metadata | https://ftp.ebi.ac.uk/pub/databases/emtest/empiar |

### Usage

To use the tools and scripts in this repository, you just need to clone it and ensure the config.ini file is properly configured as described above.

#### Executing the scripts:

Execute the scripts independently in the following recommended order:
##### EMPIAR mapping
```
fetch_empiar.py: python fetch_empiar.py -w <output_dir_to_store_annotated_empiar_files> -f <path_to_empiar_metadata_files>
```
##### Publication mapping
```
fetch_pubmed.py: python fetch_pubmed.py -w <output_dir_to_store_annotated_pubmed_files> -f <path_to_emdb_metadata_files>
```
##### Protein, complexes and ligands mapping
```
added_annotations.py: python added_annotations.py -w <output_dir_to_store_added_annotations> -f <path_to_emdb_metadata_files> --all -t <number_of_threads>
```
##### AlphaFold DB mapping
```
fetch_afdb.py: python fetch_afdb.py -w <output_dir_to_store_annotated_alphafdb_files>
```
##### Write files
```
write_xml.py: python write_xml.py <output_dir_to_store_EMICSS_xml_files>
```

### Further information

For more information about EMICSS, visit the official EMICSS website (https://www.ebi.ac.uk/emdb/emicss). This page provides detailed information about the EMDB/EMICSS project.
