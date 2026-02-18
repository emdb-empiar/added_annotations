# Added annotations: EMICSS (**E**MDB **In**tegration with **C**omplexes, **S**tructures and **S**equences)

This repository provides tools and scripts for extracting and adding annotations to EMDB entries, which are used to enhance the metadata associated with EM datasets.

### Table of Contents

* Installation
* Docker Installation
* Configuration
* Usage
* Docker Usage
* Contributing
* License

### Installation

To install the necessary dependencies, run: 
pip install -r requirements.txt

### Docker Installation

You can also run the scripts using Docker, which provides a containerized environment with all dependencies pre-installed.

#### Building the Docker Image

```bash
docker build -t added-annotations .
```

This will create a Docker image with Python 3.8, BLAST+, and all required Python dependencies.

### Configuration

The repository uses a config.ini file for configuration, which is not included in the repository. This file should be created in the root directory of the project with the following structure:

```
[file_paths]
uniprot_tab = <path_to_file>/uniprot.tsv
CP_ftp = <path_to_file>/complextab
components_cif = <path_to_file>/components.cif
chem_comp_list = <path_to_file>/chem_comp_list.xml
pmc_ftp_gz = <path_to_file>/PMID_PMCID_DOI.csv.gz
pmc_ftp = <path_to_file>/PMID_PMCID_DOI.csv
emdb_pubmed = <path_to_file>/emdb_pubmed.log
emdb_orcid = <path_to_file>/emdb_orcid.log
assembly_ftp = <path_to_file>/assembly/
BLAST_DB = <path_to_file>/ncbi-blast-2.13.0+/database/uniprot_sprot
BLASTP_BIN = blastp
sifts_GO = <path_to_file>/pdb_chain_go.csv
GO_obo = <path_to_file>/go.obo
GO_interpro = /nfs/ftp/pub/databases/GO/goa/external2go/interpro2go
sifts = <path_to_file>/split_xml/
alphafold_ftp = <path_to_file>/accession_ids.txt
rfam_ftp = <path_to_file>/rfam_files_combined.txt

[api]
pmc = https://www.ebi.ac.uk/europepmc/webservices/rest/searchPOST
```

#### Docker Configuration

When using Docker, the config.ini file should use container paths. An example configuration file is provided in `config.ini.docker-example`. Create your own config file on the host machine with the following structure:

```
[file_paths]
CP_ftp = /data/cpx/
components_cif = /data/components.cif
pmc_ftp_gz = /data/pmc/PMID_PMCID_DOI.csv.gz
pmc_ftp = /data/pmc/PMID_PMCID_DOI.csv
assembly_ftp = /data/pdbe/assembly/
BLAST_DB = /data/uniprotkb_swissprot
BLASTP_BIN = blastp
sifts_GO = /data/pdbe/go/pdb_chain_go.csv
GO_obo = /data/go.obo
emdb_empiar_list = /data/emdb_empiar.json
sifts = /data/sifts/
alphafold_ftp = /data/accession_ids.txt
uniprot_tab = /data/uniprot.tsv

[api]
pmc = https://www.ebi.ac.uk/europepmc/webservices/rest/searchPOST

[params]
minimal_map_fragment_length = 15
```

**Note:** The paths in the Docker config should match the container mount points (e.g., `/data/...`), not the host paths.

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

### Docker Usage

When running the scripts in Docker, you need to mount your data directories and config file as read-only volumes. The general pattern is:

```bash
docker run --rm \
  -v /path/on/host/config.ini:/config/config.ini:ro \
  -v /path/on/host/data:/data:ro \
  -v /path/on/host/output:/output \
  added-annotations python <script_name.py> <arguments>
```

#### Docker Volume Mounting

- `-v /path/on/host/config.ini:/config/config.ini:ro` - Mount your config file as read-only
- `-v /path/on/host/data:/data:ro` - Mount your data directory containing all required files (cpx, components.cif, etc.) as read-only
- `-v /path/on/host/output:/output` - Mount output directory for writing results (read-write)

**Important:** 
- Use `:ro` flag for read-only mounts on data and config to prevent accidental modifications
- Ensure your config.ini uses container paths (e.g., `/data/...`) that match your volume mounts
- Map all directories referenced in your config.ini file to appropriate container paths

#### Running Scripts in Docker

Execute the scripts independently in the following recommended order:

##### EMPIAR mapping
```bash
docker run --rm \
  -v /path/on/host/config.ini:/config/config.ini:ro \
  -v /path/on/host/empiar_metadata:/empiar_metadata:ro \
  -v /path/on/host/output:/output \
  added-annotations python fetch_empiar.py -w /output -f /empiar_metadata
```

##### Publication mapping
```bash
docker run --rm \
  -v /path/on/host/config.ini:/config/config.ini:ro \
  -v /path/on/host/emdb_metadata:/emdb_metadata:ro \
  -v /path/on/host/output:/output \
  added-annotations python fetch_pubmed.py -w /output -f /emdb_metadata
```

##### Protein, complexes and ligands mapping
```bash
docker run --rm \
  -v /path/on/host/config.ini:/config/config.ini:ro \
  -v /path/on/host/data:/data:ro \
  -v /path/on/host/emdb_metadata:/emdb_metadata:ro \
  -v /path/on/host/output:/output \
  added-annotations python AddedAnnotations.py -w /output -f /emdb_metadata --all -t 4
```

##### AlphaFold DB mapping
```bash
docker run --rm \
  -v /path/on/host/config.ini:/config/config.ini:ro \
  -v /path/on/host/data:/data:ro \
  -v /path/on/host/output:/output \
  added-annotations python fetch_afdb.py -w /output
```

##### Generate Europe PMC Links
```bash
docker run --rm \
  -v /path/on/host/config.ini:/config/config.ini:ro \
  -v /path/on/host/output:/output \
  added-annotations python generate_eupmc_links.py
```

##### Compare Release
```bash
docker run --rm \
  -v /path/on/host/config.ini:/config/config.ini:ro \
  -v /path/on/host/latest:/latest:ro \
  -v /path/on/host/previous:/previous:ro \
  added-annotations python compare_release.py /latest /previous
```

##### Write XML files
```bash
docker run --rm \
  -v /path/on/host/output:/output \
  added-annotations python write_xml.py /output
```

### Further information

For more information about EMICSS, visit the official EMICSS website (https://www.ebi.ac.uk/emdb/emicss). This page provides detailed information about the EMDB/EMICSS project.
