**Project Overview**

This repository contains the core codebase and scripts developed as part of my M.Sc. research at the Hebrew University of Jerusalem, focusing on LLM-based metagenomic analysis. The project explores the use of large language models (LLMs) for analyzing DNA sequencing data at scale, enabling functional interpretation, classification, and visualization of metagenomic reads across diverse organisms.

**Related Publication**

For more details on the methodology and results, see our paper:
"Fast and accurate taxonomic domain assignment of short metagenomic reads using BBERT"
Available on bioRxiv: https://www.biorxiv.org/cgi/content/short/2025.09.07.674730v1

**Main Components**

* Python scripts for read generation, annotation, classification, and statistical analysis
* Modeling pipelines that adapt LLM architectures (e.g., BERT-style models) to biological sequences
* Visualization tools for embeddings, confusion types (TP/FN), and model interpretability
* Custom datasets generated from NCBI RefSeq genomes, including bacterial and eukaryotic organisms

**TL;DR**

The Data Generation pipeline programmatically collects RefSeq chromosome-level genomes across diverse taxa, enriches them with taxonomy, filters out genera used in training, downloads genome and annotation files, and attaches file paths and metadata. It then simulates organism-specific reads using a modified version of InSilicoSeq and combines the outputs into mixed datasets (all-organism, bacteria-only, eukaryote-only) for LLM-based metagenomic analysis.


**Data Generation Pipeline**

This pipeline builds a large, taxonomically diverse genome set from NCBI RefSeq, simulates reads using a customized version of InSilicoSeq, and produces mixed datasets for downstream BBERT analysis.

Overview of Pipeline Steps
1. Collect & Annotate RefSeq Assemblies

Script: all_refseq_genomes.py

Downloads assembly_summary.txt from NCBI RefSeq across multiple domains: bacteria, archaea, fungi, protozoa, plants, invertebrates, and vertebrates.

Filters to include only chromosome or complete genome entries.

Fetches taxonomic lineage via NCBI Entrez.

Outputs:

all_assemblies.csv (raw)

chromosome_assemblies_with_taxonomy.csv (with lineage)

2. Flag Training Set Overlap (Bacteria)

Script: add_train_flag.py

Adds an in_train flag for bacterial genomes based on whether their GCF/GCA IDs appear in a reference training abundance file.

Output: CSV with in_train column added.

3. Exclude Training Genera & Select Unique Genus

Script: exclude_train_unique_genus.py

Excludes any genus that appears in the training set (in_train=True).

From the remaining set, keeps one representative genome per genus.

Output: chromosome_assemblies_unique_genus_non_train.csv

4. Download Genomes & Annotations

Script: genomes_refseq_download.py

Downloads *_genomic.fna.gz and annotation files (.gtf preferred, .gff as fallback) for selected assemblies.

Skips archaea if needed.

Genomes are stored in subfolders under genomes_from_ncbi_unique_genus/<GCF_ID>/.

5. Attach Paths & Compute Stats

Script: add_paths_and_stats.py

Adds absolute paths for each genome and annotation file.

Computes total genome length, number of sequences, and number of headers.

Categorizes organisms by domain (e.g., bacteria, eukaryotes).

Outputs:

chromosome_assemblies_unique_genus_with_paths.csv (master)

Filtered subsets: genomes_bacteria_unique_genus.csv, genomes_eukaryotes_unique_genus.csv, etc.

6. Simulate Reads

Script Chain:

InSilicoSeq/run_read_generator.sh
→ calls InSilicoSeq/read_generator_one_org.py
→ uses modified iss/generator.py

Simulates paired-end reads for each genome using InSilicoSeq, modified for project-specific needs.

Parameters: 35,000 reads per direction

Outputs are per-organism paired FASTQ files (R1/R2).

7. Create Mixed Datasets

Script: mix_reads_new.py

Combines reads across organisms into final mixed datasets:

All organisms, Bacteria-only, Eukaryote-only

Datasets are structured and labeled for use in downstream classification tasks.
