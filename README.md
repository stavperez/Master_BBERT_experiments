**Project Overview**

This repository contains the full codebase and data pipeline supporting the BBERT project, developed as part of my M.Sc. research project at the Hebrew University of Jerusalem. The goal of this work is to enable fast and accurate domain-level classification of short metagenomic reads-specifically, distinguishing whether a read originates from a bacterial or eukaryotic organism.

To achieve this, we introduce BBERT, a domain-classification model built by fine-tuning a BERT-style language model on DNA sequences. This repository includes the tools for generating the training and evaluation datasets, simulating reads from real RefSeq genomes, training the BBERT model, and evaluating its performance and interpretability.

**Related Publication**

Fast and accurate taxonomic domain assignment of short metagenomic reads using BBERT
bioRxiv preprint, 2025.09.07.674730v1

This repository provides the code and data used in the publication, including all scripts for genome processing, read simulation, model training, and analysis.

**Main Components**

* Scripts for genome collection, taxonomic filtering, and synthetic read simulation from NCBI RefSeq

* A modified version of InSilicoSeq tailored for read generation with custom per-domain parameters

* Python code for BBERT model training and evaluation, including tokenization and classification pipelines

* Tools for analyzing model predictions, visualizing embedding space, and assessing domain-level generalization

**TL;DR**

* The Data Generation Pipeline builds a high-quality, taxonomically annotated synthetic dataset from RefSeq genomes, simulates millions of short reads from bacterial and eukaryotic organisms, and uses these reads to train a BERT-based classifier (BBERT) to assign each read to its correct domain. The result is a fast and effective tool for large-scale preprocessing of metagenomic data.

**Data Generation Pipeline**

This pipeline constructs the synthetic dataset used to train and evaluate BBERT. It spans genome selection, filtering, downloading, read simulation, and dataset construction.

1. Collect and Annotate RefSeq Assemblies
   
   * Script: all_refseq_genomes.py
   * Downloads assembly_summary.txt from NCBI RefSeq for multiple domains: bacteria, archaea, fungi, protozoa, plants, invertebrates, and vertebrates
   * Filters for entries labeled “chromosome” or “complete genome”
   * Fetches full taxonomic lineage using NCBI Entrez
   * Outputs: all_assemblies.csv (raw data), chromosome_assemblies_with_taxonomy.csv (with taxonomy)

2. Flag Training Set Overlap (Bacteria)
   * Script: add_train_flag.py
   * Adds a boolean in_train column to bacterial genomes by checking if their GCF/GCA accession appears in a predefined training set
   * Output: CSV with an added in_train column

3. Exclude Training Genera and Select Unique Genomes
   * Script: exclude_train_unique_genus.py
   * Removes all genomes from genera that appear in the training set
   * From the remaining set, selects one representative genome per genus
   * Output: chromosome_assemblies_unique_genus_non_train.csv
     
4. Download Genomes and Annotations
   * Script: genomes_refseq_download.py
   * Downloads genomic FASTA (.fna.gz) and annotation files (.gtf.gz or .gff.gz) for each selected assembly
   * Files are saved under: genomes_from_ncbi_unique_genus/<GCF_ID>/. Note: Archaea are skipped by default in this pipeline.

5. Attach File Paths and Compute Genome Stats
   * Script: add_paths_and_stats.py
   * Adds full paths for each genome and annotation file
   * Computes basic statistics such as genome length and sequence count
   * Adds a category column based on domain (bacteria, eukaryote, archaea)
   * Outputs: chromosome_assemblies_unique_genus_with_paths.csv
   * Subsets: genomes_bacteria_unique_genus.csv, genomes_eukaryotes_unique_genus.csv (attached)

6. Simulate Reads from Each Genome
   * Scripts: InSilicoSeq/run_read_generator.sh → calls InSilicoSeq/read_generator_one_org.py → uses modified iss/generator.py
   * Generates paired-end reads per genome using a customized version of InSilicoSeq. Parameters: 35,000 reads per direction
   * Output: Per-organism FASTQ files (R1/R2), labeled by source

7. Create Mixed Datasets
   * Script: mix_reads_new.py
   * Combines per-organism read files into: mixed_all.fastq, mixed_bacteria.fastq, mixed_eukaryotes.fastq
   * Datasets are balanced and ready for BBERT training and evaluation
