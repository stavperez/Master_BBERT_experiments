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
* The Assembly & Alignment Comparison Pipeline evaluates the impact of BBERT‑classified reads on downstream genome assembly and alignment. It assembles the filtered reads into contigs, aligns them back to high‑quality reference genomes, and computes genome‑level coverage, mismatch, and gap statistics. This produces clear, quantitative and visual comparisons of contig quality and error profiles between BBERT‑processed and unfiltered data, demonstrating the model’s effect on assembly accuracy and reliability.
* The Model Error Analysis & Embedding Exploration Pipeline investigates where BBERT succeeds and fails, and why. It combines multiple analytical directions:
  * Confusion Type Annotation: Each read is matched to its true label and evaluated against BBERT’s prediction, producing a full per-read confusion classification (TP, FP, FN, TN) used throughout downstream analysis.
  * Dimensional Reduction - Embeddings Visualization: BBERT’s internal embeddings are projected into 2D space using PCA and t-SNE, revealing structure in how the model organizes reads. These visualizations expose clusters by organism, function, domain, and error type.
  * Phylogenetic Bias Analysis: Model errors are analyzed across taxonomic groups. False negative rates (FNR) are computed by phylum, class, and order, and tested for statistical significance to detect systematic taxonomic biases.
  * Functional Annotation & Metadata Enrichment: Each read is enriched with biological context—such as plasmid origin, coding status, inferred reading frame, and gene function—by integrating genome annotations and BLASTX results.
  * Functional Error Enrichment (TP vs FN): GO terms and product annotations are compared between true positives and false negatives to identify specific functions that are systematically harder for the model. This reveals model weaknesses at the functional level.

<details> <summary><strong>Data Generation Pipeline</strong></summary>

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
</details>

<details> <summary><strong>Assembly & Alignment Comparison Pipeline</strong></summary>

This section describes the pipeline used to evaluate BBERT’s effect on downstream assembly and alignment quality. We compare contigs generated from BBERT-classified reads to those from unfiltered reads, analyzing coverage, error rates, and alignment patterns. The methodology is described in more detail in our paper.

1. Run BBERT on Mixed Reads
   * Script: run_score_256_orgs.sh
   * BBERT is applied to paired-end reads from 256 diverse organisms, producing domain classification scores (bacterial vs eukaryotic) per read.

3. Classify and Split Reads
   * Scripts: run_split_fastq_new.sh → split_fastq_by_score_new.py
   * Reads are split into bacterial and eukaryotic subsets based on BBERT predictions, using both strands for robust labeling.

4. Assemble Classified Reads
   * Script: run_spades_assembly_new.sh
   * SPAdes is used to assemble the classified bacterial reads into contigs.

5. Prepare Reference Database & Align Contigs
   * Scripts:
     * create_master_blast_db.sh - builds a BLAST DB from 512 bacterial reference genomes
     * blast_alignment_assembly_contigs_new.sh – aligns assembled contigs back to the references
     * blast_error_alignment_counts.sh – annotates alignments with detailed error types (mismatches, gaps, runs)

6. Compare Assemblies
   * Scripts:
     * compare_contigs_pairwise.py - compares contigs from BBERT-classified vs unfiltered assemblies
     * compare_genome_covrage_mismatches_gaps.py – computes genome-level stats: coverage, error rates, error distributions
     * Results include: coverage gains/losses, error type ratios, and contig-specific differences

7. Visualize Comparison Results
   * Scripts:
     * plot_alinment_comparison.py - generates genome-wide scatterplots and histograms of coverage and error improvements
     * plot_single_org_alignment_windows.py – zooms into specific genomic regions to visualize local alignment quality in classified vs unfiltered assemblies
</details>

<details> <summary><strong>Model Error Analysis & Embedding Exploration Pipeline</strong></summary>

This part of the project aims to understand BBERT’s predictions by analyzing where the model succeeds (True Positives) and fails (False Negatives), and why. We explore taxonomic, functional, and sequence-level patterns behind BBERT's decisions using a variety of tools: confusion matrices, dimensionality reduction, GO term enrichment, and phylogenetic profiling.

1. Confusion Type Annotation
   * Goal: Infer the true label for each read and determine whether BBERT's classification is a TP, FP, FN, or TN.
   * Scripts:
     * scripts/error_analysis_new.py - Parses the original FASTQ files to extract the true source of each read (from the Source: tag). Then compares this to BBERT’s predicted label and generates a confusion type (TP, FP, FN, TN) per read. Outputs a CSV with true label, predicted label, and confusion type.
     * local_scripts/confusion_matrix_for_bbert_results.py - Computes and visualizes confusion matrices, either from raw FASTQ files or from the annotated confusion results.

2. Dimensional Reduction - Embeddings Visualization
   * This stage analyzes BBERT’s internal representations by projecting high-dimensional embeddings into 2D space, enabling visualization of patterns across organisms, functions, and error types. The reduced embeddings are used to interpret how reads are grouped or separated by the model.
   * Scripts:
     * scripts/to_embeddings.py, scripts/to_embeddings_2.py - Prepare balanced FASTQ subsets for embedding analysis:
       1. Equal TP/FN reads per reading frame
       2. Balanced sets for coding vs. non-coding reads
     * local_scripts/plot_embeddings/plot_tSNE.py - Visualizes embeddings from individual organisms using PCA and t-SNE, Highlights how reads cluster within a single genome.
     * local_scripts/plot_embeddings/plot_embeddings_many_organisms.py - Projects and compares embeddings across many pairs of organisms, generating PCA and t-SNE plots. Color-coded by organism and domain to reveal inter-species and inter-domain separability.
     * local_scripts/plot_embeddings/plot_embeddings_go_terms.py - Merges embeddings with GO Slim categories, confusion labels, and reading frame annotations.
     * Generates:
       1. PCA/t-SNE plots by GO Slim category
       2. Plots colored by TP/FN confusion type
       3. Joint visualization of confusion type + reading frame using colors and markers
     * local_scripts/plot_embeddings/data_analyze.py - Analyzes how confusion types are distributed across reading frames and coding status. Provides a 3D breakdown (CDS × Confusion × Frame) to quantify model weaknesses.

3. Phylogenetic Bias Analysis
   * Goal: Quantify how BBERT’s error rates (especially FN) vary across taxonomic groups.
   * Scripts:
     * local_scripts/phylogenetic_analysis/plot_phylogenetic_tree.py - Builds a phylogenetic tree from species taxonomies using ETE3 and overlays confusion statistics (e.g., FNR) as colored highlights. Can export trees to Newick and iTOL formats with labeled clades and gradients.
     * local_scripts/phylogenetic_analysis/bias_analysis.py - Computes False Negative Rate (FNR) per taxonomic level (phylum, class, etc.), and performs Shapiro tests for normality, Kruskal-Wallis tests to detect differences between groups, Dunn’s post-hoc tests to find significantly biased taxonomic clades

4. Functional Annotation & Enrichment
   * Goal: Add biological context (e.g. plasmid origin, GO annotations, coding status) to each read, and then analyze if errors correlate with function.
   * Scripts:
     * scripts/add_plasmid_info.py → add_plasmid_info_map_creation.py - Scans .fna.gz files to identify plasmid vs. chromosome sequences per organism. Builds a dictionary and annotates each read as is_plasmid.
     * scripts/add_reading_frame_from_blastx.py - Converts FASTQ → FASTA, aligns against known protein sequences using BLASTX, and assigns a likely reading frame to each read based on the top hit.
     * scripts/extract_read_metadata.py - Parses encoded read metadata from the FASTQ headers (start, end, CDS, etc.) and saves structured CSVs for downstream enrichment.
     * scripts/run_annotate_reads_array.sh → annotate_single_organism.py - For each organism:
       1. Builds an interval tree from its GTF annotation
       2. Matches each read to the closest gene/CDS using coordinates
       3. Assigns gene/function/CDS overlap to each read
     * scripts/aggregate_and_merge_functionality_results.py - Combines all per-organism annotations and merges with BBERT prediction results. Output is a unified dataset with classification + metadata per read.
     * scripts/generalize_GO_ontology.py -
       1. Converts specific GO terms into generalized categories using the GO DAG (up to depth 3)
       2. Maps GO terms to human-readable names and to GO Slim categories
       3. Adds GO_terms_generalized, GO_names_generalized, and GO_slim_categories columns

5. Functional Error Enrichment (TP vs FN)
   * Goal: Identify GO terms or product annotations that are significantly enriched among BBERT’s errors (FN), using fold-change and statistical tests.
   * Scripts:
     * local_scripts/functionality_go_analysis/analyze_functionality.py - For each functional field (GO_terms, product, etc.):
       1. Counts TP vs. FN reads per term
       2. Computes normalized proportions
       3. Calculates log2 fold-change (FN/TP)
       4. Applies Fisher’s exact test
       5. Saves volcano plot input CSVs
     * local_scripts/functionality_go_analysis/Fisher's_test.py - Loads volcano CSVs and:
       1. Adds odds ratios and corrected p-values
       2. Plots: Volcano plots (log2FC vs p-value), Bar plots of top enriched terms, Histograms of odds ratio distributions
     * local_scripts/functionality_go_analysis/plot_analyze_functionality.py - Adds Mann–Whitney U tests alongside Fisher's test. Produces comparative volcano plots for both methods to validate significance.
     * local_scripts/functionality_go_analysis/dist.py - Analyzes and visualizes the distribution of log2 fold-change for each functional category (e.g., GO term, product). Highlights top and bottom 2.5% terms that most distinguish FN from TP.

</details>
