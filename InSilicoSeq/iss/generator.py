#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gc
import logging
import os
import random
import sys
import re

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
import gffutils


from iss import abundance, download, util
from iss.error_models import basic, kde, perfect
from iss.util import load, rev_comp


def simulate_reads(
    record,
    error_model,
    n_pairs,
    cpu_number,
    forward_handle,
    reverse_handle,
    mutations_handle,
    sequence_type,
    gc_bias=False,
    mode="default",
):
    """Simulate reads from one genome (or sequence) according to an ErrorModel

    This function makes use of the `simulate_read` function to simulate reads
    and save them in a fastq file

    Args:
        record (SeqRecord): sequence or genome of reference
        error_model (ErrorModel): an ErrorModel
        n_pairs (int): the number of reads to generate
        cpu_number (int): an int indentifying the cpu that is used by the
            function. Is used for naming the output file
        forward_handle (file): a file handle to write the forward reads to
        reverse_handle (file): a file handle to write the reverse reads to
        mutations_handle (file): a file handle to write the mutations to
        sequencing_type (str): metagenomics or amplicon sequencing used
        gc_bias (bool): if set, the function may skip a read due to abnormal
            GC content

    """
    logger = logging.getLogger(__name__)

    # load the record from disk if mode is memmap
    if mode == "memmap":
        record_mmap = load(record)
        record = record_mmap

    logger.debug("Cpu #%s: Generating %s read pairs" % (cpu_number, n_pairs))

    for forward_record, reverse_record, mutations in reads_generator(
        n_pairs, record, error_model, cpu_number, gc_bias, sequence_type
    ):
        SeqIO.write(forward_record, forward_handle, "fastq-sanger")
        SeqIO.write(reverse_record, reverse_handle, "fastq-sanger")
        write_mutations(mutations, mutations_handle)


def reads_generator(n_pairs, record, error_model, cpu_number, gc_bias,
                    sequence_type, cds_regions):
    logger = logging.getLogger(__name__)

    i = 0
    while i < n_pairs:
        try:
#            forward, reverse, mutations = simulate_read(record,
#                                                        error_model, i,
#                                                        cpu_number, sequence_type, cds_regions)

            forward, reverse, mutations = simulate_read_150_length(record,
                                                        error_model, i,
                                                        cpu_number,
                                                        sequence_type,
                                                        cds_regions)


        except AssertionError:
#            logger.warning("%s shorter than read length for this ErrorModel" % record.id)
#            logger.warning("Skipping %s. You will have less reads than specified" % record.id)
            break
        else:
            if gc_bias:
                stiched_seq = forward.seq + reverse.seq
                gc_content = gc_fraction(stiched_seq)
                if 40 < gc_content < 60:
                    yield (forward, reverse, mutations)
                    i += 1
                elif np.random.rand() < 0.90:
                    yield (forward, reverse, mutations)
                    i += 1
                else:
                    continue
            else:
                yield (forward, reverse, mutations)
                i += 1



def simulate_read_150_length(record, error_model, i, cpu_number,
                            sequence_type, cds_regions):
    """From a read pair from one genome (or sequence) according to an
    ErrorModel

    Each read is a SeqRecord object
    returns a tuple containing the forward and reverse read.

    Args:
        record (SeqRecord): sequence or genome of reference
        error_model (ErrorModel): an ErrorModel class
        i (int): a number identifying the read
        cpu_number (int): cpu number. Is added to the read id.
        sequence_type (str): metagenomics or amplicon sequencing used

    Returns:
        tuple: tuple containg a forward read and a reverse read
    """
    logger = logging.getLogger(__name__)
    sequence = record.seq
    seq_len = len(sequence)
    header = record.id
    description = record.description
    chromosome = 0
    match = re.search(r"chromosome\s+(\d+)", description, re.IGNORECASE)
    if match:
        chromosome = match.group(1)

    read_length = error_model.read_length
    # read_length = 150

    if error_model.fragment_length is not None and error_model.fragment_sd is not None:
        fragment_length = int(np.random.normal(error_model.fragment_length, error_model.fragment_sd))
        insert_size = fragment_length - (read_length * 2)
    else:
        insert_size = error_model.random_insert_size()
        fragment_length = insert_size + (read_length * 2)

    # generate the forward read
    try:  # a ref sequence has to be longer than 2 * read_length + i_size
        assert read_length < len(record.seq)
        # assign the start position of the forward read
        # if sequence_type == metagenomics, get a random start position
        # if sequence_type == amplicon, start position is the start of the read
        if sequence_type == "metagenomics":
            forward_start = random.randrange(0, seq_len - fragment_length)
        elif sequence_type == "amplicon":
            forward_start = 0
        else:
            raise RuntimeError(f"sequence type '{sequence_type}' is not supported")
    except AssertionError:
        raise
    except ValueError as e:
        logger.debug("%s shorter than template length for this ErrorModel:%s" % (record.id, e))
        forward_start = max(0, random.randrange(0, len(record.seq) - read_length))

    forward_end = forward_start + read_length

    bounds = (forward_start, forward_end)

    forward_substring = sequence[forward_start:forward_end]
    forward_read_info = {"start": forward_start, "end": forward_end,
                         "strand": "+",
                 "chrom": str(header) + "_" + str(chromosome)}

    forward_read = Seq(str(sequence[forward_start:forward_end]))
    forward_info = update_read_description2(Seq(forward_substring),
                                           forward_read_info,
                                             cds_regions)

    # create a perfect read
    forward = SeqRecord(
        forward_read, id="%s_%s_%s/1" % (
            header, i, cpu_number), description=forward_info
    )
    forward.annotations["mutations"] = []
    forward.annotations["original"] = str(forward.seq)

    # add the indels, the qual scores and modify the record accordingly
    forward = error_model.introduce_indels(forward, "forward", sequence, bounds)
    forward = error_model.introduce_error_scores(forward, "forward")
    forward = error_model.mut_sequence(forward, "forward")

    trim_start = 75
    trim_end = len(forward.seq) - 75

    if trim_end > trim_start:
        # Extract and trim the quality scores if they exist
        if "phred_quality" in forward.letter_annotations:
            trimmed_phred_quality = forward.letter_annotations[
                                        "phred_quality"][trim_start:trim_end]
        else:
            trimmed_phred_quality = None  # No quality scores to trim

        # Clear all existing letter annotations to allow sequence modification
        forward.letter_annotations = {}

        # Trim the sequence
        forward.seq = Seq(str(forward.seq[trim_start:trim_end]))

        # Restore the trimmed quality scores (if they exist)
        if trimmed_phred_quality is not None:
            forward.letter_annotations[
                "phred_quality"] = trimmed_phred_quality
        # Update "start" and "end" in the description
        forward_read_info["start"] += trim_start
        forward_read_info["end"] -= 75
        forward.description = update_read_description2(forward.seq,
                                                       forward_read_info,
                                                       cds_regions)

    # generate the reverse read
    # assign start position reverse read
    # if sequence_type == metagenomics, get a start position based on insert_size
    # if sequence_type == amplicon, start position is the end of the read
    if sequence_type == "metagenomics":
        reverse_start = forward_end + insert_size
        reverse_end = reverse_start + read_length
    elif sequence_type == "amplicon":
        reverse_start = seq_len - read_length
        reverse_end = reverse_start + read_length
    else:
        raise ValueError(f"Sequence type {sequence_type} not known")
    if reverse_end > len(record.seq):
        # we use random insert when the modelled template length distribution
        # is too large
        reverse_end = random.randrange(read_length, len(record.seq))
        reverse_start = reverse_end - read_length

    bounds = (reverse_start, reverse_end)

    reverse_read_info = {"start": reverse_start, "end": reverse_end,
                       "strand": "+",
                 "chrom": str(header) + "_" + str(chromosome)}
    reverse_substring = sequence[reverse_start:reverse_end]
    rev_comp_seq = rev_comp(str(reverse_substring))

    backward_info = update_read_description2(Seq(rev_comp_seq),
                                             reverse_read_info,
                                             cds_regions)
    # create a perfect read
    reverse = SeqRecord(
        Seq(rev_comp_seq),
        id="%s_%s_%s/2" % (header, i, cpu_number),
        description=backward_info,
    )
    reverse.annotations["mutations"] = []
    reverse.annotations["original"] = str(reverse.seq)

    # add the indels, the qual scores and modify the record accordingly
    reverse = error_model.introduce_indels(reverse, "reverse", sequence, bounds)
    reverse = error_model.introduce_error_scores(reverse, "reverse")
    reverse = error_model.mut_sequence(reverse, "reverse")
    if trim_end > trim_start:
        # Extract and trim the quality scores if they exist
        if "phred_quality" in reverse.letter_annotations:
            trimmed_phred_quality_reverse = reverse.letter_annotations[
                                                "phred_quality"][
                                            trim_start:trim_end]
        else:
            trimmed_phred_quality_reverse = None  # No quality scores to trim

        # Clear all existing letter annotations to allow sequence modification
        reverse.letter_annotations = {}

        # Trim the sequence
        reverse.seq = Seq(str(reverse.seq[trim_start:trim_end]))

        # Restore the trimmed quality scores (if they exist)
        if trimmed_phred_quality_reverse is not None:
            reverse.letter_annotations[
                "phred_quality"] = trimmed_phred_quality_reverse

        # Update read description AFTER all modifications
        reverse_read_info["start"] += trim_start
        reverse_read_info["end"] -= 75
        updated_reverse_description = update_read_description2(reverse.seq,
                                                               reverse_read_info,
                                                               cds_regions)

        # Update the reverse read with the new description
        reverse.description = updated_reverse_description

    combined_mutations = forward.annotations["mutations"] + reverse.annotations["mutations"]

    return (forward, reverse, combined_mutations)



def simulate_read(record, error_model, i, cpu_number, sequence_type, cds_regions):
    """From a read pair from one genome (or sequence) according to an
    ErrorModel

    Each read is a SeqRecord object
    returns a tuple containing the forward and reverse read.

    Args:
        record (SeqRecord): sequence or genome of reference
        error_model (ErrorModel): an ErrorModel class
        i (int): a number identifying the read
        cpu_number (int): cpu number. Is added to the read id.
        sequence_type (str): metagenomics or amplicon sequencing used

    Returns:
        tuple: tuple containg a forward read and a reverse read
    """
    logger = logging.getLogger(__name__)
    sequence = record.seq
    seq_len = len(sequence)
    header = record.id
    description = record.description
    chromosome = 0
    match = re.search(r"chromosome\s+(\d+)", description, re.IGNORECASE)
    if match:
        chromosome = match.group(1)

    read_length = error_model.read_length
    # read_length = 150

    if error_model.fragment_length is not None and error_model.fragment_sd is not None:
        fragment_length = int(np.random.normal(error_model.fragment_length, error_model.fragment_sd))
        insert_size = fragment_length - (read_length * 2)
    else:
        insert_size = error_model.random_insert_size()
        fragment_length = insert_size + (read_length * 2)

    # generate the forward read
    try:  # a ref sequence has to be longer than 2 * read_length + i_size
        assert read_length < len(record.seq)
        # assign the start position of the forward read
        # if sequence_type == metagenomics, get a random start position
        # if sequence_type == amplicon, start position is the start of the read
        if sequence_type == "metagenomics":
            forward_start = random.randrange(0, seq_len - fragment_length)
        elif sequence_type == "amplicon":
            forward_start = 0
        else:
            raise RuntimeError(f"sequence type '{sequence_type}' is not supported")
    except AssertionError:
        raise
    except ValueError as e:
        logger.debug("%s shorter than template length for this ErrorModel:%s" % (record.id, e))
        forward_start = max(0, random.randrange(0, len(record.seq) - read_length))

    forward_end = forward_start + read_length
    bounds = (forward_start, forward_end)

    forward_substring = sequence[forward_start:forward_end]
    forward_read_info = {"start": forward_start, "end": forward_end,
                         "strand": "+",
                 "chrom": str(header) + "_" + str(chromosome)}

    forward_read = Seq(str(sequence[forward_start:forward_end]))
    forward_info = update_read_description2(Seq(forward_substring),
                                           forward_read_info,
                                             cds_regions)

    # create a perfect read
    forward = SeqRecord(
        forward_read, id="%s_%s_%s/1" % (
            header, i, cpu_number), description=forward_info
    )
    forward.annotations["mutations"] = []
    forward.annotations["original"] = str(forward.seq)

    # add the indels, the qual scores and modify the record accordingly
    forward = error_model.introduce_indels(forward, "forward", sequence, bounds)
    forward = error_model.introduce_error_scores(forward, "forward")
    forward = error_model.mut_sequence(forward, "forward")

    # generate the reverse read
    # assign start position reverse read
    # if sequence_type == metagenomics, get a start position based on insert_size
    # if sequence_type == amplicon, start position is the end of the read
    if sequence_type == "metagenomics":
        reverse_start = forward_end + insert_size
        reverse_end = reverse_start + read_length
    elif sequence_type == "amplicon":
        reverse_start = seq_len - read_length
        reverse_end = reverse_start + read_length
    else:
        raise ValueError(f"Sequence type {sequence_type} not known")
    if reverse_end > len(record.seq):
        # we use random insert when the modelled template length distribution
        # is too large
        reverse_end = random.randrange(read_length, len(record.seq))
        reverse_start = reverse_end - read_length
    bounds = (reverse_start, reverse_end)

    backward_read_info = {"start": reverse_start, "end": reverse_end,
                       "strand": "+",
                 "chrom": str(header) + "_" + str(chromosome)}
    reverse_substring = sequence[reverse_start:reverse_end]
    rev_comp_seq = rev_comp(str(reverse_substring))

    backward_info = update_read_description2(Seq(rev_comp_seq),
                                             backward_read_info,
                                             cds_regions)
    # create a perfect read
    reverse = SeqRecord(
        Seq(rev_comp_seq),
        id="%s_%s_%s/2" % (header, i, cpu_number),
        description=backward_info,
    )
    reverse.annotations["mutations"] = []
    reverse.annotations["original"] = str(reverse.seq)

    # add the indels, the qual scores and modify the record accordingly
    reverse = error_model.introduce_indels(reverse, "reverse", sequence, bounds)
    reverse = error_model.introduce_error_scores(reverse, "reverse")
    reverse = error_model.mut_sequence(reverse, "reverse")
    combined_mutations = forward.annotations["mutations"] + reverse.annotations["mutations"]

    return (forward, reverse, combined_mutations)


def to_fastq(generator, output):
    """Write reads to a fastq file

    Take a generator or a list containing read pairs (tuples) and write them
        in two fastq files: output_R1.fastq and output_R2.fastq

    Args:
        generator (generator): a read generator (or list)
        output (string): the output files prefix
    """
    logger = logging.getLogger(__name__)
    # define name of output files
    output_forward = output + "_R1.fastq"
    output_reverse = output + "_R2.fastq"

    try:
        f = open(output_forward, "a")
        r = open(output_reverse, "a")
    except PermissionError as e:
        logger.error("Failed to open output file(s): %s" % e)
        sys.exit(1)
    else:
        with f, r:
            for read_tuple in generator:
                SeqIO.write(read_tuple[0], f, "fastq-sanger")
                SeqIO.write(read_tuple[1], r, "fastq-sanger")


def worker_iterator(work, error_model, cpu_number, worker_prefix, seed, sequence_type, gc_bias):
    """A utility function to run the reads simulation of each record in a loop for a specific cpu"""
    logger = logging.getLogger(__name__)
    try:
        forward_handle = open(f"{worker_prefix}_R1.fastq", "w")
        reverse_handle = open(f"{worker_prefix}_R2.fastq", "w")
        mutation_handle = open(f"{worker_prefix}.vcf", "w")
    except PermissionError as e:
        logger.error("Failed to write temporary output file(s): %s" % e)
        sys.exit(1)

    if seed is not None:
        random.seed(seed + cpu_number)
        np.random.seed(seed + cpu_number)

    with forward_handle, reverse_handle, mutation_handle:
        for record, n_pairs, mode in work:
            simulate_reads(
                record=record,
                n_pairs=n_pairs,
                mode=mode,
                error_model=error_model,
                cpu_number=cpu_number,
                forward_handle=forward_handle,
                reverse_handle=reverse_handle,
                mutations_handle=mutation_handle,
                sequence_type=sequence_type,
                gc_bias=gc_bias,
            )


def generate_work_divider(
    fasta_file, readcount_dic, abundance_dic, n_reads, coverage, coverage_file, error_model, output, chunk_size
):
    """Yields a list of tuples containing the records and the number of reads to generate for each record

    Yields:
        list[tuple[SeqIO.Record, int, str]]: a list of tuples containing the records,
            the number of reads to generate for each record and the mode
    """
    logger = logging.getLogger(__name__)

    current_chunk = 0
    total_reads_generated = 0
    total_reads_generated_unrounded = 0

    chunk_work = []

    for record in fasta_file:
        # generate reads for records
        if readcount_dic is not None:
            if record.id not in readcount_dic:
                logger.warning(f"Record {record.id} not found in readcount file")
                continue
            n_pairs_unrounded = readcount_dic[record.id] / 2
        elif abundance_dic is not None:
            if record.id not in abundance_dic:
                logger.warning(f"Record {record.id} not found in abundance file")
                continue
            record_abundance = abundance_dic[record.id]

            genome_size = len(record.seq)

            if coverage or coverage_file:
                record_coverage = record_abundance
            else:
                record_coverage = abundance.to_coverage(
                    n_reads,
                    record_abundance,
                    error_model.read_length,
                    genome_size,
                )
            n_pairs_unrounded = ((record_coverage * len(record.seq)) / error_model.read_length) / 2
        else:
            raise RuntimeError("No readcount or abundance file provided")

        # check that the rounding does not cause to drop read pairs
        n_pairs = round(n_pairs_unrounded)
        total_reads_generated_unrounded += n_pairs_unrounded
        total_reads_generated += n_pairs
        if round(total_reads_generated_unrounded) > total_reads_generated:
            logger.debug("Adding a pair to correct rounding error")
            n_pairs += 1
            total_reads_generated += 1

        logger.debug("Will generate %s read pairs for %s" % (n_pairs, record.id))

        if n_pairs == 0:
            continue

        # due to a bug in multiprocessing
        # https://bugs.python.org/issue17560
        # we can't send records taking more than 2**31 bytes
        # through serialisation.
        # In those cases we use memmapping
        if sys.getsizeof(str(record.seq)) >= 2**31 - 1:
            logger.warning("record %s unusually big." % record.id)
            logger.warning("Using a memory map.")
            mode = "memmap"

            record_mmap = "%s.memmap" % output
            if os.path.exists(record_mmap):
                os.unlink(record_mmap)
            util.dump(record, record_mmap)
            del record
            record = record_mmap
            gc.collect()
        else:
            mode = "default"

        n_pairs_remaining = n_pairs
        while n_pairs_remaining > 0:
            chunk_remaining = chunk_size - current_chunk

            if n_pairs_remaining <= chunk_remaining:
                # Record fits in the current chunk
                chunk_work.append((record, n_pairs_remaining, mode))
                n_pairs_added = n_pairs_remaining
            else:
                # Record does not fit in the current chunk
                chunk_work.append((record, chunk_remaining, mode))
                n_pairs_added = chunk_remaining

            n_pairs_remaining -= n_pairs_added
            current_chunk += n_pairs_added

            if current_chunk == chunk_size:
                yield chunk_work
                chunk_work = []
                current_chunk = 0

    if chunk_work:
        # Yield the last (not full) chunk
        yield chunk_work


def load_error_model(mode, seed, model, fragment_length, fragment_length_sd, store_mutations):
    """
    Load the error model based on the specified mode and parameters.

    Args:
        mode (str): The mode of the error model. Possible values are 'kde', 'basic', and 'perfect'.
        seed (int): The random seed to use for generating random numbers.
        model (str): The model to use for the error model. Only applicable for the 'kde' mode.
        fragment_length (float): The mean fragment length for the error model.
        fragment_length_sd (float): The standard deviation of the fragment length for the error model.

    Returns:
        err_mod: The loaded error model based on the specified mode and parameters.
    """
    logger = logging.getLogger(__name__)

    logger.info("Using %s ErrorModel" % mode)

    precomputed_error_models = {
        "hiseq": os.path.join(os.path.dirname(__file__), "profiles/HiSeq"),
        "novaseq": os.path.join(os.path.dirname(__file__), "profiles/NovaSeq"),
        "miseq": os.path.join(os.path.dirname(__file__), "profiles/miSeq_0.npz"),
        "miseq-20": os.path.join(os.path.dirname(__file__), "profiles/miSeq_20.npz"),
        "miseq-24": os.path.join(os.path.dirname(__file__), "profiles/miSeq_24.npz"),
        "miseq-28": os.path.join(os.path.dirname(__file__), "profiles/miSeq_28.npz"),
        "miseq-32": os.path.join(os.path.dirname(__file__), "profiles/miSeq_32.npz"),
        "miseq-36": os.path.join(os.path.dirname(__file__), "profiles/miSeq_36.npz"),
        "nextseq": os.path.join(os.path.dirname(__file__), "profiles/nextSeq.npz"),
    }

    if fragment_length is not None and fragment_length_sd is not None:
        logger.info(
            f"Using custom fragment length {fragment_length} and default fragment length sd {fragment_length_sd}"
        )
    elif bool(fragment_length) ^ bool(fragment_length_sd):
        logger.error("fragment_length and fragment_length_sd must be specified together")
        sys.exit(1)

    if seed:
        logger.info("Setting random seed to %i" % seed)
        random.seed(seed)
        np.random.seed(seed)
    if mode == "kde":
        if model is None:
            logger.error("--model is required in --mode kde")
            sys.exit(1)
        elif model.lower() in precomputed_error_models:
            npz = precomputed_error_models[model.lower()]
        else:
            npz = model
        err_mod = kde.KDErrorModel(npz, fragment_length, fragment_length_sd, store_mutations)
    elif mode == "basic":
        if model is not None:
            logger.warning("--model %s will be ignored in --mode %s" % (model, mode))

        err_mod = basic.BasicErrorModel(fragment_length, fragment_length_sd, store_mutations)
    elif mode == "perfect":
        if model is not None:
            logger.warning("--model %s will be ignored in --mode %s" % (model, mode))

        err_mod = perfect.PerfectErrorModel(fragment_length, fragment_length_sd)

    return err_mod


def load_genomes(genomes, draft, ncbi, n_genomes_ncbi, output, n_genomes):
    """Load genomes from different sources and concatenate them into a single file.

    Args:
        genomes (list): List of paths to genome files.
        draft (list): List of paths to draft genome files.
        ncbi (bool): Flag indicating whether to download genomes from NCBI.
        n_genomes_ncbi (list): List of tuples specifying the NCBI genome IDs and the number of genomes to download.
        output (str): Path to the output file.
        n_genomes (int): Number of genomes to select randomly from the concatenated genome file.

    Returns:
        tuple: A tuple containing the number of genomes and the path to the concatenated genome file.
    """

    logger = logging.getLogger(__name__)
    if not (genomes or draft or ncbi):
        logger.error("One of --genomes/-g, --draft, --ncbi/-k is required")
        sys.exit(1)

    genome_files = []
    if genomes:
        genome_files.extend(genomes)

    if draft:
        genome_files.extend(draft)

    if ncbi and n_genomes_ncbi:
        util.genome_file_exists(output + "_ncbi_genomes.fasta")
        if len(*ncbi) != len(*n_genomes_ncbi):
            logger.error(
                "--ncbi and --n_genomes_ncbi of unequal lengths. \
                Aborting"
            )
            sys.exit(1)

        for g, n in zip(*ncbi, *n_genomes_ncbi):
            genomes_ncbi = download.ncbi(g, n, output + "_ncbi_genomes.fasta")
        genome_files.append(genomes_ncbi)

    if ncbi and not n_genomes_ncbi:
        logger.error("--ncbi/-k requires --n_genomes_ncbi/-U. Aborting.")
        sys.exit(1)

    genome_file = output + ".iss.tmp.genomes.fasta"
    util.concatenate(genome_files, output=genome_file)

    # for n_genomes we use reservoir sampling to draw random genomes
    # from the concatenated genome file. We then override the file.
    if n_genomes and not draft and not ncbi:
        genome_count = util.count_records(genome_file)
        genome_files = [genome for genome in util.reservoir(SeqIO.parse(genome_file, "fasta"), genome_count, n_genomes)]
        SeqIO.write(genome_files, genome_file, "fasta")

    if os.stat(genome_file).st_size == 0:
        logger.error("Genome(s) file seems empty: %s" % genome_file)
        sys.exit(1)

    try:
        f = open(genome_file, "r")
        with f:  # count the number of records
            genome_list = util.count_records(f)
    except IOError as e:
        logger.error("Failed to open genome(s) file:%s" % e)
        sys.exit(1)

    return genome_list, genome_file


def load_readcount_or_abundance(
    readcount_file,
    abundance_file,
    coverage_file,
    coverage,
    abundance_distribution,
    draft,
    genome_list,
    genome_file,
    n_reads,
    output,
    error_model,
):
    """Load readcount or abundance information based on the provided input parameters.

    Args:
        readcount_file (str): Path to the readcount file.
        abundance_file (str): Path to the abundance file.
        coverage_file (str): Path to the coverage file.
        coverage (str): Coverage distribution type.
        abundance_distribution (str): Abundance distribution type.
        draft (str): Draft mode.
        genome_list (list): List of genomes.
        genome_file (str): Path to the genome file.
        n_reads (int): Number of reads.
        output (str): Output file path.
        error_model (object): Error model object.

    Returns:
        tuple: A tuple containing the readcount dictionary and abundance dictionary.
    """
    logger = logging.getLogger(__name__)
    abundance_dispatch = {
        "uniform": abundance.uniform,
        "halfnormal": abundance.halfnormal,
        "exponential": abundance.exponential,
        "lognormal": abundance.lognormal,
        "zero_inflated_lognormal": abundance.zero_inflated_lognormal,
    }
    readcount_dic = None
    abundance_dic = None
    if readcount_file:
        logger.info("Using readcount file:%s" % readcount_file)
        logger.warning("--readcount_file disables --n_reads, n_reads will be calculated from the readcount file")
        if draft:
            raise RuntimeError("readcount_file is only supported using --genomes, not --draft")
        readcount_dic = abundance.parse_readcount_file(readcount_file)

    elif abundance_file:
        logger.info("Using abundance file:%s" % abundance_file)
        if draft:
            abundance_dic_short = abundance.parse_abundance_file(abundance_file)
            complete_genomes_dic = {k: v for k, v in abundance_dic_short.items() if k not in draft}
            draft_dic = abundance.expand_draft_abundance(abundance_dic_short, draft)
            abundance_dic = {**complete_genomes_dic, **draft_dic}
        else:
            abundance_dic = abundance.parse_abundance_file(abundance_file)
    elif coverage_file:
        logger.warning("--coverage_file is an experimental feature")
        logger.warning("--coverage_file disables --n_reads")
        logger.info("Using coverage file:%s" % coverage_file)
        if draft:
            coverage_dic = abundance.parse_abundance_file(coverage_file)
            complete_genomes_dic = {k: v for k, v in coverage_dic.items() if k not in draft}
            draft_dic = abundance.expand_draft_abundance(coverage_dic, draft, mode="coverage")
            abundance_dic = {**complete_genomes_dic, **draft_dic}
        else:
            abundance_dic = abundance.parse_abundance_file(coverage_file)
    elif coverage in abundance_dispatch:
        logger.info("Using %s coverage distribution" % coverage)
        if draft:
            abundance_dic = abundance.draft(
                genome_list,
                draft,
                abundance_dispatch[abundance_distribution],
                output,
                mode="coverage",
            )
        else:
            abundance_dic = abundance_dispatch[coverage](genome_list)
        if n_reads:
            n_reads = util.convert_n_reads(n_reads)
            logger.info("scaling coverage to %s reads" % n_reads)
            abundance_dic = abundance.coverage_scaling(n_reads, abundance_dic, genome_file, error_model.read_length)
        abundance.to_file(abundance_dic, output, mode="coverage")
    elif abundance_distribution in abundance_dispatch:
        logger.info("Using %s abundance distribution" % abundance_distribution)
        if draft:
            abundance_dic = abundance.draft(
                genome_list,
                draft,
                abundance_dispatch[abundance_distribution],
                output,
            )
        else:
            abundance_dic = abundance_dispatch[abundance_distribution](genome_list)
            abundance.to_file(abundance_dic, output)
    else:
        logger.error("Could not get abundance, or coverage or readcount information")
        sys.exit(1)

    return readcount_dic, abundance_dic


def write_mutations(mutations, mutations_handle):
    """Write mutations to a file

    Args:
        mutations (list): List of mutations.
        mutations_handle (file): File handle to write the mutations to.
    """
    for vcf_dict in mutations:
        mutations_handle.write(
            "\t".join(
                [
                    str(vcf_dict["id"]),
                    str(vcf_dict["position"] + 1),  # vcf files have 1-based index
                    ".",
                    vcf_dict["ref"],
                    str(vcf_dict["alt"]),
                    str(vcf_dict["quality"]),
                    "",
                    "",
                ]
            )
            + "\n"
        )

def create_or_load_gff_db(gff_file):
    """
    Create a GFF database from a GFF file if it doesn't exist, or load it from disk if it does.

    Args:
        gff_file (str): Path to the GFF file.

    Returns:
        gffutils.FeatureDB: The GFF database object.
    """
    dbfn = gff_file + ".db"  # Database file

    # Check if the database already exists
    if os.path.exists(dbfn):
        # Load the existing database
        db = gffutils.FeatureDB(dbfn)
        print(f"Loaded existing GFF database from {dbfn}")
    else:
        # Create a new GFF database and save it
        db = gffutils.create_db(
            gff_file,
            dbfn=dbfn,
            force=True,
            keep_order=True,
            merge_strategy='merge',
            disable_infer_genes=True  # Disable gene inference for GFF files
        )
        print(f"Created and saved new GFF database to {dbfn}")

    return db


def create_or_load_gtf_db(gtf_file):
    """
    Create a GFF database from a GTF file if it doesn't exist, or load it from disk if it does.

    Args:
        gtf_file (str): Path to the GTF file.

    Returns:
        gffutils.FeatureDB: The GFF database object.
    """
    dbfn = gtf_file + ".db"  # Database file

    # Check if the database already exists
    if os.path.exists(dbfn):
        # Load the existing database
        db = gffutils.FeatureDB(dbfn)
        print(f"Loaded existing GFF database from {dbfn}")
    else:
        # Create a new GFF database and save it
        db = gffutils.create_db(gtf_file, dbfn=dbfn, force=True,
                                keep_order=True, merge_strategy='merge',
                                sort_attribute_values='gene_id')
        print(f"Created and saved new GFF database to {dbfn}")

    return db

def create_or_load_gff_db2(gtf_file):
    """
    Create a GFF database from a GTF file if it doesn't exist, or delete it and re-create it.
    """
    dbfn = gtf_file + ".db"  # Database file
    # Remove the DB if it already exists so that a fresh one is created.
    if os.path.exists(dbfn):
        os.remove(dbfn)
        print(f"Removed existing GFF database at {dbfn}")
    db = gffutils.create_db(gtf_file, dbfn=dbfn, force=True,
                            keep_order=True, merge_strategy='merge',
                            sort_attribute_values='gene_id')
    print(f"Created and saved new GFF database to {dbfn}")
    return db



def query_cds_regions(db):
    """
    Query the CDS regions from the GFF database.

    Args:
        db (gffutils.FeatureDB): The GFF database object.

    Returns:
        dict: A dictionary mapping chromosome/plasmid to a list of CDS regions (start, end, strand).
    """
    cds_regions = {}
    for feature in db.features_of_type("CDS"):
        chrom = feature.chrom
        start = feature.start
        end = feature.end
        strand = feature.strand
        if chrom not in cds_regions:
            cds_regions[chrom] = []
        # # add only if on "+"
        # if strand == "+":
        cds_regions[chrom].append((start, end, strand))
    return cds_regions


def update_read_description(read, read_info, cds_regions):
    """
    Update the description of a read with information from the CDS regions.
    read_info = {"start": forward_start, "end": forward_end, "strand": 1,
                 "chrom": header}
    Args:
        read (SeqRecord): The read to update.
        cds_regions (dict): Dictionary containing CDS regions for each sequence.

    Returns:
        SeqRecord: The updated read with the new description.
    """
    is_cds = False
    overlap_percentage = 0
    # Extract the start and end of the read
    read_start = read_info["start"]
    read_end = read_info["end"]
    chrom = read_info["chrom"]
    strand = read_info["strand"]

    # Calculate the total number of bases in the read
    read_length = read_end - read_start

    # Variable to accumulate the total overlap length
    total_overlap = 0
    chrom_ = chrom.split("_")[0] + "_" + chrom.split("_")[1]
    if chrom_ in cds_regions:
        # Loop through all CDS regions to calculate overlap
        for cds_start, cds_end, cds_strand in cds_regions[chrom_]:
            # Check for overlap between read and CDS region
            overlap_start = max(read_start, cds_start)
            overlap_end = min(read_end, cds_end)

            # If there is an overlap (start < end), add it to the total overlap
            if overlap_start < overlap_end:
                total_overlap += (overlap_end - overlap_start)

    # Calculate the percentage of overlap
    if read_length > 0:
        overlap_percentage = (total_overlap / read_length) * 100
        is_cds = overlap_percentage > 0
    else:
        overlap_percentage = 0

    # Build the new description {read_info["strand"]}
    cds_info = 'CDS' if is_cds else 'Non-CDS'
    new_description = f"Start:{read_start}\nEnd:{read_end}\nStrand" \
                      f":{strand}\nChromosome:{chrom_}\n" \
                      f"CDS_info:{cds_info}_" \
                      f"{round(overlap_percentage, 3)}%_coding"

    return new_description


def update_read_description2(read, read_info, cds_regions):
    """
    Update the description of a read with CDS overlap information.
    """
    read_start = read_info["start"]
    read_end = read_info["end"]
    strand = read_info["strand"]
    # Cache parsed chromosome if possible
    chrom_parts = read_info["chrom"].split("_", 2)
    chrom_key = f"{chrom_parts[0]}_{chrom_parts[1]}" if len(
        chrom_parts) >= 2 else read_info["chrom"]

    read_length = read_end - read_start
    total_overlap = 0

    # Assume cds_regions[chrom_key] is pre-sorted by start position
    if chrom_key in cds_regions:
        for cds_start, cds_end, _ in cds_regions[chrom_key]:
            # If the CDS is completely before the read, skip it
            if cds_end < read_start:
                continue
            # Break early if CDS start is beyond the read
            if cds_start > read_end:
                break
            # Calculate overlap
            overlap = min(read_end, cds_end) - max(read_start, cds_start)
            if overlap > 0:
                total_overlap += overlap

    overlap_percentage = (
                                     total_overlap / read_length) * 100 if read_length > 0 else 0
    cds_info = 'CDS' if overlap_percentage > 0 else 'Non-CDS'

    new_description = (
        f"Start:{read_start}\n"
        f"End:{read_end}\n"
        f"Strand:{strand}\n"
        f"Chromosome:{chrom_key}\n"
        f"CDS_info:{cds_info}_{round(overlap_percentage, 3)}%_coding"
    )
    return new_description


def simulate_and_add_info_to_reads(gtf_file, record,error_model,
                                    n_pairs,
                                    cpu_number,
                                    forward_handle,
                                    reverse_handle,
                                    mutations_handle,
                                    sequence_type="metagenomics",
                                    gc_bias=False,
                                    mode="default",
):
    """
    Simulate reads and add detailed information to the description of each read.

    Args:
        gtf_file (str): Path to the GTF file.
        genome_file (str): Path to the genome FASTA file.
        output_dir (str): Directory to save the FASTQ files.
        error_model (ErrorModel): The error model to simulate the reads.
        n_pairs (int): Number of read pairs to generate.
        cpu_number (int): CPU number (used in read ID).
        sequence_type (str): Type of sequencing (default is "metagenomics").
    """
    # Parse the GTF file to get CDS information
    db = create_or_load_gff_db2(gtf_file)
    cds_regions = query_cds_regions(db)

    # Read the genome
    genome_records = SeqIO.parse(record, "fasta")
    n_pairs_per_chrom = round(n_pairs / len(list(genome_records)))

    genome_records = SeqIO.parse(record, "fasta")

    for genome_seqrecord in genome_records:
        # Generate reads for each genome sequence
        generated_reads = reads_generator(n_pairs_per_chrom, genome_seqrecord,
                                       error_model,
                         cpu_number, gc_bias, sequence_type, cds_regions)
        for forward_record, reverse_record, mutations in generated_reads:
            # Update the descriptions of the forward and reverse reads
            # forward_record = update_read_description(forward_record,
            #                                          cds_regions)
            # reverse_record = update_read_description(reverse_record,
            #                                          cds_regions)

            # Write the updated reads to FASTQ files
            SeqIO.write(forward_record, forward_handle, "fastq-sanger")
            SeqIO.write(reverse_record, reverse_handle, "fastq-sanger")
            # Write mutations if needed
            write_mutations(mutations, mutations_handle)
