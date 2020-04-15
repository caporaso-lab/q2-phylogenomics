# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._prinseq import _prinseq_defaults
from ._filter import _filter_defaults


_cut_defaults = {
    'cores': 1,
    'adapter_f': ['GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'],
    'adapter_r': ['AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
    'front_f': None,
    'front_r': None,
    'anywhere_f': None,
    'anywhere_r': None,
    'error_rate': 0.1,
    'indels': True,
    'times': 3,
    'overlap': 3,
    'match_read_wildcards': False,
    'match_adapter_wildcards': True,
    'minimum_length': 1,
    'discard_untrimmed': False,
}

_map_defaults = {
    'mismatches_per_seed': 1,
    'ceil_coefficient': 0.3,
    'n_threads': 1,
    'mapped_only': True,
}

_pileup_defaults = {
    'min_mapq': 10,
    'max_depth': 300,
}

_consensus_defaults = {
    'min_depth': 10
}


def filter_clean_consensus(
        ctx,
        demultiplexed_sequences,
        alignment_ref,
        filter_ref=None,
        enable_cutadapt=True,
        enable_prinseq=True,
        cutadapt_cores=_cut_defaults['cores'],
        cutadapt_adapter_f=_cut_defaults['adapter_f'],
        cutadapt_front_f=_cut_defaults['front_f'],
        cutadapt_anywhere_f=_cut_defaults['anywhere_f'],
        cutadapt_adapter_r=_cut_defaults['adapter_r'],
        cutadapt_front_r=_cut_defaults['front_r'],
        cutadapt_anywhere_r=_cut_defaults['anywhere_r'],
        cutadapt_error_rate=_cut_defaults['error_rate'],
        cutadapt_indels=_cut_defaults['indels'],
        cutadapt_times=_cut_defaults['times'],
        cutadapt_overlap=_cut_defaults['overlap'],
        cutadapt_match_read_wildcards=_cut_defaults['match_read_wildcards'],
        cutadapt_match_adapter_wildcards=_cut_defaults[
            'match_adapter_wildcards'],
        cutadapt_minimum_length=_cut_defaults['minimum_length'],
        cutadapt_discard_untrimmed=_cut_defaults['discard_untrimmed'],
        bowtie2_n_threads=_filter_defaults['n_threads'],
        bowtie2_mode=_filter_defaults['mode'],
        bowtie2_sensitivity=_filter_defaults['sensitivity'],
        bowtie2_ref_gap_open_penalty=_filter_defaults['ref_gap_open_penalty'],
        bowtie2_ref_gap_ext_penalty=_filter_defaults['ref_gap_ext_penalty'],
        bowtie2_exclude_seqs=_filter_defaults['exclude_seqs'],
        bowtie2_mismatches_per_seed=_map_defaults['mismatches_per_seed'],
        bowtie2_ceil_coefficient=_map_defaults['ceil_coefficient'],
        bowtie2_mapped_only=_map_defaults['mapped_only'],
        prinseq_trim_qual_right=_prinseq_defaults['trim_qual_right'],
        prinseq_trim_qual_type=_prinseq_defaults['trim_qual_type'],
        prinseq_trim_qual_window=_prinseq_defaults['trim_qual_window'],
        prinseq_min_qual_mean=_prinseq_defaults['min_qual_mean'],
        prinseq_min_len=_prinseq_defaults['min_len'],
        prinseq_lc_method=_prinseq_defaults['lc_method'],
        prinseq_lc_threshold=_prinseq_defaults['lc_threshold'],
        prinseq_derep=_prinseq_defaults['derep'],
        samtools_min_mapq=_pileup_defaults['min_mapq'],
        samtools_max_depth=_pileup_defaults['max_depth'],
        consensus_min_depth=_consensus_defaults['min_depth']
        ):

    # Quality Control
    if filter_ref is not None:
        filter = ctx.get_action('phylogenomics', 'filter_paired')
        filtered_seqs, = filter(
            demultiplexed_sequences=demultiplexed_sequences,
            database=filter_ref,
            n_threads=bowtie2_n_threads,
            mode=bowtie2_mode,
            sensitivity=bowtie2_sensitivity,
            ref_gap_open_penalty=bowtie2_ref_gap_open_penalty,
            ref_gap_ext_penalty=bowtie2_ref_gap_ext_penalty,
            exclude_seqs=bowtie2_exclude_seqs)
    else:
        filtered_seqs = demultiplexed_sequences

    if enable_cutadapt:
        cut = ctx.get_action('cutadapt', 'trim_paired')
        cut_seqs, = cut(
            demultiplexed_sequences=filtered_seqs,
            cores=cutadapt_cores,
            adapter_f=cutadapt_adapter_f,
            front_f=cutadapt_front_f,
            anywhere_f=cutadapt_anywhere_f,
            adapter_r=cutadapt_adapter_r,
            front_r=cutadapt_front_r,
            anywhere_r=cutadapt_anywhere_r,
            error_rate=cutadapt_error_rate,
            indels=cutadapt_indels,
            times=cutadapt_times,
            overlap=cutadapt_overlap,
            match_read_wildcards=cutadapt_match_read_wildcards,
            match_adapter_wildcards=cutadapt_match_adapter_wildcards,
            minimum_length=cutadapt_minimum_length,
            discard_untrimmed=cutadapt_discard_untrimmed)
    else:
        cut_seqs = filtered_seqs

    if enable_prinseq:
        clean = ctx.get_action('phylogenomics', 'prinseq_paired')
        clean_seqs, = clean(
            demultiplexed_sequences=cut_seqs,
            trim_qual_right=prinseq_trim_qual_right,
            trim_qual_type=prinseq_trim_qual_type,
            trim_qual_window=prinseq_trim_qual_window,
            min_qual_mean=prinseq_min_qual_mean,
            min_len=prinseq_min_len,
            lc_method=prinseq_lc_method,
            lc_threshold=prinseq_lc_threshold,
            derep=prinseq_derep)
    else:
        clean_seqs = cut_seqs

    # Map reads and build consensus
    map_reads = ctx.get_action('phylogenomics', 'map_paired_reads')
    sort_aln = ctx.get_action('phylogenomics', 'sort_alignment_maps')
    dedup = ctx.get_action('phylogenomics', 'remove_duplicates')
    make_pileup = ctx.get_action('phylogenomics', 'make_pileups')
    make_consensus = ctx.get_action('phylogenomics', 'consensus_sequence')
    mapped_reads, = map_reads(
        demux=clean_seqs,
        database=alignment_ref,
        mismatches_per_seed=bowtie2_mismatches_per_seed,
        ceil_coefficient=bowtie2_ceil_coefficient,
        n_threads=bowtie2_n_threads,
        mapped_only=bowtie2_mapped_only)
    mapped_reads, = sort_aln(mapped_reads)
    mapped_reads, = dedup(mapped_reads)
    pileup, = make_pileup(
        sorted=mapped_reads,
        reference=alignment_ref,
        min_mapq=samtools_min_mapq,
        max_depth=samtools_min_mapq)
    table, features, = make_consensus(pileup)

    return table, features, filtered_seqs, clean_seqs
