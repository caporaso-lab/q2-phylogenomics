# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import skbio
import pandas as pd
import os.path
import hashlib

from q2_types.per_sample_sequences import (
    # SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt)
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.feature_data import (DNASequencesDirectoryFormat,
                                   DNAFASTAFormat)

from ._util import run_command
from ._format import (
        BAMFilesDirFmt, PileUpFilesDirFmt,
        BAMFormat, PileUpTSVFormat)


def map_paired_reads(demux: SingleLanePerSamplePairedEndFastqDirFmt,
                     database: Bowtie2IndexDirFmt,
                     mismatches_per_seed: int = 1,
                     ceil_coefficient: float = 0.3,
                     n_threads: int = 1,
                     mapped_only: bool = True) -> BAMFilesDirFmt:
    result = BAMFilesDirFmt()
    ceil_func = '0,%f' % ceil_coefficient
    df = demux.manifest.view(pd.DataFrame)
    for name, fwd, rev in df.itertuples():
        sam = tempfile.NamedTemporaryFile()
        bowtie2_cmd = ['bowtie2', '--threads', str(n_threads), '--phred33',
                       '-N', str(mismatches_per_seed), '--n-ceil', ceil_func,
                       '-x', str(database.path / database.get_basename()),
                       '-1', fwd, '-2', rev, '-S', sam.name]
        run_command(bowtie2_cmd)

        dest = str(result.path / name) + '.bam'
        fixmate_cmd = ['samtools', 'fixmate', '-m']
        if mapped_only:
            fixmate_cmd.append('-r')
        fixmate_cmd.extend([sam.name, dest])
        run_command(fixmate_cmd)

    return result


# def map_single_reads(demux: SingleLanePerSampleSingleEndFastqDirFmt,
#                      database: Bowtie2IndexDirFmt,
#                      mismatches_per_seed: int = 1,
#                      ceil_coefficient: float = 0.3,
#                      n_threads: int = 1,
#                      mapped_only: bool = True) -> BAMFilesDirFmt:
#     result = BAMFilesDirFmt()
#     ceil_func = '0,%f' % ceil_coefficient
#     df = demux.manifest.view(pd.DataFrame)
#     for name, fwd in df.itertuples():
#         sam = tempfile.NamedTemporaryFile()
#         bowtie2_cmd = ['bowtie2', '--threads', str(n_threads), '--phred33',
#                        '-N', str(mismatches_per_seed), '--n-ceil', ceil_func,
#                        '-x', str(database.path / database.get_basename()),
#                        '-U', fwd, '-S', sam.name]
#         run_command(bowtie2_cmd)

#         # need to figure out if this is still the best
#         # way to remove unmapped reads
#         dest = str(result.path / name) + '.bam'
#         fixmate_cmd = ['samtools', 'fixmate', '-m']
#         if mapped_only:
#             fixmate_cmd.append('-r')
#         fixmate_cmd.extend([sam.name, dest])
#         run_command(fixmate_cmd)

#     return result


def sort_alignment_maps(unsorted: BAMFilesDirFmt) -> BAMFilesDirFmt:
    result = BAMFilesDirFmt()
    for path, view in unsorted.bams.iter_views(BAMFormat):
        samtools_sort_cmd = ['samtools', 'sort',
                             '-o', str(result.path / path.stem) + '.bam',
                             str(unsorted.path / path)]
        run_command(samtools_sort_cmd)
    return result


def remove_duplicates(sorted: BAMFilesDirFmt) -> BAMFilesDirFmt:
    result = BAMFilesDirFmt()
    for path, view in sorted.bams.iter_views(BAMFormat):
        samtools_markdup_cmd = ['samtools', 'markdup', '-r',
                                str(sorted.path / path),
                                str(result.path / path)]
        run_command(samtools_markdup_cmd)
    return result


def make_pileups(sorted: BAMFilesDirFmt,
                 reference: DNAFASTAFormat,
                 min_mapq: int = 10,
                 max_depth: int = 300) -> PileUpFilesDirFmt:
    result = PileUpFilesDirFmt()
    for path, view in sorted.bams.iter_views(BAMFormat):
        samtools_mpileup_cmd = ['samtools', 'mpileup',
                                '-f', str(reference),
                                '-q', str(min_mapq),
                                '-d', str(max_depth),
                                '-BA',
                                str(sorted.path / path),
                                '-o', str(result.path / path.stem) + '.tsv']
        run_command(samtools_mpileup_cmd)
    return result


def consensus_sequence(pileups: PileUpFilesDirFmt,
                       min_depth: int = 10,
                       ) -> (pd.DataFrame, DNASequencesDirectoryFormat):
    features = DNASequencesDirectoryFormat()
    # what's the right way to get this file handle?
    sequences_fp = str(features.path / 'dna-sequences.fasta')
    table_data = {}
    sample_ids = set()
    feature_ids = set()
    with open(sequences_fp, 'w') as sequences_f, \
            tempfile.TemporaryDirectory() as tmpdirname:
        for path, view in pileups.pileups.iter_views(PileUpTSVFormat):
            ivar_cmd = ['ivar', 'consensus',
                        '-p', 'consensus-genome',
                        '-n', 'N',
                        '-m', str(min_depth)]
            run_command(ivar_cmd, cwd=tmpdirname,
                        stdin=open(str(pileups.path / path)))

            sample_id = path.stem
            if sample_id in table_data:
                raise ValueError('Sample %s represented by more than one '
                                 'pileup file.' % sample_id)
            sample_ids.add(sample_id)
            table_data[sample_id] = []

            fasta_fp = os.path.join(tmpdirname, 'consensus-genome.fa')
            try:
                feature = skbio.DNA.read(fasta_fp, format='fasta')
            except skbio.io.FASTAFormatError:
                # no genome was assembled from this sample
                continue
            else:
                # is the following consistent with how we hash feature ids
                # in other spots in q2?
                feature_id = hashlib.md5(
                    str(feature).encode('utf-8')).hexdigest()
                table_data[sample_id].append(feature_id)

                if feature_id not in feature_ids:
                    feature_ids.add(feature_id)
                    feature.metadata['id'] = feature_id
                    feature.write(sequences_f, format='fasta')

    # i'm sure there's a better way to do this...
    df = pd.DataFrame(columns=feature_ids, index=sample_ids).fillna(False)
    for sid, fids in table_data.items():
        for fid in fids:
            df[fid][sid] = True
    return df, features
