# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import hashlib

import qiime2
from qiime2.plugin.testing import TestPluginBase
import pandas as pd

from q2_phylogenomics._format import (
    SAMFormat, SAMFilesDirFmt, PileUpFilesDirFmt, PileUpTSVFormat)


class TestMapPaired(TestPluginBase):
    package = 'q2_phylogenomics.tests'

    def setUp(self):
        super().setUp()

        self.demuxed_art = qiime2.Artifact.load(
            self.get_data_path('paired-end.qza'))
        self.paired_mapped_unsorted = qiime2.Artifact.load(
            self.get_data_path('paired-end-mapped-unsorted.qza'))
        self.indexed_genome = qiime2.Artifact.load(
            self.get_data_path('sars2-indexed.qza'))
        self.sorted_alignment_maps = qiime2.Artifact.load(
            self.get_data_path('sorted-alignment-maps.qza'))
        self.pileups = qiime2.Artifact.load(
            self.get_data_path('pileups.qza')
        )

    def test_map_paired_mapped_only(self):
        obs_art, = self.plugin.methods['map_paired_reads'](
            self.demuxed_art, self.indexed_genome)
        obs = obs_art.view(SAMFilesDirFmt)
        exp = [('sample_a.sam', 10,
                ('SARS2:6:73:567:7631#0', 'SARS2:6:73:233:3421#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0')),
               ('sample_b.sam', 10,
                ('SARS2:6:73:941:1973#0', 'SARS2:6:73:552:2457#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0')),
               ('sample_c.sam', 10,
                ('SARS2:6:73:231:3321#0', 'SARS2:6:73:552:2457#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0'))]
        obs_sams = obs.sams.iter_views(SAMFormat)
        for (obs_fn, obs_sam), (exp_fn, num_records, mapped_ids, unmapped_ids)\
                in zip(obs_sams, exp):
            self.assertEqual(str(obs_fn), exp_fn)
            with open(str(obs_sam)) as sam_f:
                obs_mapped_ids = [line.split('\t')[0]
                                  for line in sam_f
                                  if not line.startswith('@')]
            self.assertEqual(len(obs_mapped_ids), num_records)
            for e in mapped_ids:
                self.assertIn(e, obs_mapped_ids)
            for e in unmapped_ids:
                self.assertNotIn(e, obs_mapped_ids)

    def test_map_paired_not_mapped_only(self):
        obs_art, = self.plugin.methods['map_paired_reads'](
            self.demuxed_art, self.indexed_genome, mapped_only=False)
        obs = obs_art.view(SAMFilesDirFmt)
        exp = [('sample_a.sam', 12,
                ('SARS2:6:73:567:7631#0', 'SARS2:6:73:233:3421#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0')),
               ('sample_b.sam', 12,
                ('SARS2:6:73:941:1973#0', 'SARS2:6:73:552:2457#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0')),
               ('sample_c.sam', 12,
                ('SARS2:6:73:231:3321#0', 'SARS2:6:73:552:2457#0'),
                ('SARS2:6:73:356:9806#0', 'SARS2:6:73:356:9806#0'))]
        obs_sams = obs.sams.iter_views(SAMFormat)
        for (obs_fn, obs_sam), (exp_fn, num_records, mapped_ids, unmapped_ids)\
                in zip(obs_sams, exp):
            self.assertEqual(str(obs_fn), exp_fn)
            with open(str(obs_sam)) as sam_f:
                obs_mapped_ids = [line.split('\t')[0]
                                  for line in sam_f
                                  if not line.startswith('@')]
            self.assertEqual(len(obs_mapped_ids), num_records)
            for e in mapped_ids:
                self.assertIn(e, obs_mapped_ids)
            for e in unmapped_ids:
                self.assertIn(e, obs_mapped_ids)

    def test_map_paired_alt_ceil_coefficient(self):
        obs_art, = self.plugin.methods['map_paired_reads'](
            self.demuxed_art, self.indexed_genome, ceil_coefficient=0.5)
        obs = obs_art.view(SAMFilesDirFmt)
        obs_sams = obs.sams.iter_views(SAMFormat)
        for _, obs_sam in obs_sams:
            with open(str(obs_sam)) as sam_f:
                self.assertIn('--n-ceil 0,0.5', sam_f.read())

    def test_sort_alignment_maps(self):
        obs_art, = self.plugin.methods['sort_alignment_maps'](
            self.paired_mapped_unsorted)
        obs = obs_art.view(SAMFilesDirFmt)
        exp_mapped_positions = \
            [1, 1, 192, 211, 402, 421, 612, 631, 823, 841]
        obs_sams = obs.sams.iter_views(SAMFormat)
        for _, obs_sam in obs_sams:
            with open(str(obs_sam)) as sam_f:
                obs_mapped_positions = \
                    [int(line.split('\t')[3])
                     for line in sam_f
                     if not line.startswith('@')]
                self.assertEqual(obs_mapped_positions,
                                 exp_mapped_positions)

    def test_remove_duplicates(self):
        obs_art, = self.plugin.methods['remove_duplicates'](
            self.sorted_alignment_maps)
        obs = obs_art.view(SAMFilesDirFmt)
        obs_sams = obs.sams.iter_views(SAMFormat)
        for _, obs_sam in obs_sams:
            with open(str(obs_sam)) as sam_f:
                obs_mapped_ids = [line.split('\t')[0]
                                  for line in sam_f
                                  if not line.startswith('@')]
            # one occurance of duplicate alignment is retained...
            self.assertIn('NB501727:157:HFWHJBGXF:2:22105:18312:6802',
                          obs_mapped_ids)
            # and the other isn't
            self.assertNotIn('NB501727:157:HFWHJBGXF:3:23610:2922:9731',
                             obs_mapped_ids)

    def test_make_pileups(self):
        obs_art, = self.plugin.methods['make_pileups'](
            self.sorted_alignment_maps, self.indexed_genome
        )
        obs = obs_art.view(PileUpFilesDirFmt)
        obs_pileups = obs.pileups.iter_views(PileUpTSVFormat)
        for _, obs_pileup in obs_pileups:
            obs_df = pd.read_csv(str(obs_pileup), header=None, sep='\t', )
            # expected values are derived from running samtools
            # directly on these input files
            self.assertEqual(obs_df.shape, (345, 6))
            self.assertEqual(list(obs_df.iloc[:, 0]), ['NC_045512.2'] * 345)

    def test_consensus_sequence_min_depth_1(self):
        # this min depth value gets actual sequence with the test data
        obs_table_art, obs_seq_art = self.plugin.methods['consensus_sequence'](
            self.pileups, min_depth=1,
        )
        obs_table = obs_table_art.view(pd.DataFrame)

        # table tests
        # two different genomes across four samples
        self.assertEqual(obs_table.shape, (4, 2))

        # expected sample ids
        self.assertEqual(
                set(['sample_1', 'sample_2', 'sample_3', 'empty']),
                set(obs_table.index))

        # expected feature ids
        seq1_md5 = 'e8b3172acb8547d54deb27e85b596233'  # in samples 1 & 2
        seq2_md5 = '940f1f1bb24a34601d20afbac3147543'  # in sample 3
        self.assertEqual(set([seq1_md5, seq2_md5]), set(obs_table.columns))

        # expected total counts
        self.assertEqual(obs_table.loc['sample_1'].sum(), 1)
        self.assertEqual(obs_table.loc['sample_2'].sum(), 1)
        self.assertEqual(obs_table.loc['sample_3'].sum(), 1)
        self.assertEqual(obs_table.loc['empty'].sum(), 0)
        self.assertEqual(
            obs_table.loc[:, seq1_md5].sum(), 2)
        self.assertEqual(
            obs_table.loc[:, seq2_md5].sum(), 1)

        self.assertTrue(obs_table[seq1_md5]['sample_1'])
        self.assertFalse(obs_table[seq2_md5]['sample_1'])

        self.assertTrue(obs_table[seq1_md5]['sample_2'])
        self.assertFalse(obs_table[seq2_md5]['sample_2'])

        self.assertFalse(obs_table[seq1_md5]['sample_3'])
        self.assertTrue(obs_table[seq2_md5]['sample_3'])

        self.assertFalse(obs_table[seq1_md5]['empty'])
        self.assertFalse(obs_table[seq2_md5]['empty'])

        # sequence collection tests
        self.obs_seq = obs_seq_art.view(pd.Series)
        self.assertEqual(len(self.obs_seq), 2)
        # confirm expected sequences by computing their md5 here
        # and comparing to the expected sequence ids (which were
        # determined independently)
        self.assertEqual(
            hashlib.md5(str(self.obs_seq[seq1_md5]).
                        encode('utf-8')).hexdigest(),
            seq1_md5)
        self.assertEqual(
            hashlib.md5(str(self.obs_seq[seq2_md5]).
                        encode('utf-8')).hexdigest(),
            seq2_md5)

    def test_consensus_sequence_min_depth_default(self):
        # this min depth value results in sequences that are
        # 100% N with the test data
        obs_table_art, obs_seq_art = self.plugin.methods['consensus_sequence'](
            self.pileups,
        )
        obs_table = obs_table_art.view(pd.DataFrame)

        # table tests
        # one genome across four samples
        self.assertEqual(obs_table.shape, (4, 1))

        # expected sample ids
        self.assertEqual(
                set(['sample_1', 'sample_2', 'sample_3', 'empty']),
                set(obs_table.index))

        # expected feature ids
        seq1_md5 = '53b4ce3236e629aebc69f8a2b5abc96b'  # in samples 1-3 (all N)
        self.assertEqual([seq1_md5], obs_table.columns)

        # expected total counts
        self.assertEqual(obs_table.loc['sample_1'].sum(), 1)
        self.assertEqual(obs_table.loc['sample_2'].sum(), 1)
        self.assertEqual(obs_table.loc['sample_3'].sum(), 1)
        self.assertEqual(obs_table.loc['empty'].sum(), 0)
        self.assertEqual(
            obs_table.loc[:, seq1_md5].sum(), 3)

        self.assertTrue(obs_table[seq1_md5]['sample_1'])
        self.assertTrue(obs_table[seq1_md5]['sample_2'])
        self.assertTrue(obs_table[seq1_md5]['sample_3'])
        self.assertFalse(obs_table[seq1_md5]['empty'])

        # sequence collection tests
        self.obs_seq = obs_seq_art.view(pd.Series)
        self.assertEqual(len(self.obs_seq), 1)
        # confirm expected sequences by computing their md5 here
        # and comparing to the expected sequence ids (which were
        # determined independently)
        self.assertEqual(
            hashlib.md5(str(self.obs_seq[seq1_md5]).
                        encode('utf-8')).hexdigest(),
            seq1_md5)


if __name__ == '__main__':
    unittest.main()
