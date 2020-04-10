# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import (
    Choices,
    Plugin,
    Citations,
    Range,
    Int,
    Str,
    List,
    Bool,
    Float
)
from q2_types.feature_table import FeatureTable, PresenceAbsence
from q2_types.feature_data import FeatureData, Sequence
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    SequencesWithQuality,
    PairedEndSequencesWithQuality,
)

import q2_phylogenomics
import q2_phylogenomics._prinseq
import q2_phylogenomics._filter
import q2_phylogenomics._assemble
from q2_types.bowtie2 import Bowtie2Index
from q2_types.feature_data import DNASequencesDirectoryFormat
from q2_phylogenomics._format import (GenBankFormat, GenBankDirFmt,
                                      BAMFormat, SAMFormat,
                                      BAMFilesDirFmt, SAMFilesDirFmt,
                                      PileUpTSVFormat, PileUpFilesDirFmt,
                                      FASTAFilesDirFmt)
from q2_phylogenomics._types import (AlignmentMap, PileUp, ConsensusSequences,
                                     ReferenceSequence)


citations = Citations.load('citations.bib', package='q2_phylogenomics')

plugin = Plugin(
    name='phylogenomics',
    version=q2_phylogenomics.__version__,
    website='https://github.com/qiime2/q2-phylogenomics',
    package='q2_phylogenomics',
    description='A QIIME 2 plugin for phylogenomics analyses.',
    short_description='A QIIME 2 plugin for phylogenomics analyses.',
)

plugin.register_formats(GenBankFormat, GenBankDirFmt, citations=[])
plugin.register_formats(BAMFormat, SAMFormat, BAMFilesDirFmt, SAMFilesDirFmt,
                        PileUpTSVFormat, PileUpFilesDirFmt,
                        citations=[])
plugin.register_formats(FASTAFilesDirFmt)

plugin.register_semantic_types(AlignmentMap, PileUp, ConsensusSequences,
                               ReferenceSequence)

# before release we want to use GenBank format for this,
# but I think it's broken with skbio < 0.5.6 - I get this
# error when trying to load a genbank file:
# ValueError: cannot set WRITEABLE flag to True of this array
# plugin.register_semantic_type_to_format(ReferenceSequence, GenBankDirFmt)
plugin.register_semantic_type_to_format(ReferenceSequence,
                                        DNASequencesDirectoryFormat)
plugin.register_semantic_type_to_format(SampleData[PileUp], PileUpFilesDirFmt)
plugin.register_semantic_type_to_format(SampleData[AlignmentMap],
                                        BAMFilesDirFmt)
plugin.register_semantic_type_to_format(SampleData[ConsensusSequences],
                                        FASTAFilesDirFmt)

importlib.import_module('q2_phylogenomics._transformers')

prinseq_input = {'demultiplexed_sequences': 'The sequences to be trimmed.'}
prinseq_output = {'trimmed_sequences': 'The resulting trimmed sequences.'}

prinseq_parameters = {
    'trim_qual_right': Int % Range(1, None),
    'trim_qual_type': Str % Choices(['min', 'mean', 'max', 'sum']),
    'trim_qual_window': Int % Range(1, None),
    'min_qual_mean': Int % Range(1, None),
    'min_len': Int % Range(1, None),
    'lc_method': Str % Choices(['dust', 'entropy']),
    'lc_threshold': Int % Range(0, 100),
    'derep': List[Str % Choices(list('12345'))]}

prinseq_parameter_descriptions = {
    'trim_qual_right': 'Trim sequence by quality score from the 3\'-end with '
                       'this threshold score.',
    'trim_qual_type': 'Type of quality score calculation to use. Allowed '
                      'options are min, mean, max and sum.',
    'trim_qual_window': 'The sliding window size used to calculate quality '
                        'score by type. To stop at the first base that fails '
                        'the rule defined, use a window size of 1.',
    'min_qual_mean': 'Filter sequence with quality score mean below '
                     'min_qual_mean.',
    'min_len': 'Filter sequence shorter than min_len.',
    'lc_method': 'Method to filter low complexity sequences.',
    'lc_threshold': 'The threshold value used to filter sequences by sequence '
                    'complexity. The dust method uses this as maximum allowed '
                    'score and the entropy method as minimum allowed value.',
    'derep': 'Type of duplicates to filter. Use integers for multiple '
             'selections (e.g. 124 to use type 1, 2 and 4). The order does '
             'not matter. Option 2 and 3 will set 1 and option 5 will set 4 '
             'as these are subsets of the other option.\n\n1 (exact '
             'duplicate), 2 (5\' duplicate), 3 (3\' duplicate), 4 (reverse '
             'complement exact duplicate), 5 (reverse complement 5\'/3\' '
             'duplicate).'
}

plugin.methods.register_function(
    function=q2_phylogenomics._prinseq.prinseq_single,
    inputs={'demultiplexed_sequences': SampleData[SequencesWithQuality]},
    parameters=prinseq_parameters,
    outputs=[('trimmed_sequences', SampleData[SequencesWithQuality])],
    input_descriptions=prinseq_input,
    parameter_descriptions=prinseq_parameter_descriptions,
    output_descriptions=prinseq_output,
    name='Filter and trim demultiplexed single-end sequences with PRINSEQ.',
    description='Filter and trim demultiplexed single-end FASTQ sequences '
                'based on quality scores using PRINSEQ-lite.',
    citations=[citations['schmieder_prinseq']]
)

plugin.methods.register_function(
    function=q2_phylogenomics._prinseq.prinseq_paired,
    inputs={
        'demultiplexed_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters=prinseq_parameters,
    outputs=[('trimmed_sequences', SampleData[PairedEndSequencesWithQuality])],
    input_descriptions=prinseq_input,
    parameter_descriptions=prinseq_parameter_descriptions,
    output_descriptions=prinseq_output,
    name='Filter and trim demultiplexed paired-end sequences with PRINSEQ.',
    description='Filter and trim demultiplexed paired-end FASTQ sequences '
                'based on quality scores using PRINSEQ-lite.',
    citations=[citations['schmieder_prinseq']]
)

filter_input = {'demultiplexed_sequences': 'The sequences to be trimmed.',
                'database': 'Bowtie2 indexed database.'}
filter_output = {'filtered_sequences': 'The resulting filtered sequences.'}

filter_parameters = {
    'n_threads': Int % Range(1, None),
    'mode': Str % Choices(['local', 'global']),
    'sensitivity': Str % Choices([
        'very-fast', 'fast', 'sensitive', 'very-sensitive']),
    'ref_gap_open_penalty': Int % Range(1, None),
    'ref_gap_ext_penalty': Int % Range(1, None),
    'exclude_seqs': Bool,
}

filter_parameter_descriptions = {
    'n_threads': 'Number of alignment threads to launch.',
    'mode': 'Bowtie2 alignment settings. See bowtie2 manual for more details.',
    'sensitivity': 'Bowtie2 alignment sensitivity. See bowtie2 manual for '
                   'details.',
    'ref_gap_open_penalty': 'Reference gap open penalty.',
    'ref_gap_ext_penalty': 'Reference gap extend penalty.',
    'exclude_seqs': 'Exclude sequences that align to reference. Set this '
                    'option to False to exclude sequences that do not align '
                    'to the reference database.'
}

filter_citations = [citations['langmead2012fast'],
                    citations['heng2009samtools']]

filter_description = (
    'Filter out (or keep) sequences that align to reference database, using '
    'bowtie2 and samtools. This method can be used to filter out human DNA '
    'sequences and other contaminant in any FASTQ sequence data (e.g., '
    'shotgun genome or amplicon sequence data), or alternatively (when '
    'exclude_seqs is False) to only keep sequences that do align to the '
    'reference.')

plugin.methods.register_function(
    function=q2_phylogenomics._filter.filter_single,
    inputs={'demultiplexed_sequences': SampleData[SequencesWithQuality],
            'database': Bowtie2Index},
    parameters=filter_parameters,
    outputs=[('filtered_sequences', SampleData[SequencesWithQuality])],
    input_descriptions=filter_input,
    parameter_descriptions=filter_parameter_descriptions,
    output_descriptions=filter_output,
    name='Filter single-end sequences by alignment to reference database.',
    description=filter_description,
    citations=filter_citations
)

plugin.methods.register_function(
    function=q2_phylogenomics._filter.filter_paired,
    inputs={
        'demultiplexed_sequences': SampleData[PairedEndSequencesWithQuality],
        'database': Bowtie2Index},
    parameters=filter_parameters,
    outputs=[
        ('filtered_sequences', SampleData[PairedEndSequencesWithQuality])],
    input_descriptions=filter_input,
    parameter_descriptions=filter_parameter_descriptions,
    output_descriptions=filter_output,
    name='Filter paired-end sequences by alignment to reference database.',
    description=filter_description,
    citations=filter_citations
)

plugin.methods.register_function(
    function=q2_phylogenomics._filter.bowtie2_build,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={'n_threads': Int % Range(1, None)},
    outputs=[('database', Bowtie2Index)],
    input_descriptions={
        'sequences': 'Reference sequences used to build bowtie2 index.'},
    parameter_descriptions={'n_threads': 'Number of threads to launch'},
    output_descriptions={'database': 'Bowtie2 index.'},
    name='Build bowtie2 index from reference sequences.',
    description='Build bowtie2 index from reference sequences.',
    citations=[citations['langmead2012fast']]
)

map_paired_reads_input_descriptions = {
    'demux': 'The demultiplexed sequences to map to the reference.',
    'database': 'The reference sequence(s).'
}

map_paired_reads_parameter_descriptions = {
    'mismatches_per_seed': 'Max mismatches allowed in seed alignment.',
    'ceil_coefficient': 'Coefficient used to specify the bowtie '
                        'function L(0,x) for max number of non-A/C/G/T '
                        'characters allowed in an alignment.',
    'n_threads': 'Number of alignment threads to launch.',
    'mapped_only': 'Retain only records for reads that were mapped '
                   'to the database in the output files.'
}

map_paired_reads_output_descriptions = {
    'alignment_maps': 'Results of mapping reads in each input sample '
                      'to the provided database.'
}

plugin.methods.register_function(
    function=q2_phylogenomics._assemble.map_paired_reads,
    inputs={'demux': SampleData[PairedEndSequencesWithQuality],
            'database': Bowtie2Index},
    parameters={'mismatches_per_seed': Int % Range(0, 1, inclusive_end=True),
                'ceil_coefficient': Float,
                'n_threads': Int % Range(1, None),
                'mapped_only': Bool},
    outputs=[('alignment_maps', SampleData[AlignmentMap])],
    input_descriptions=map_paired_reads_input_descriptions,
    parameter_descriptions=map_paired_reads_parameter_descriptions,
    output_descriptions=map_paired_reads_output_descriptions,
    name='Map paired end reads.',
    description='Map paired end reads to a database.',
    citations=[citations['langmead2012fast']]
)

sort_alignment_maps_input_descriptions = {
    'unsorted': 'The unsorted alignment maps.'
}

sort_alignment_maps_output_descriptions = {
    'sorted': 'The sorted alignment maps.'
}

plugin.methods.register_function(
    function=q2_phylogenomics._assemble.sort_alignment_maps,
    inputs={'unsorted': SampleData[AlignmentMap]},
    parameters={},
    outputs=[('sorted', SampleData[AlignmentMap])],
    input_descriptions=sort_alignment_maps_input_descriptions,
    parameter_descriptions={},
    output_descriptions=sort_alignment_maps_output_descriptions,
    name='Sort alignment maps.',
    description='Sort alignment maps by reference start position.',
    citations=[citations['heng2009samtools']]
)

remove_duplicates_input_descriptions = {
    'sorted': 'The sorted alignment maps.'
}

remove_duplicates_output_descriptions = {
    'duplicate_filtered': 'The sorted and filtered alignment maps.'
}

plugin.methods.register_function(
    function=q2_phylogenomics._assemble.remove_duplicates,
    inputs={'sorted': SampleData[AlignmentMap]},
    parameters={},
    outputs=[('duplicate_filtered', SampleData[AlignmentMap])],
    input_descriptions=remove_duplicates_input_descriptions,
    parameter_descriptions={},
    output_descriptions=remove_duplicates_output_descriptions,
    name='Remove duplicates.',
    description='Remove duplicate reads from alignment maps.',
    citations=[citations['heng2009samtools']]
)

make_pileup_input_descriptions = {
    'sorted': 'Sorted alignment maps.',
    'reference': 'The reference sequence'
}

make_pileup_parameter_descriptions = {
    'min_mapq': 'The minimum mapQ to consider an alignment.',
    'max_depth': 'The max per-file depth.'
}

make_pileup_output_descriptions = {
    'pileups': 'The resulting PileUp data.'
}

plugin.methods.register_function(
    function=q2_phylogenomics._assemble.make_pileups,
    inputs={'sorted': SampleData[AlignmentMap],  # need a sorted property?
            # the following should become type ReferenceSequence
            # or somehow be integrated with the Bowtie Index
            'reference': Bowtie2Index},
    parameters={'min_mapq': Int % Range(0, None),
                'max_depth': Int % Range(1, None)},
    outputs=[('pileups', SampleData[PileUp])],
    input_descriptions=make_pileup_input_descriptions,
    parameter_descriptions=make_pileup_parameter_descriptions,
    output_descriptions=make_pileup_output_descriptions,
    name='Create PileUp files',
    description='Create PileUp Files from sorted alignment maps.',
    citations=[citations['heng2009samtools']]
)

consensus_sequence_input_descriptions = {
    'pileups': 'The PileUp data.'
}

consensus_sequence_output_descriptions = {
    'table': 'Table describing which consensus sequences are '
             'observed in each sample.',
    'consensus_sequences': 'Mapping of consensus sequence identifiers '
                           'to consensus sequences.'
}

plugin.methods.register_function(
    function=q2_phylogenomics._assemble.consensus_sequence,
    inputs={'pileups': SampleData[PileUp]},
    parameters={},
    outputs=[('table', FeatureTable[PresenceAbsence]),
             ('consensus_sequences', FeatureData[Sequence])],
    input_descriptions=consensus_sequence_input_descriptions,
    parameter_descriptions={},
    output_descriptions=consensus_sequence_output_descriptions,
    name='',
    description='',
    citations=[citations['Grubaugh2019ivar']]
)
