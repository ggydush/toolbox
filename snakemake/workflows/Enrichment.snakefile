from collections import namedtuple
import os
import re
import shutil

import pandas as pd
import numpy as np


# Software constants
JAVA_PARAMS = "-Xmx14g -Djava.io.tmpdir=./tmp/"
PICARD = f"picard {JAVA_PARAMS}"
FGBIO = f"fgbio {JAVA_PARAMS}"
GATK = f"gatk {JAVA_PARAMS}"
BWA = "bwa"
SUMMARIZE_MUTANTS = "python /xchip/bloodbiopsy/Greg/scripts/get_mutant_metrics.py"

# Reference constants
REF = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
DBSNP = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.dbsnp.vcf"
KNOWN_INDELS = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.known_indels.vcf"
VARIANT_GOLD_STANDARD = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.variantEvalGoldStandard.vcf"

# Snakemake Setup
def parse_metadata_file(metadata_file):
    """
    Expects metadata to have the following columns:

    bam_file sample_id library bait_intervals target_intervals mut_maf
    Optional columns: het_vcf, normal_id
    """
    return pd.read_csv(metadata_file, sep="\t").set_index("sample_id")


def get_info_from_metadata(column_name):
    """ Extract info from metadata """
    return lambda wildcards: metadata.loc[wildcards.sample][column_name]


META_FILE = "sample_metadata.txt"
metadata = parse_metadata_file(META_FILE)
samples = metadata.index
workdir: "../"

# Rules
rule all:
    input:
        "summary_metrics.txt"


rule CollectMultipleMetrics:
    input:
        get_info_from_metadata("bam_file")
    output:
        "CollectMultipleMetrics/{sample}.alignment_summary_metrics",
        "CollectMultipleMetrics/{sample}.insert_size_metrics",
        "CollectMultipleMetrics/{sample}.quality_distribution_metrics"
    resources:
        mem = 12,
        runtime = 16
    params:
        prefix = "CollectMultipleMetrics/{sample}",
        reference = REF
    shell:
        """
        {PICARD} CollectMultipleMetrics \
            I={input} \
            O={params.prefix} \
            R={params.reference} \
            PROGRAM=CollectAlignmentSummaryMetrics \
            PROGRAM=CollectInsertSizeMetrics
        """


rule CollectHsMetrics:
    input:
        get_info_from_metadata("bam_file")
    output:
        hs_metrics = "CollectHsMetrics/{sample}.hs_metrics.txt",
        per_target = "CollectHsMetrics/{sample}.per_target.txt"
    resources:
        mem = 12,
        runtime = 16
    params:
        reference = REF,
        bait = get_info_from_metadata("bait_intervals"),
        target = get_info_from_metadata("target_intervals")
    shell:
        """
        {PICARD} CollectHsMetrics \
            I={input} \
            O={output.hs_metrics} \
            R={params.reference} \
            BI={params.bait} \
            TI={params.target} \
            PER_TARGET_COVERAGE={output.per_target}
        """


rule AddReadGroups:
    input:
        get_info_from_metadata("bam_file")
    output:
        temp("{sample}_{library}_addedRG.bam")
    params:
        rgsm = get_info_from_metadata("library"),
        rglb = get_info_from_metadata("library"),
        rgid = lambda wildcards: wildcards.sample
    resources:
        mem = 8,
        runtime = 24
    shell:
        """
        {PICARD} AddOrReplaceReadGroups \
            INPUT={input} \
            OUTPUT={output} \
            RGID={params.rgid} \
            RGLB={params.rglb} \
            RGPL=illumina \
            RGPU=NA \
            RGSM={params.rgsm} \
        """


def create_library_sample_mapping(metadata):
    sample_mapping = metadata.reset_index().groupby('library').sample_id.agg(list).to_dict()
    library_mapping = metadata['library'].to_dict()
    return sample_mapping, library_mapping

sample_mapping, library_mapping = create_library_sample_mapping(metadata)


rule MergeBams:
    input:
        lambda wildcards: expand("{sample}_{{library}}_addedRG.bam", sample=sample_mapping[wildcards.library])
    output:
        temp("{library}_merged.bam")
    run:
        inputs = " ".join(["INPUT={}".format(sample) for sample in input])
        shell("{PICARD} MergeSamFiles {inputs} OUTPUT={output}")


rule GroupReadsByUmi:
    input:
        "{library}_merged.bam"
    output:
        bam = temp("{library}.GroupedByUmi.bam"),
        histogram = "GroupReadsByUmi/{library}_umiHistogram.txt"
    resources:
        mem = 16,
        runtime = 32
    shell:
        """
        {FGBIO} GroupReadsByUmi \
            -i {input} \
            -o {output.bam} \
            -f {output.histogram} \
            --strategy=Paired \
        """


rule SplitGroupReadsByUmi:
    input:
        lambda wildcards: expand("{library}.GroupedByUmi.bam", library=library_mapping[wildcards.sample])
    output:
        "GroupReadsByUmi/{sample}.GroupByUmi.bam"
    params:
        rgid = lambda wildcards: wildcards.sample
    resources:
        mem = 8,
        runtime = 12
    shell:
        "samtools view -bhr {params.rgid} {input} > {output}"


rule CollectDuplexMetrics:
    input:
        "GroupReadsByUmi/{sample}.GroupByUmi.bam"
    output:
        duplex_family_sizes = "DuplexMetrics/{sample}.duplex_family_sizes.txt",
        duplex_yield_metrics = "DuplexMetrics/{sample}.duplex_yield_metrics.txt",
        family_sizes = "DuplexMetrics/{sample}.family_sizes.txt",
        umi_counts = "DuplexMetrics/{sample}.umi_counts.txt"
    resources:
        mem = 12,
        runtime = 12
    params:
        interval = get_info_from_metadata("bait_intervals"),
        out_prefix = "DuplexMetrics/{sample}"
    shell:
        """
        {FGBIO} CollectDuplexSeqMetrics \
            -i {input} \
            -l {params.interval} \
            -a 2 \
            -b 2 \
            -o {params.out_prefix}
        """


rule CallDuplexConsensusReads:
    input:
        "GroupReadsByUmi/{sample}.GroupByUmi.bam"
    output:
        temp("{sample}.DSC.consensus.bam")
    resources:
        mem = 12,
        runtime = 12
    shell:
        """
        {FGBIO} CallDuplexConsensusReads \
            -i {input} \
            -o {output} \
            --trim true
        """

rule AlignConsensusReads:
    input:
        "{sample}.DSC.consensus.bam"
    output:
        fastq1 = temp("{sample}_DSC_1.fq"),
        fastq2 = temp("{sample}_DSC_2.fq"),
        tmp = temp("{sample}.DSC.aligned_tmp.bam"),
        bam = temp("{sample}.DSC.consensus_realigned.bam")
    params:
        reference = get_info_from_metadata("reference")
    resources:
        mem = 12,
        runtime = 16
    shell:
        """
        {PICARD} SamToFastq \
            I={input} \
            FASTQ={output.fastq1} \
            SECOND_END_FASTQ={output.fastq2} &&

        {BWA} mem \
            -K 100000000 \
            -t 16 \
            {params.reference} {output.fastq1} {output.fastq2} > {output.tmp} &&

        {PICARD} MergeBamAlignment \
            ALIGNED={output.tmp} \
            UNMAPPED={input} \
            O={output.bam} \
            R={params.reference}
        """

rule IndelRealignmentDSC:
    input:
        bam = "{sample}.DSC.consensus_realigned.bam",
        targets = get_info_from_metadata("bait_intervals")
    output:
        bam = temp("{sample}.DSC.consensus_indel_realigned.bam"),
        index = temp("{sample}.DSC.consensus_indel_realigned.bai")
    params:
        reference = get_info_from_metadata("reference")
    resources:
        mem = 16,
        runtime = 4
    shell:
        """
        samtools index {input.bam} &&

        {GATK} -T IndelRealigner \
            -I {input.bam} \
            -o {output.bam} \
            -R {params.reference} \
            -allowPotentiallyMisencodedQuals \
            -targetIntervals {input.targets} \
            -model USE_READS \
            -known {DBSNP} \
            -known {KNOWN_INDELS} \
            -known {VARIANT_GOLD_STANDARD}
        """


rule FilterConsensusReads:
    input:
        "{sample}.DSC.consensus_indel_realigned.bam"
    output:
        bam = "FilterConsensusReads/{sample}.dsc_consensus.bam",
        index = "FilterConsensusReads/{sample}.dsc_consensus.bai"
    params:
        reference = get_info_from_metadata("reference")
    shell:
        """
        {FGBIO} FilterConsensusReads \
            -i {input} \
            -r {params.reference} \
            -o {output.bam} \
            -M 2 \
            -N 20 \
            --max-read-error-rate 0.01 \
            --max-base-error-rate 0.05 \
            --max-no-call-fraction 0.05 \
            --min-mean-base-quality 40 \
            --require-single-strand-agreement true
        """


rule CallMolecularConsensusReads:
    input:
        "GroupReadsByUmi/{sample}.GroupByUmi.bam"
    output:
        temp("{sample}.SSC.consensus.bam")
    resources:
        mem = 12,
        runtime = 12
    shell:
        """
        {FGBIO} CallMolecularConsensusReads \
            -i {input} \
            -o {output} \
            -M 1
        """


rule AlignMolecularConsensusReads:
    input:
        "{sample}.SSC.consensus.bam"
    output:
        fastq1 = temp("{sample}_SSC_1.fq"),
        fastq2 = temp("{sample}_SSC_2.fq"),
        tmp = temp("{sample}.SSC.aligned_tmp.bam"),
        bam = "{sample}.SSC.consensus_realigned.bam"
    params:
        reference = REF
    resources:
        mem = 16,
        runtime = 8
    shell:
        """
        {PICARD} SamToFastq \
            I={input} \
            FASTQ={output.fastq1} \
            SECOND_END_FASTQ={output.fastq2} &&

        {BWA} mem \
            -K 100000000 \
            -t 16 \
            {params.reference} {output.fastq1} {output.fastq2} > {output.tmp} &&

        {PICARD} MergeBamAlignment \
            ALIGNED={output.tmp} \
            UNMAPPED={input} \
            O={output.bam} \
            R={params.reference}
        """


rule IndelRealignmentSSC:
    input:
        bam = "{sample}.SSC.consensus_realigned.bam",
        targets = get_info_from_metadata("bait_intervals")
    output:
        bam = temp("{sample}.SSC.consensus_indel_realigned.bam"),
        index = temp("{sample}.SSC.consensus_indel_realigned.bai")
    params:
        reference = get_info_from_metadata("reference")
    resources:
        mem = 16,
        runtime = 8
    shell:
        """
        samtools index {input.bam} &&

        {GATK} -T IndelRealigner \
            -I {input.bam} \
            -o {output.bam} \
            -R {params.reference} \
            -allowPotentiallyMisencodedQuals \
            -targetIntervals {input.targets} \
            -model USE_READS \
            -known {DBSNP} \
            -known {KNOWN_INDELS} \
            -known {VARIANT_GOLD_STANDARD}
        """


rule FilterMolecularConsensusReads:
    input:
        "{sample}.SSC.consensus_indel_realigned.bam"
    output:
        bam = "FilterMolecularConsensusReads/{sample}.ssc_consensus.bam",
        index = "FilterMolecularConsensusReads/{sample}.ssc_consensus.bai"
    params:
        reference = get_info_from_metadata("reference")
    shell:
        """
        {FGBIO} FilterConsensusReads \
            -i {input} \
            -r {params.reference} \
            -o {output.bam} \
            -M 2 \
            -N 20 \
            --max-read-error-rate 0.01 \
            --max-base-error-rate 0.05 \
            --max-no-call-fraction 0.05 \
            --min-mean-base-quality 40 \
            --require-single-strand-agreement true
        """


def get_normal_and_het(metadata):
    def value_if_not_null(row, column):
        if column in row:
            if row[column] == row[column]:
                return row[column]
        return False
    res = {}
    for sample_id, row in metadata.iterrows():
        normal = value_if_not_null(row, 'normal_id') or sample_id
        het = value_if_not_null(row, 'het_vcf') or None
        res[sample_id] = {'normal_id': normal, 'het_vcf': het}
    return res

germline_info = get_normal_and_het(metadata)

rule MiredasDetectFingerprint:
    input:
        bam = "FilterConsensusReads/{sample}.dsc_consensus.bam",
        ssc_bam = "FilterMolecularConsensusReads/{sample}.ssc_consensus.bam",
        normal = lambda wildcards: f"FilterConsensusReads/{germline_info[wildcards.sample]['normal_id']}.dsc_consensus.bam"
    output:
        detect = "Miredas/{sample}.miredas_detect_fingerprint.txt",
        families = "Miredas/{sample}.miredas_detect_families.txt"
    resources:
        mem = 4,
        runtime = 2
    params:
        mut_maf = get_info_from_metadata("mut_maf")
    run:
        cmd = """
                miredas DetectFingerprint \
                    -b {input.bam} \
                    -s {input.ssc_bam} \
                    -m {params.mut_maf} \
                    --min_max_read_position 13 \
                    --max_alt_bases 2 \
                    -o {output.detect} \
                    --mutant_families_outfile {output.families} \
                """
        if input.bam != input.normal:
            cmd += "--normal_bam {input.normal}"
        shell(cmd)


rule MiredasCollectErrorMetrics:
    input:
        bam = "FilterConsensusReads/{sample}.dsc_consensus.bam",
        het_vcf = lambda wildcards: germline_info[wildcards.sample]['het_vcf'] or "FilterConsensusReads/{sample}.dsc_consensus.bam",
        normal = lambda wildcards: f"FilterConsensusReads/{germline_info[wildcards.sample]['normal_id']}.dsc_consensus.bam"
    output:
        "Miredas/{sample}._error_metrics.txt"
    params:
        sample_id = "{sample}",
        reference = REF,
        mut_maf = get_info_from_metadata("mut_maf"),
        interval_list = get_info_from_metadata("bait_intervals")
    resources:
        mem = 8,
        runtime = 96
    run:
        cmd = """
                miredas CollectErrorMetrics \
                    -b {input.bam} \
                    -s {params.sample_id} \
                    -m {params.mut_maf} \
                    -i {params.interval_list} \
                    -r {params.reference} \
                    --min_max_read_position 12 \
                    --max_alt_bases 2 \
                    -o Miredas/{params.sample_id}. \
                """
        if input.bam != input.het_vcf:
            cmd += " --het_vcf {input.het_vcf}"
        if input.bam != input.normal:
            cmd += " --normal_bam {input.normal}"
        shell(cmd)


rule IntervalListToBed:
    input:
        get_info_from_metadata("bait_intervals")
    output:
        temp("{sample}_bait_intervals.bed")
    shell:
        "grep -vE '^@' {input} | cut -f1-3 > {output}"


rule SummarizeConsensusSsc:
    input:
        bam = "FilterMolecularConsensusReads/{sample}.ssc_consensus.bam",
    output:
        temp("{sample}.summarize_mutants_ssc.txt"),
    params:
        mut_maf = get_info_from_metadata("mut_maf"),
    shell:
        """
        {SUMMARIZE_MUTANTS} \
            -b {input.bam} \
            -m {params.mut_maf} \
            -o {output}
        """

rule SscReadAlleleFraction:
    input:
        ssc_file = "{sample}.summarize_mutants_ssc.txt"
    output:
        output = "SscReadAlleleFraction/{sample}.ssc_read_af.txt"
    resources:
        mem = 12,
        runtime = 2
    run:
        ssc_df = pd.read_table(input.ssc_file)
        outdf = ssc_df.groupby(['target_site', 'molecule_class']).family_size.sum().unstack().replace(np.nan, 0).reset_index()
        if 'OTHER' not in outdf.columns:
            outdf['OTHER'] = 0
        outdf['VRF'] = outdf['ALT'] / outdf[['ALT', 'REF', 'OTHER']].sum(axis=1)
        outdf.to_csv(output.output, sep='\t', index=None)


rule AggregateAndSummarize:
    input:
        expand("CollectMultipleMetrics/{sample}.alignment_summary_metrics", sample=samples),
        expand("CollectMultipleMetrics/{sample}.insert_size_metrics", sample=samples),
        expand("CollectMultipleMetrics/{sample}.quality_distribution_metrics", sample=samples),
        expand("CollectHsMetrics/{sample}.hs_metrics.txt", sample=samples),
        expand("DuplexMetrics/{sample}.duplex_family_sizes.txt", sample=samples),
        expand("DuplexMetrics/{sample}.duplex_yield_metrics.txt", sample=samples),
        expand("DuplexMetrics/{sample}.family_sizes.txt", sample=samples),
        expand("DuplexMetrics/{sample}.umi_counts.txt", sample=samples),
        expand("FilterConsensusReads/{sample}.dsc_consensus.bam", sample=samples),
        expand("FilterMolecularConsensusReads/{sample}.ssc_consensus.bam", sample=samples),
        expand("Miredas/{sample}.miredas_detect_fingerprint.txt", sample=samples),
        expand("Miredas/{sample}._error_metrics.txt", sample=samples),
        expand("SscReadAlleleFraction/{sample}.ssc_read_af.txt", sample=samples)
    output:
        "summary_metrics.txt"
    run:
        def create_metadata_file(input_arr):
            """ Loop through input array (list of files) and create summary dataframe that
            shows where each file is located """
            results_dict = {}
            for filename in input_arr:
                # files must be in form sample.file_id(ext)?
                for sample_id in samples:
                    if sample_id in filename:
                        corrected = sample_id.replace('.', '_')
                        filename_corrected = filename.replace(sample_id, corrected)
                file_id = re.findall('\.([a-zA-Z_]+)\.?', filename_corrected)[0]
                file_path = os.path.join(os.getcwd(), filename)
                if file_id in results_dict:
                    results_dict[file_id].append(file_path)
                else:
                    results_dict[file_id] = [file_path]
            df = pd.DataFrame.from_dict(results_dict, orient='columns')
            df.insert(loc=0, column='sample_id', value=samples)
            return df


        def get_columns_from_metrics(filename, columns, row=0, **kwargs):
            """ Extract data from metrics files by specifying file, column, and row """
            if next(open(filename)).startswith('#'):
                kwargs['comment'] = '#'

            df = pd.read_csv(filename, sep='\t', **kwargs)

            # Special cases
            if row == -1:
                row = len(df) - 1
            if columns == '*':
                return df.loc[row, :].to_dict()

            return df.loc[row, columns].to_dict()

        Metric = namedtuple('Metric', 'file_id row columns shorthand')
        metrics = [
            Metric(file_id='alignment_summary_metrics',
                   row=2,
                   columns=['TOTAL_READS', 'PF_READS_ALIGNED', 'PF_MISMATCH_RATE', 'PF_INDEL_RATE'],
                   shorthand='asm'),

            Metric(file_id='hs_metrics',
                   row=0,
                   columns=['PCT_SELECTED_BASES', 'MEAN_TARGET_COVERAGE'],
                   shorthand='hsm'),

            Metric(file_id='duplex_yield_metrics',
                   row=-1,
                   columns='*',
                   shorthand='dym'),

            Metric(file_id='_error_metrics',
                   row=0,
                   columns=['perc_error_free_positions', 'errors_per_base_sequenced'],
                   shorthand='mirerr')
        ]
        metrics_df = create_metadata_file(input)
        summary_df = pd.DataFrame()
        summary_df['sample_id'] = metrics_df['sample_id']
        for idx, row in metrics_df.iterrows():
            for metric in metrics:
                metric_dict = get_columns_from_metrics(row[metric.file_id], metric.columns, metric.row)
                metric_dict = {f'{metric.shorthand}.{k.lower()}': v for k, v in metric_dict.items()}
                for k, v in metric_dict.items():
                    summary_df.loc[idx, k] = v
        summary_df = pd.concat([summary_df, metrics_df], axis=1)
        summary_df.to_csv(str(output), sep='\t', index=None)

        # Clean up the mess
        shutil.rmtree("tmp", ignore_errors=True)
