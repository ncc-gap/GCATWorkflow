#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class Compatible(stage_task.Stage_task):
    def __init__(self, params):
        super().__init__(params)
        self.shell_script_template = """#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -o errexit
set -o nounset
set -o pipefail
set -x

mkdir -p $(dirname {OUTPUT_FILE_PREFIX})
/usr/bin/java \\
  {MULTIPLE_METRICS_JAVA_OPTION} \\
  -jar {GATK_JAR} \\
  CollectMultipleMetrics \\
  -I {INPUT_CRAM} \\
  -O {OUTPUT_FILE_PREFIX} \\
  -R {REFERENCE} \\
  --PROGRAM CollectAlignmentSummaryMetrics \\
  --PROGRAM CollectInsertSizeMetrics \\
  --PROGRAM QualityScoreDistribution \\
  --PROGRAM MeanQualityByCycle \\
  --PROGRAM CollectBaseDistributionByCycle \\
  --PROGRAM CollectGcBiasMetrics \\
  --PROGRAM CollectSequencingArtifactMetrics \\
  --PROGRAM CollectQualityYieldMetrics {MULTIPLE_METRICS_OPTION}
"""

class Parabricks(stage_task.Stage_task):
    def __init__(self, params):
        super().__init__(params)
        self.shell_script_template = """#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -o errexit
set -o nounset
set -o pipefail
set -x

mkdir -p {OUTPUT_DIR}
{PBRUN} collectmultiplemetrics \\
    --ref {REFERENCE} \\
    --bam {INPUT_CRAM} \\
    --out-qc-metrics-dir {OUTPUT_DIR} {MULTIPLE_METRICS_OPTION}

mv {OUTPUT_DIR}/base_distribution_by_cycle.pdf {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.base_distribution_by_cycle.pdf
mv {OUTPUT_DIR}/base_distribution_by_cycle.png {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.base_distribution_by_cycle.txt
mv {OUTPUT_DIR}/base_distribution_by_cycle.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.base_distribution_by_cycle_metrics
mv {OUTPUT_DIR}/gcbias.pdf {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.gc_bias.pdf
mv {OUTPUT_DIR}/gcbias_0.png {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.gc_bias_0.png
mv {OUTPUT_DIR}/gcbias_detail.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.gc_bias.detail_metrics
mv {OUTPUT_DIR}/gcbias_summary.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.gc_bias.summary_metrics
mv {OUTPUT_DIR}/insert_size.pdf {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.insert_size_histogram.pdf
mv {OUTPUT_DIR}/insert_size.png {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.insert_size_histogram.png
mv {OUTPUT_DIR}/insert_size.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.insert_size_metrics
mv {OUTPUT_DIR}/mean_quality_by_cycle.pdf {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.quality_by_cycle.pdf
mv {OUTPUT_DIR}/mean_quality_by_cycle.png {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.quality_by_cycle.png
mv {OUTPUT_DIR}/mean_quality_by_cycle.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.quality_by_cycle_metrics
mv {OUTPUT_DIR}/quality_yield.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.quality_yield_metrics
mv {OUTPUT_DIR}/qualityscore.pdf {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.quality_distribution.pdf
mv {OUTPUT_DIR}/qualityscore.png {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.quality_distribution.png
mv {OUTPUT_DIR}/qualityscore.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.quality_distribution_metrics
mv {OUTPUT_DIR}/sequencingArtifact.bait_bias_detail_metrics.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.bait_bias_detail_metrics
mv {OUTPUT_DIR}/sequencingArtifact.bait_bias_summary_metrics.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.bait_bias_summary_metrics
mv {OUTPUT_DIR}/sequencingArtifact.error_summary_metrics.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.error_summary_metrics
mv {OUTPUT_DIR}/sequencingArtifact.pre_adapter_detail_metrics.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.pre_adapter_detail_metrics
mv {OUTPUT_DIR}/sequencingArtifact.pre_adapter_summary_metrics.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.pre_adapter_summary_metrics
mv {OUTPUT_DIR}/alignment.txt {OUTPUT_DIR}/{SAMPLE}.collect_multiple_metrics.alignment_summary_metrics
"""

STAGE_NAME = "collect_multiple_metrics"

def _compatible(input_bams, gcat_conf, run_conf, sample_conf):
    
    CONF_SECTION = "gatk_%s_compatible" % (STAGE_NAME)
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Compatible(params)
    
    output_files = {}
    for sample in sample_conf.multiple_metrics:
        output_prefix = "%s/summary/%s/%s.collect_multiple_metrics" % (run_conf.project_root, sample, sample)
        output_files[sample] = [
            output_prefix + ".base_distribution_by_cycle.pdf",
            output_prefix + ".gc_bias.pdf",
            output_prefix + ".insert_size_histogram.pdf",
            output_prefix + ".quality_by_cycle.pdf",
            output_prefix + ".quality_distribution.pdf",
            output_prefix + ".alignment_summary_metrics",
        ]
        arguments = {
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_FILE_PREFIX": output_prefix,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "GATK_JAR": gcat_conf.get(CONF_SECTION, "gatk_jar"),
            "MULTIPLE_METRICS_OPTION": gcat_conf.get(CONF_SECTION, "multiple_metrics_option") + " " + gcat_conf.get(CONF_SECTION, "multiple_metrics_threads_option"),
            "MULTIPLE_METRICS_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "multiple_metrics_java_option"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)
    
    return output_files

def _parabricks(input_bams, gcat_conf, run_conf, sample_conf):

    CONF_SECTION = STAGE_NAME

    image = gcat_conf.safe_get(CONF_SECTION, "image", "")
    singularity_option = gcat_conf.safe_get(CONF_SECTION, "singularity_option", "")
    if image != "":
        image = gcat_conf.path_get(CONF_SECTION, "image")
        singularity_option = gcat_conf.get(CONF_SECTION, "singularity_option")

    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": image,
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": singularity_option
    }
    stage_class = Parabricks(params)
    
    output_files = {}
    for sample in sample_conf.multiple_metrics:
        output_prefix = "%s/summary/%s/%s.collect_multiple_metrics" % (run_conf.project_root, sample, sample)
        output_files[sample] = [
            output_prefix + ".base_distribution_by_cycle.pdf",
            output_prefix + ".gc_bias.pdf",
            output_prefix + ".insert_size_histogram.pdf",
            output_prefix + ".quality_by_cycle.pdf",
            output_prefix + ".quality_distribution.pdf",
            output_prefix + ".alignment_summary_metrics"
        ]
        input_real_path = ""
        if not os.path.islink(input_bams[sample]):
            input_real_path = input_bams[sample]
        else:
            for path in sample_conf.bam_import_src[sample]:
                if not os.path.islink(path) and path.split(".")[-1] in ["bam", "cram"]:
                    input_real_path = path

        arguments = {
            "INPUT_CRAM": input_real_path,
            "OUTPUT_DIR": "%s/summary/%s" % (run_conf.project_root, sample),
            "SAMPLE": sample,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "MULTIPLE_METRICS_OPTION": gcat_conf.get(CONF_SECTION, "multiple_metrics_option") + " " + gcat_conf.get(CONF_SECTION, "multiple_metrics_threads_option"),
            "PBRUN": gcat_conf.get(CONF_SECTION, "pbrun"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)
    
    return output_files

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    if gcat_conf.safe_get(STAGE_NAME, "gpu_support", "False").lower() == "true":
        return _parabricks(input_bams, gcat_conf, run_conf, sample_conf)
    return _compatible(input_bams, gcat_conf, run_conf, sample_conf)

