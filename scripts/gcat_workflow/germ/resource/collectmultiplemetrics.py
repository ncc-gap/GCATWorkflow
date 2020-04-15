#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class CollectMultipleMetrics(stage_task.Stage_task):
    def __init__(self, params):
        super().__init__(params)
        self.shell_script_template = """
#!/bin/bash
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
java \\
    -jar /gatk/gatk.jar \\
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
    --PROGRAM CollectQualityYieldMetrics

"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(input_bams, gcat_conf, run_conf, sample_conf):
    
    STAGE_NAME = "gatk-collect-multiple-metrics"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = CollectMultipleMetrics(params)
    
    output_files = []
    for sample in sample_conf.multiple_metrics:
        output_prefix = "summary/%s/%s.collect-multiple-metrics" % (sample, sample)
        output_files.append(output_prefix + ".gc_bias.pdf")
        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_FILE_PREFIX":  "%s/%s" % (run_conf.project_root, output_prefix),
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    
    return output_files
