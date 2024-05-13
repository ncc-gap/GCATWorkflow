#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class Hs_metrics(stage_task.Stage_task):
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

mkdir -p $(dirname {OUTPUT_FILE})
/usr/bin/java \\
  {HS_METRICS_JAVA_OPTION} \\
  -jar {GATK_JAR} \\
  CollectHsMetrics \\
  -I {INPUT_CRAM} \\
  -O {OUTPUT_FILE} \\
  -R {REFERENCE} \\
  --BAIT_INTERVALS {BAIT_INTERVALS} \\
  --TARGET_INTERVALS {TARGET_INTERVALS} {HS_METRICS_OPTION}
"""


STAGE_NAME = "collect_hs_metrics"

def configure(input_bams, gcat_conf, run_conf, sample_conf):

    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Hs_metrics(params)
    
    output_files = {}
    for sample in sample_conf.hs_metrics:
        output_txt = "%s/summary/%s/%s.collect_hs_metrics.txt" % (run_conf.project_root, sample, sample)
        output_files[sample] = output_txt
        arguments = {
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_FILE":  output_txt,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "BAIT_INTERVALS": gcat_conf.path_get(CONF_SECTION, "bait_intervals"),
            "TARGET_INTERVALS": gcat_conf.path_get(CONF_SECTION, "target_intervals"),
            "GATK_JAR": gcat_conf.get(CONF_SECTION, "gatk_jar"),
            "HS_METRICS_OPTION": gcat_conf.get(CONF_SECTION, "hs_metrics_option") + " " + gcat_conf.get(CONF_SECTION, "hs_metrics_threads_option"),
            "HS_METRICS_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "hs_metrics_java_option"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)
    
    return output_files

