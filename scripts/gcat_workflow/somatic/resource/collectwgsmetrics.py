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

mkdir -p $(dirname {OUTPUT_FILE})
/usr/bin/java \\
  {WGS_METRICS_JAVA_OPTION} \\
  -jar {GATK_JAR} \\
  CollectWgsMetrics \\
  -I {INPUT_CRAM} \\
  -O {OUTPUT_FILE} \\
  -R {REFERENCE} {WGS_METRICS_OPTION}
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

mkdir -p $(dirname {OUTPUT_FILE})
{PBRUN} bammetrics \\
  --ref {REFERENCE} \\
  --bam {INPUT_CRAM} \\
  --out-metrics-file {OUTPUT_FILE} {WGS_METRICS_OPTION}
"""

STAGE_NAME = "collect_wgs_metrics"

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
    for sample in sample_conf.wgs_metrics:
        output_txt = "%s/summary/%s/%s.collect_wgs_metrics.txt" % (run_conf.project_root, sample, sample)
        output_files[sample] = output_txt
        arguments = {
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_FILE":  output_txt,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "GATK_JAR": gcat_conf.get(CONF_SECTION, "gatk_jar"),
            "WGS_METRICS_OPTION": gcat_conf.get(CONF_SECTION, "wgs_metrics_option") + " " + gcat_conf.get(CONF_SECTION, "wgs_metrics_threads_option"),
            "WGS_METRICS_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "wgs_metrics_java_option"),
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
    for sample in sample_conf.wgs_metrics:
        output_txt = "%s/summary/%s/%s.collect_wgs_metrics.txt" % (run_conf.project_root, sample, sample)
        output_files[sample] = output_txt

        input_real_path = ""
        if not os.path.islink(input_bams[sample]):
            input_real_path = input_bams[sample]
        else:
            for path in sample_conf.bam_import_src[sample]:
                if not os.path.islink(path) and path.split(".")[-1] in ["bam", "cram"]:
                    input_real_path = path

        arguments = {
            "INPUT_CRAM": input_real_path,
            "OUTPUT_FILE": output_txt,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "WGS_METRICS_OPTION": gcat_conf.get(CONF_SECTION, "wgs_metrics_option") + " " + gcat_conf.get(CONF_SECTION, "wgs_metrics_threads_option"),
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

