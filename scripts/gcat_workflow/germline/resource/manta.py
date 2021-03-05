#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class Manta(stage_task.Stage_task):
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

if [ -f {OUTPUT_DIR}/runWorkflow.py ]; then
    rm -rf {OUTPUT_DIR}/*
fi

python /manta/bin/configManta.py \\
    --bam {INPUT_CRAM} \\
    --referenceFasta {REFERENCE} \\
    --runDir {OUTPUT_DIR} {MANTA_CONFIG_OPTION}

python {OUTPUT_DIR}/runWorkflow.py {MANTA_WORKFLOW_OPTION}

"""

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    
    STAGE_NAME = "manta"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Manta(params)
    
    output_files = {}
    for sample in sample_conf.manta:
        output_dir = "%s/manta/%s" % (run_conf.project_root, sample)

        output_files[sample] = []
        output_files[sample].append(output_dir + "/results/variants/candidateSV.vcf.gz")
        output_files[sample].append(output_dir + "/results/variants/candidateSV.vcf.gz.tbi")

        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_DIR": output_dir,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "MANTA_CONFIG_OPTION": gcat_conf.get(CONF_SECTION, "manta_config_option"),
            "MANTA_WORKFLOW_OPTION": gcat_conf.get(CONF_SECTION, "manta_workflow_option") + " " + gcat_conf.get(CONF_SECTION, "manta_workflow_threads_option"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    
    return output_files
