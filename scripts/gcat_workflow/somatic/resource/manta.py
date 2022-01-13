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

if [ "{INPUT_NORMAL_CRAM}" != "" ]; then
    python /manta/bin/configManta.py \\
      --normalBam {INPUT_NORMAL_CRAM} \\
      --tumorBam {INPUT_TUMOR_CRAM} \\
      --referenceFasta {REFERENCE} \\
      --runDir {OUTPUT_DIR} {MANTA_CONFIG_OPTION}
else
    python /manta/bin/configManta.py \\
      --tumorBam {INPUT_TUMOR_CRAM} \\
      --referenceFasta {REFERENCE} \\
      --runDir {OUTPUT_DIR} {MANTA_CONFIG_OPTION}
fi

python {OUTPUT_DIR}/runWorkflow.py {MANTA_WORKFLOW_OPTION}

rm -rf {OUTPUT_DIR}/workspace/
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
    for (tumor, normal) in sample_conf.manta:
        output_dir = "%s/manta/%s" % (run_conf.project_root, tumor)
        
        output_files[tumor] = []
        output_files[tumor].append(output_dir + "/results/variants/candidateSV.vcf.gz")
        output_files[tumor].append(output_dir + "/results/variants/candidateSV.vcf.gz.tbi")
        
        input_normal_cram = ""
        if normal != None:
            input_normal_cram = input_bams[normal]

        arguments = {
            "SAMPLE": tumor,
            "INPUT_TUMOR_CRAM": input_bams[tumor],
            "INPUT_NORMAL_CRAM": input_normal_cram,
            "OUTPUT_DIR":  "%s/manta/%s" % (run_conf.project_root, tumor),
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "MANTA_CONFIG_OPTION": gcat_conf.get(CONF_SECTION, "manta_config_option"),
            "MANTA_WORKFLOW_OPTION": gcat_conf.get(CONF_SECTION, "manta_workflow_option") + " " + gcat_conf.get(CONF_SECTION, "manta_workflow_threads_option"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if tumor in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[tumor]
        if normal in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[normal]
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = tumor)
    
    return output_files
