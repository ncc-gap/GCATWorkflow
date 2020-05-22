#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

class IR_count(stage_task.Stage_task):
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

intron_retention_utils simple_count \
  {INPUT_BAM} \
  {OUTPUT_DIR}/{SAMPLE}.genomonIR.result.txt \
  {OPTION}
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(input_bams, gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "intron_retention"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": gcat_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = IR_count(params)
    
    output_files = []
    for sample in sample_conf.ir_count:
        output_dir = "%s/ir_count/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)  
        output_file = "ir_count/{sample}/{sample}.genomonIR.result.txt".format(sample = sample)
        output_files.append(output_file)
        
        arguments = {
            "SAMPLE": sample,
            "INPUT_BAM": input_bams[sample],
            "OUTPUT_DIR": output_dir,
            "OPTION": gcat_conf.get(SECTION_NAME, "params")
        }
       
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)

    return output_files
