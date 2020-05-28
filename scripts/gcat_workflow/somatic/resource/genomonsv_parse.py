#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class genomonsv_parse(stage_task.Stage_task):
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

genomonsv parse {input_bam} {output_prefix} {param}
"""

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    
    STAGE_NAME = "genomonsv-parse"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = genomonsv_parse(params)

    samples = []
    for (tumor, normal, panel) in sample_conf.genomon_sv:
        samples.append(tumor)
        if normal != None:
            samples.append(normal)
        if panel != None:
            samples.extend(sample_conf.control_panel[panel])
    samples = list(set(samples))

    output_files = {}
    for sample in samples:
        output_dir = '%s/genomonsv/%s' % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)
        output_files[sample] = "%s/%s.junction.clustered.bedpe.gz.tbi"% (output_dir, sample)

        arguments = {
            "input_bam": input_bams[sample],
            "output_prefix": "%s/%s"% (output_dir, sample),
            "param": gcat_conf.get(CONF_SECTION, "params"),
        }
        
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    
    return output_files
