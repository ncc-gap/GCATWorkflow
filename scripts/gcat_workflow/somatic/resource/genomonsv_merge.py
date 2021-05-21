#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class GenomonSV_merge(stage_task.Stage_task):
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

GenomonSV merge {control_info} {merge_output_file} {param}
"""

def configure(gcat_conf, run_conf, sample_conf):
    
    STAGE_NAME = "genomonsv_merge"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = GenomonSV_merge(params)
    
    bedpe_dir = run_conf.project_root + '/genomonsv/non_matched_control_panel'
    control_info_dir = run_conf.project_root + '/genomonsv/control_panel'
    os.makedirs(bedpe_dir, exist_ok=True)
    os.makedirs(control_info_dir, exist_ok=True)

    output_files = {}
    for (tumor, normal, panel) in sample_conf.genomon_sv:
        if panel == None:
            continue
        if panel in output_files:
            continue
            
        output_file = "%s/%s.merged.junction.control.bedpe.gz.tbi" % (bedpe_dir, panel)
        output_files[panel] = output_file
        control_conf = "%s/%s.control_info.txt" % (control_info_dir, panel)
        with open(control_conf,  "w") as out_handle:
            for sample in sample_conf.control_panel[panel]:
                out_handle.write(sample+ "\t"+ run_conf.project_root + "/genomonsv/"+ sample +"/"+ sample+ "\n")

        arguments = {
            "control_info": control_conf,
            "merge_output_file": output_file,
            "param": gcat_conf.get(CONF_SECTION, "params"),
        }
       
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = panel)
    
    return output_files
