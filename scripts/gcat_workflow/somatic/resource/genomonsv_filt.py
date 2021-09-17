#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class GenomonSV_filt(stage_task.Stage_task):
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

GenomonSV filt {input_bam} {output_prefix} {reference_genome} {param}
mv {output_prefix}.genomonSV.result.txt {output_prefix}.genomonSV.result.txt.tmp
echo -e "{meta_info}" > {output_prefix}.genomonSV.result.txt
cat {output_prefix}.genomonSV.result.txt.tmp >> {output_prefix}.genomonSV.result.txt
rm -rf {output_prefix}.genomonSV.result.txt.tmp
sv_utils filter {output_prefix}.genomonSV.result.txt {output_prefix}.genomonSV.result.filt.txt.tmp {simple_repeat_file} {sv_utils_param} 

mv {output_prefix}.genomonSV.result.filt.txt.tmp {output_prefix}.genomonSV.result.filt.txt

bgzip {BGZIP_OPTION} {output_prefix}.genomonSV.result.txt
bgzip {BGZIP_OPTION} {output_prefix}.genomonSV.result.filt.txt
"""

def configure(input_bams, sv_merged, gcat_conf, run_conf, sample_conf):
    
    STAGE_NAME = "genomonsv_filt"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = GenomonSV_filt(params)
    
    output_files = {}
    for (tumor, normal, panel) in sample_conf.genomon_sv:
        output_prefix = "{root}/genomonsv/{sample}/{sample}".format(root = run_conf.project_root, sample=tumor)
        output_files[tumor] = output_prefix + ".genomonSV.result.filt.txt.gz"

        filt_param = ""
        if normal != None:
            filt_param = filt_param + " --matched_control_bam " + input_bams[normal]

        if panel != None:
            filt_param = filt_param + " --non_matched_control_junction " + sv_merged[panel]
            if normal != None:
                filt_param = filt_param + " --matched_control_label " + normal

        filt_param = filt_param.lstrip(' ') + ' ' + gcat_conf.get(CONF_SECTION, "params")
        
        simple_repeat_file_path =  gcat_conf.path_get(CONF_SECTION, "simple_repeat_file")
        simple_repeat_file = ""
        if simple_repeat_file_path != "":
            simple_repeat_file = "--simple_repeat_file " + simple_repeat_file_path

        arguments = {
            "input_bam": input_bams[tumor],
            "output_prefix": output_prefix,
            "meta_info": "# genomon_sv: %s" % (gcat_conf.path_get(CONF_SECTION, "image")),
            "reference_genome": gcat_conf.path_get(CONF_SECTION, "reference"),
            "param": filt_param,
            "sv_utils_param": gcat_conf.get(CONF_SECTION, "sv_utils_params"),
            "simple_repeat_file": simple_repeat_file,
            "BGZIP_OPTION": gcat_conf.get(CONF_SECTION, "bgzip_option") + " " + gcat_conf.get(CONF_SECTION, "bgzip_threads_option"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if tumor in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[tumor]
        if normal != None and normal in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[normal]
        if simple_repeat_file_path != "":
            singularity_bind.append(simple_repeat_file_path)

        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = tumor)
    
    return output_files
