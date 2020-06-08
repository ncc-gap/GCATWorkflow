#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

OUTPUT_FORMAT = "fusionfusion/{sample}/{sample}.genomonFusion.result.filt.txt"

class Fusionfusion(stage_task.Stage_task):
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

if [ "_{MERGED_COUNT}" != "_" ]; then
    OPTION="{OPTION} --pooled_control_file {MERGED_COUNT}"
fi
/usr/local/bin/fusionfusion --star {INPUT} --out {OUTPUT_DIR} --reference_genome {REFERENCE} {OPTION}

mv {OUTPUT_DIR}/star.fusion.result.txt {OUTPUT_DIR}/{SAMPLE}.star.fusion.result.txt
mv {OUTPUT_DIR}/fusion_fusion.result.txt {OUTPUT_DIR}/{SAMPLE}.genomonFusion.result.txt

/usr/local/bin/fusion_utils filt {OUTPUT_DIR}/{SAMPLE}.genomonFusion.result.txt {OUTPUT_DIR}/{SAMPLE}.fusion_fusion.result.filt.txt {FILT_OPTION}

mv {OUTPUT_DIR}/{SAMPLE}.fusion_fusion.result.filt.txt {OUTPUT_DIR}/{SAMPLE}.genomonFusion.result.filt.txt
"""

def configure(input_counts, input_merges, gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "fusionfusion"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": gcat_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Fusionfusion(params)
    output_files = {}
    for (sample, panel) in sample_conf.fusionfusion:
        output_dir = "%s/fusionfusion/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)
        output_files[sample] = run_conf.project_root + "/" + OUTPUT_FORMAT.format(sample = sample)
        
        merged_count = ""
        if panel != None:
            merged_count = "%s/%s" % (run_conf.project_root, input_merges[panel])
            
        arguments = {
            "SAMPLE": sample,
            "INPUT": "%s" % (input_counts[sample]),
            "OUTPUT_DIR": output_dir,
            "OPTION": gcat_conf.get(SECTION_NAME, "fusionfusion_option"),
            "FILT_OPTION": gcat_conf.get(SECTION_NAME, "filt_option"),
            "REFERENCE": gcat_conf.get(SECTION_NAME, "reference"),
            "MERGED_COUNT": merged_count
        }
       
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
        
        singularity_bind.append(os.path.dirname(gcat_conf.path_get(SECTION_NAME, "reference")))
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)

    return output_files
