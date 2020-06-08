#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

class Star_fusion(stage_task.Stage_task):
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

/usr/local/src/STAR-Fusion/STAR-Fusion --genome_lib_dir {REFERENCE} \
  -J {CHIMERIC_OUT_JUNCTION} \
  --output_dir {OUTPUT_DIR} \
  {STAR_FUSION_OPTIONS}

rm -rf {OUTPUT_DIR}/_starF_checkpoints
"""

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "star_fusion"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": gcat_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Star_fusion(params)
    
    output_files = {}
    for sample in sample_conf.star_fusion:
        output_dir = "%s/star_fusion/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)    
        output_files[sample] = "%s/star-fusion.fusion_predictions.abridged.tsv" % (output_dir)

        arguments = {
            "CHIMERIC_OUT_JUNCTION": input_bams[sample],
            "OUTPUT_DIR": output_dir,
            "REFERENCE": gcat_conf.path_get(SECTION_NAME, "star_genome"),
            "STAR_FUSION_OPTIONS": gcat_conf.get(SECTION_NAME, "star_fusion_option")
        }
       
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
        
        singularity_bind.append(os.path.dirname(gcat_conf.path_get(SECTION_NAME, "star_genome")))
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)

    return output_files
