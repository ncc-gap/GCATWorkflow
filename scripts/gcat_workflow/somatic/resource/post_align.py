#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

OUTPUT_FORMAT = "cram/{sample}/{sample}.markdup.cram"
BAM_POSTFIX = ".markdup.cram"
BAI_POSTFIX = ".markdup.cram.crai"
    
class PostBwa(stage_task.Stage_task):
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

/tools/samtools-1.9/samtools view \\
  {SAMTOOLS_VIEW_OPTION} \\
  -C \\
  -T {REFERENCE} \\
  -o {OUTPUT_CRAM} \\
  {INPUT_BAM}

/tools/samtools-1.9/samtools index \\
  {SAMTOOLS_INDEX_OPTION} \\
  {OUTPUT_CRAM}
  
rm -f {INPUT_BAM}
rm -f {INPUT_BAI}
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(aligned_bams, gcat_conf, run_conf, sample_conf):

    STAGE_NAME = "post-bwa-alignment-parabricks"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = PostBwa(params)
     
    output_crams = {}
    for sample in aligned_bams:
        output_dir = "%s/cram/%s" % (run_conf.project_root, sample)
        output_crams[sample] = "%s/%s.markdup.cram" % (output_dir, sample)
        
        arguments = {
            "INPUT_BAM": aligned_bams[sample],
            "INPUT_BAI": aligned_bams[sample] + ".bai",
            "OUTPUT_CRAM": output_crams[sample],
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "SAMTOOLS_VIEW_OPTION": gcat_conf.get(CONF_SECTION, "samtools_view_option"),
            "SAMTOOLS_INDEX_OPTION": gcat_conf.get(CONF_SECTION, "samtools_index_option")
        }
        
        singularity_bind = [
            run_conf.project_root,
            os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference")),
        ]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_crams
