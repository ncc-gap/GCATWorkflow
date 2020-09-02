#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

class Join(stage_task.Stage_task):
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

{RM_FASTQS}
{RM_BAMS}
touch {JOIN_FILE}
"""

def configure(samples, fastq_files, bam_files, gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "join"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": "",
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": ""
    }
    stage_class = Join(params)
    remove_fastq = gcat_conf.get(SECTION_NAME, "remove_fastq").lower() == "true"
    remove_bam = gcat_conf.get(SECTION_NAME, "remove_bam").lower() == "true" 
    #print([remove_fastq,remove_bam])
    output_files = {}
    for sample in samples:
        output_dir = "%s/join/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)
        output_file = "%s/join.txt" % (output_dir)
        output_files[sample] = output_file

        list_fastq = []
        if sample in fastq_files:
            for li in fastq_files[sample]:
                list_fastq.extend(li)
        rm_fastqs = ""
        if remove_fastq and len(list_fastq) > 0:
            rm_fastqs = "rm -f " +  " ".join(list_fastq)

        list_bam = []
        if sample in bam_files:
            list_bam.append(bam_files[sample])
            list_bam.append(bam_files[sample] + ".bai")
        rm_bams = ""
        if remove_bam and len(list_bam) > 0:
            rm_bams = "rm -f " +  " ".join(list_bam)

        arguments = {
            "RM_FASTQS": rm_fastqs,
            "RM_BAMS": rm_bams,
            "JOIN_FILE": output_file
        }
    
        singularity_bind = [run_conf.project_root]
    
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)

    return output_files
