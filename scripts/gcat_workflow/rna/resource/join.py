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

if [ "{BAM_TOCRAM}" = "T" ]
then
  samtools view -C {OPTION} -T {REFERENCE} {INPUT_BAM} -o {OUTPUT_CRAM}
  samtools index {OUTPUT_CRAM}
done

{RM_BAMS}
touch {JOIN_FILE}
"""

def configure(samples, import_samples, fastq_files, bam_files, gcat_conf, run_conf, sample_conf):
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
    remove_bam = gcat_conf.get(SECTION_NAME, "remove_bam").lower() == "true" 
    bam_tocram = gcat_conf.get(SECTION_NAME, "bam_tocram").lower() == "true" 

    output_files = {}
    for sample in samples:
        output_dir = "%s/join/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)
        output_file = "%s/join.txt" % (output_dir)
        output_files[sample] = output_file

        list_bam = []
        if sample in bam_files:
            list_bam.append(bam_files[sample])
            list_bam.append(bam_files[sample] + ".bai")
        rm_bams = ""
        if remove_bam and len(list_bam) > 0:
            rm_bams = "rm -f " +  " ".join(list_bam)

        flg_bamtocram = "F"
        reference  = gcat_conf.path_get(SECTION_NAME, "reference")
        input_bam = ""
        output_cram = ""
        if bam_tocram and sample in bam_files and not sample in import_samples:
            flg_bamtocram = "T"
            output_dir = "%s/cram/%s" % (run_conf.project_root, sample)
            os.makedirs(output_dir, exist_ok=True)
            input_bam = bam_files[sample]
            (basename, ext) = os.path.splitext(os.path.basename(input_bam))
            output_cram = "%s/%s.cram" % (output_dir, basename)

        arguments = {
            "BAM_TOCRAM": flg_bamtocram,
            "INPUT_BAM": input_bam,
            "OUTPUT_CRAM": output_cram,
            "REFERENCE": reference,
            "OPTION": " ".join([
                gcat_conf.get(SECTION_NAME, "samtools_option"),
                gcat_conf.get(SECTION_NAME, "samtools_threads_option"),
            ]),
            "RM_BAMS": rm_bams,
            "JOIN_FILE": output_file
        }
    
        singularity_bind = [run_conf.project_root, reference]
    
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)

    return output_files
