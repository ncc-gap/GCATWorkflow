#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

BAM_POSTFIX = ".Aligned.sortedByCoord.out.bam"
BAI_POSTFIX = ".Aligned.sortedByCoord.out.bam.bai"
CHIMERIC_JUNCTION_POSTFIX = ".Chimeric.out.junction"
CHIMERIC_SAM_POSTFIX = ".Chimeric.out.sam"
SJ_TAB_POSTFIX = ".SJ.out.tab.gz"

OUTPUT_BAM_FORMAT = "star/{sample}/{sample}" + BAM_POSTFIX
OUTPUT_CHIMERIC_JUNCTION_FORMAT = "star/{sample}/{sample}" + CHIMERIC_JUNCTION_POSTFIX
OUTPUT_CHIMERIC_SAM_FORMAT = "star/{sample}/{sample}" + CHIMERIC_SAM_POSTFIX
OUTPUT_SJ_TAB_FORMAT = "star/{sample}/{sample}" + SJ_TAB_POSTFIX

class Star_align(stage_task.Stage_task):
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

OUTPUT_PREF={OUTPUT_DIR}/{SAMPLE}
rm -rf {OUTPUT_DIR}/*

# cat fastq
{cat_fastq}

# STAR
if [ -e {FASTQ2} ]; then
  /usr/local/bin/STAR \
    --genomeDir {STAR_REFERENCE} \
    --readFilesIn {FASTQ1} {FASTQ2} \
    --outFileNamePrefix ${{OUTPUT_PREF}}. \
    {STAR_OPTION}
else
  /usr/local/bin/STAR \
    --genomeDir {STAR_REFERENCE} \
    --readFilesIn {FASTQ1} \
    --outFileNamePrefix ${{OUTPUT_PREF}}. \
    {STAR_OPTION}
fi

gzip ${{OUTPUT_PREF}}.SJ.out.tab
rm ${{OUTPUT_PREF}}.SJ.out.tab

# sort
/usr/local/bin/samtools sort \
  -T ${{OUTPUT_PREF}}.Aligned.sortedByCoord.out \
  -@ 6 -m 3G \
  ${{OUTPUT_PREF}}.Aligned.out.bam \
  -O bam > ${{OUTPUT_PREF}}.Aligned.sortedByCoord.out.bam

rm ${{OUTPUT_PREF}}.Aligned.out.bam

# index
/usr/local/bin/samtools index \
  ${{OUTPUT_PREF}}.Aligned.sortedByCoord.out.bam

# remove fastq
{remove_fastq}
"""

def configure(gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "star_alignment"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": gcat_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Star_align(params)
    
    output_bams = {}
    for sample in sample_conf.fastq:
        output_dir = "%s/star/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)
        output_bams[sample] = "%s/%s.Aligned.sortedByCoord.out.bam" % (output_dir, sample)

        paired = len(sample_conf.fastq[sample]) > 1
        cat_fastq = ""
        remove_fastq = ""
        fastq2 = ""
        if len(sample_conf.fastq[sample][0]) == 1:
            fastq1 = sample_conf.fastq[sample][0][0]
            if paired:
                fastq2 = sample_conf.fastq[sample][1][0]
        else:
            cat_fastq += "cat {FASTQ1s} > {OUTPUT_DIR}/1_1_temp.fastq\n".format(
                FASTQ1s = " ".join(sample_conf.fastq[sample][0]),
                OUTPUT_DIR = output_dir
            )
            fastq1 = "{OUTPUT_DIR}/1_1_temp.fastq".format(OUTPUT_DIR = output_dir)
            remove_fastq += "rm -f %s\n" %(fastq1)
            if paired:
                cat_fastq += "if [ -e {FASTQ2_1} ]; then".format(FASTQ2_1=sample_conf.fastq[sample][1][0])
                cat_fastq += "cat {FASTQ2s} > {OUTPUT_DIR}/2_1_temp.fastq\n".format(
                    FASTQ2s = " ".join(sample_conf.fastq[sample][1]),
                    OUTPUT_DIR = output_dir
                )
                cat_fastq += "fi"
                fastq2 = "{OUTPUT_DIR}/2_1_temp.fastq".format(OUTPUT_DIR = output_dir)
                remove_fastq += "rm -f %s\n" %(fastq2)

        if gcat_conf.get(SECTION_NAME, "remove_fastq").lower() == "true":
            rfastq = " ".join(sample_conf.fastq[sample][0])
            if paired:
                rfastq += " " + " ".join(sample_conf.fastq[sample][1])
            if rfastq != "":
                remove_fastq += "rm -f %s\n" % (rfastq)

        arguments = {
            "SAMPLE": sample,
            "cat_fastq": cat_fastq,
            "remove_fastq": remove_fastq,
            "FASTQ1": fastq1,
            "FASTQ2": fastq2,
            "OUTPUT_DIR": output_dir,
            "STAR_REFERENCE": gcat_conf.path_get(SECTION_NAME, "star_genome"),
            "STAR_OPTION": " ".join([
                gcat_conf.get(SECTION_NAME, "star_option"),
                gcat_conf.get(SECTION_NAME, "star_threads_option"),
            ]),
        }
        
        singularity_bind = [
            run_conf.project_root,
            gcat_conf.get(SECTION_NAME, "star_genome"),
        ] + sample_conf.fastq_src[sample][0] + sample_conf.fastq_src[sample][1]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_bams
