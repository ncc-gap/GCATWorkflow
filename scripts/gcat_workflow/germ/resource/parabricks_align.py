#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

OUTPUT_FORMAT = "cram/{sample}/{sample}.markdup.bam"

class Compatible(stage_task.Stage_task):
    def __init__(self, params):
        super().__init__(params)
        self.shell_script_template = """
#!/bin/bash
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

/tools/bwa-0.7.15/bwa mem \\
  {BWA_OPTION} \\
  -R "@RG\\tID:{SAMPLE_NAME}\\tPL:{READ_GROUP_PL}\\tLB:{READ_GROUP_LB}\\tSM:{SAMPLE_NAME}\\tPU:{READ_GROUP_PU}" \\
  {REFERENCE} \\
  {INPUT_FASTQ_1} \\
  {INPUT_FASTQ_2} \\
| /usr/bin/java \\
  {GATK_SORT_JAVA_OPTION} \\
  -jar {GATK_JAR} SortSam \\
  {GATK_SORT_OPTION} \\
  -I=/dev/stdin \\
  -O={WORK_DIR}/{SAMPLE_NAME}.bam \\
  --SORT_ORDER=coordinate

/usr/bin/java \\
  {GATK_MARKDUP_JAVA_OPTION} \\
  -jar {GATK_JAR} MarkDuplicates \\
  -I={WORK_DIR}/{SAMPLE_NAME}.bam\\
  -O={OUTPUT_BAM} \\
  -M={OUTPUT_MARKDUP_METRICS} {GATK_MARKDUP_OPTION}

rm -f {WORK_DIR}/{SAMPLE_NAME}.bam
"""

class Parabricks(stage_task.Stage_task):
    def __init__(self, params):
        super().__init__(params)
        self.shell_script_template = """
#!/bin/bash
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

{PBRUN} fq2bam \\
  --bwa-options "{BWA_OPTION}" \\
  "@RG\\tID:{SAMPLE_NAME}\\tPL:{READ_GROUP_PL}\\tLB:{READ_GROUP_LB}\\tSM:{SAMPLE_NAME}\\tPU:{READ_GROUP_PU}" \\
  --ref {REFERENCE} \\
  --in-fq {INPUT_FASTQ_1} {INPUT_FASTQ_2} \\
  --out-bam {OUTPUT_BAM}
#  --tmp-dir /scratch/tmp
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def _compatible(gcat_conf, run_conf, sample_conf):

    STAGE_NAME = "bwa-alignment-parabricks"
    CONF_SECTION = "bwa-alignment-parabricks-compatible"
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Compatible(params)
     
    output_bams = {}
    for sample in sample_conf.fastq:
        output_dir = "%s/cram/%s" % (run_conf.project_root, sample)
        output_bams[sample] = "%s/%s.markdup.bam" % (output_dir, sample)
        
        if len(sample_conf.fastq[sample][0]) == 1:
            fastq1 = sample_conf.fastq[sample][0][0]
        else:
            fastq1 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))

        if len(sample_conf.fastq[sample][1]) == 1:
            fastq2 = sample_conf.fastq[sample][1][0]
        else:
            fastq2 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))

        arguments = {
            "SAMPLE_NAME": sample,
            "INPUT_FASTQ_1": fastq1,
            "INPUT_FASTQ_2": fastq2,
            "OUTPUT_BAM": output_bams[sample],
            "OUTPUT_MARKDUP_METRICS": "%s/%s.markdup.metrics" % (output_dir, sample),
            "WORK_DIR": output_dir,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "BWA_OPTION": gcat_conf.get(CONF_SECTION, "bwa_option"),
            "READ_GROUP_PL": gcat_conf.get(CONF_SECTION, "read_group_pl"),
            "READ_GROUP_LB": gcat_conf.get(CONF_SECTION, "read_group_lb"),
            "READ_GROUP_PU": gcat_conf.get(CONF_SECTION, "read_group_pu"),
            "GATK_JAR": gcat_conf.get(CONF_SECTION, "gatk_jar"),
            "GATK_SORT_OPTION": gcat_conf.get(CONF_SECTION, "gatk_sort_option"),
            "GATK_SORT_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "gatk_sort_java_option"),
            "GATK_MARKDUP_OPTION": gcat_conf.get(CONF_SECTION, "gatk_markdup_option"),
            "GATK_MARKDUP_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "gatk_markdup_java_option"),
        }
        
        singularity_bind = [
            run_conf.project_root,
            os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference")),
        ] + sample_conf.fastq_src[sample]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_bams

def _parabricks(gcat_conf, run_conf, sample_conf):

    STAGE_NAME = "bwa-alignment-parabricks"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": "",
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": ""
    }
    stage_class = Parabricks(params)
     
    output_bams = {}
    for sample in sample_conf.fastq:
        output_dir = "%s/cram/%s" % (run_conf.project_root, sample)
        output_bams[sample] = "%s/%s.markdup.bam" % (output_dir, sample)
        
        if len(sample_conf.fastq[sample][0]) == 1:
            fastq1 = sample_conf.fastq[sample][0][0]
        else:
            fastq1 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))

        if len(sample_conf.fastq[sample][1]) == 1:
            fastq2 = sample_conf.fastq[sample][1][0]
        else:
            fastq2 = "'<cat %s'" % (" ".join(sample_conf.fastq[sample][0]))

        arguments = {
            "SAMPLE_NAME": sample,
            "INPUT_FASTQ_1": fastq1,
            "INPUT_FASTQ_2": fastq2,
            "OUTPUT_BAM": output_bams[sample],
            "PBRUN": gcat_conf.path_get(CONF_SECTION, "pbrun"),
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "BWA_OPTION": gcat_conf.get(CONF_SECTION, "bwa_option"),
            "READ_GROUP_PL": gcat_conf.get(CONF_SECTION, "read_group_pl"),
            "READ_GROUP_LB": gcat_conf.get(CONF_SECTION, "read_group_lb"),
            "READ_GROUP_PU": gcat_conf.get(CONF_SECTION, "read_group_pu"),
        }
        
        singularity_bind = []
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_bams

def configure(gcat_conf, run_conf, sample_conf):
    if gcat_conf.safe_get("bwa-alignment-parabricks", "gpu_support", "False").lower() == "true":
        return _parabricks(gcat_conf, run_conf, sample_conf)
    return _compatible(gcat_conf, run_conf, sample_conf)

