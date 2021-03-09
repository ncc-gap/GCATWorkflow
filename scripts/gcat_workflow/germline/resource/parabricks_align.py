#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

OUTPUT_FORMAT = "cram/{sample}/{sample}.markdup.bam"

class Compatible(stage_task.Stage_task):
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

rm -rf {WORK_DIR}/*

arrays=(
{ARRAYS})

SORTED_BAMS=
REMOVE_BAMS=

for i in "${{arrays[@]}}"
do
    array=(${{i[@]}})
    NUM=${{array[0]}}
    INPUT_FASTQ_1=${{array[1]}}
    INPUT_FASTQ_2=${{array[2]}}
    RG=${{array[3]}}

    /tools/bwa-0.7.15/bwa mem \\
      {BWA_OPTION} \\
      -R "$RG" \\
      {REFERENCE} \\
      $INPUT_FASTQ_1 \\
      $INPUT_FASTQ_2 \\
    | /usr/bin/java \\
      {GATK_SORT_JAVA_OPTION} \\
      -jar {GATK_JAR} SortSam \\
      {GATK_SORT_OPTION} \\
      -I=/dev/stdin \\
      -O={WORK_DIR}/{SAMPLE_NAME}_$NUM.bam \\
      --SORT_ORDER=coordinate
      
    SORTED_BAMS=$SORTED_BAMS" -I={WORK_DIR}/{SAMPLE_NAME}_$NUM.bam"
    REMOVE_BAMS=$REMOVE_BAMS" {WORK_DIR}/{SAMPLE_NAME}_$NUM.bam"
done

/usr/bin/java \\
  {GATK_MARKDUP_JAVA_OPTION} \\
  -jar {GATK_JAR} MarkDuplicates \\
  $SORTED_BAMS \\
  -O={OUTPUT_BAM} \\
  -M={OUTPUT_MARKDUP_METRICS} {GATK_MARKDUP_OPTION}

rm -f {WORK_DIR}/{SAMPLE_NAME}.bam
rm -f $REMOVE_BAMS
"""

class Parabricks(stage_task.Stage_task):
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

rm -rf {OUTPUT_DIR}/*

{PBRUN} fq2bam \\
  --ref {REFERENCE} \\
  {INPUT} \\
  --bwa-options "{BWA_OPTION}" \\
  --out-bam {OUTPUT_BAM} \\
  --tmp-dir {OUTPUT_DIR}/tmp
"""

STAGE_NAME = "bwa_alignment_parabricks"

def _compatible(gcat_conf, run_conf, sample_conf):

    CONF_SECTION = "gatk_%s_compatible" % (STAGE_NAME)
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
        os.makedirs(output_dir, exist_ok = True)
        output_bams[sample] = "%s/%s.markdup.bam" % (output_dir, sample)
        
        readgroups = open(sample_conf.readgroup[sample]).readlines()
        arrays = ""
        for i in range(len(sample_conf.fastq[sample][0])):
            arrays += '"%d %s %s %s"\n' % (i, sample_conf.fastq[sample][0][i], sample_conf.fastq[sample][1][i], readgroups[i].rstrip())
            
        arguments = {
            "SAMPLE_NAME": sample,
            "OUTPUT_BAM": output_bams[sample],
            "OUTPUT_MARKDUP_METRICS": "%s/%s.markdup.metrics" % (output_dir, sample),
            "WORK_DIR": output_dir,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "BWA_OPTION": gcat_conf.get(CONF_SECTION, "bwa_option") + " " + gcat_conf.get(CONF_SECTION, "bwa_threads_option"),
            "GATK_JAR": gcat_conf.get(CONF_SECTION, "gatk_jar"),
            "GATK_SORT_OPTION": gcat_conf.get(CONF_SECTION, "gatk_sort_option"),
            "GATK_SORT_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "gatk_sort_java_option"),
            "GATK_MARKDUP_OPTION": gcat_conf.get(CONF_SECTION, "gatk_markdup_option"),
            "GATK_MARKDUP_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "gatk_markdup_java_option"),
            "ARRAYS": arrays
        }
        
        singularity_bind = [
            run_conf.project_root,
            os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference")),
        ] + sample_conf.fastq_src[sample][0] + sample_conf.fastq_src[sample][1] + sample_conf.readgroup_src[sample]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_bams

def _parabricks(gcat_conf, run_conf, sample_conf):

    CONF_SECTION = STAGE_NAME

    image = gcat_conf.safe_get(CONF_SECTION, "image", "")
    singularity_option = gcat_conf.safe_get(CONF_SECTION, "singularity_option", "")
    if image != "":
        image = gcat_conf.path_get(CONF_SECTION, "image")
        singularity_option = gcat_conf.get(CONF_SECTION, "singularity_option")

    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": image,
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": singularity_option
    }
    stage_class = Parabricks(params)
    output_bams = {}
    for sample in sample_conf.fastq:
        output_dir = "%s/cram/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok = True)
        output_bams[sample] = "%s/%s.markdup.bam" % (output_dir, sample)
        
        readgroups = open(sample_conf.readgroup[sample]).readlines()
        input_params = ""
        bind_fastqs = []
        for i in range(len(sample_conf.fastq[sample][0])):
            fastq1 = ""
            fastq2 = ""
            for path in sample_conf.fastq_src[sample][i]:
                if not os.path.islink(path):
                    bind_fastqs.append(path)
                    if fastq1 == "":
                        fastq1 = path
                    else:
                        fastq2 = path
            input_params += ' --in-fq %s %s "%s"' % (fastq1, fastq2, readgroups[i].rstrip())
            
        arguments = {
            "SAMPLE_NAME": sample,
            "INPUT": input_params,
            "OUTPUT_BAM": output_bams[sample],
            "OUTPUT_DIR":  output_dir,
            "PBRUN": gcat_conf.get(CONF_SECTION, "pbrun"),
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "BWA_OPTION": gcat_conf.get(CONF_SECTION, "bwa_option") + " " + gcat_conf.get(CONF_SECTION, "bwa_threads_option"),
        }
        
        singularity_bind = [
            run_conf.project_root,
            os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference")),
        ] + bind_fastqs
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_bams

def configure(gcat_conf, run_conf, sample_conf):
    if gcat_conf.safe_get(STAGE_NAME, "gpu_support", "False").lower() == "true":
        return _parabricks(gcat_conf, run_conf, sample_conf)
    return _compatible(gcat_conf, run_conf, sample_conf)

