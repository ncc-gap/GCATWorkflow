#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class Melt(stage_task.Stage_task):
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

export LANG=C

mkdir -p {OUTPUT_DIR}
cd {OUTPUT_DIR}

ls {MELT_REFS}/*.zip > mei_list.txt

java \\
    {MELT_JAVA_OPTION} \\
    -jar {MELT_JAR} Single \\
    -a \\
    -h {REFERENCE} \\
    -w {OUTPUT_DIR} \\
    -t mei_list.txt \\
    -n {MELT_BED} \\
    -bamfile {INPUT_CRAM} \\
    -bowtie /tools/bowtie2-2.4.1-linux-x86_64/bowtie2 \\
    -samtools /tools/samtools-1.9/samtools {MELT_OPTION}

"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(input_bams, gcat_conf, run_conf, sample_conf):
    
    STAGE_NAME = "melt"
    CONF_SECTION = STAGE_NAME
    image = gcat_conf.safe_get(CONF_SECTION, "image", "")
    if image != "":
        image = gcat_conf.path_get(CONF_SECTION, "image")

    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": image,
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Melt(params)
    
    output_files = []
    for sample in sample_conf.melt:
        output_files.append("melt/%s/ALU.final_comp.vcf" % (sample))
        output_files.append("melt/%s/HERVK.final_comp.vcf" % (sample))
        output_files.append("melt/%s/LINE1.final_comp.vcf" % (sample))
        output_files.append("melt/%s/SVA.final_comp.vcf" % (sample))
        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_DIR":  "%s/melt/%s" % (run_conf.project_root, sample),
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "MELT_JAR": gcat_conf.get(CONF_SECTION, "melt_jar"),
            "MELT_BED": gcat_conf.get(CONF_SECTION, "melt_bed"),
            "MELT_REFS": gcat_conf.get(CONF_SECTION, "melt_refs"),
            "MELT_OPTION": gcat_conf.get(CONF_SECTION, "melt_option"),
            "MELT_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "melt_java_option")
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    
    return output_files
