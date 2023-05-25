#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class Gridss(stage_task.Stage_task):
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

output_dir=$(dirname {OUTPUT_VCF})
rm -rf ${{output_dir}}
mkdir -p ${{output_dir}}

output_pref=${{output_dir}}/{SAMPLE}

export JAVA_TOOL_OPTIONS="{GRIDSS_JAVA_OPTION}"

/opt/gridss/gridss \\
    -o {OUTPUT_VCF}  \\
    -a ${{output_pref}}.gridss-assembly.bam \\
    -r {REFERENCE}  \\
    -j /opt/gridss/{GRIDSS_JAR} \\
    -w ${{output_dir}} \\
    {GRIDSS_OPTION} \\
    {INPUT_CRAM}

bgzip {BGZIP_OPTION} {OUTPUT_VCF}
tabix {TABIX_OPTION} {OUTPUT_VCF}.gz
rm -f {OUTPUT_VCF}.idx

rm -f ${{output_pref}}.gridss-assembly.bam
rm -f ${{output_dir}}/libsswjni.so
rm -rf ${{output_dir}}/*.working/
"""

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    
    STAGE_NAME = "gridss"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Gridss(params)
    
    output_files = {}
    for sample in sample_conf.gridss:
        output_vcf = "%s/gridss/%s/%s.gridss.vcf" % (run_conf.project_root, sample, sample)
        output_files[sample] = []
        output_files[sample].append(output_vcf + ".gz")
        output_files[sample].append(output_vcf + ".gz.tbi")

        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_VCF":  output_vcf,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "GRIDSS_OPTION": gcat_conf.get(CONF_SECTION, "gridss_option") + " " + gcat_conf.get(CONF_SECTION, "gridss_threads_option"),
            "GRIDSS_JAR": gcat_conf.get(CONF_SECTION, "gridss_jar"),
            "GRIDSS_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "gridss_java_option"),
            "SAMTOOLS_OPTION": gcat_conf.get(CONF_SECTION, "samtools_option") + " " + gcat_conf.get(CONF_SECTION, "samtools_threads_option"),
            "BGZIP_OPTION": gcat_conf.get(CONF_SECTION, "bgzip_option") + " " + gcat_conf.get(CONF_SECTION, "bgzip_threads_option"),
            "TABIX_OPTION": gcat_conf.get(CONF_SECTION, "tabix_option"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)
    
    return output_files
