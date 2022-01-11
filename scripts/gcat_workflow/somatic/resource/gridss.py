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
    {INPUT_NORMAL_CRAM} {INPUT_TUMOR_CRAM}

if [ "{INPUT_NORMAL_CRAM}" != "" ]; then
    Rscript /opt/gridss/gridss_somatic_filter \\
        -i {OUTPUT_VCF} \\
        -o {OUTPUT_VCF_SOMATIC} \\
        {FILTER_OPTION} --normalordinal 1 --tumourordinal 2 --scriptdir /opt/gridss/
fi

bgzip {BGZIP_OPTION} {OUTPUT_VCF}
tabix {TABIX_OPTION} {OUTPUT_VCF}.gz
rm -f {OUTPUT_VCF}.idx

rm -f ${{output_pref}}.gridss-assembly.bam
rm -rf ${{output_pref}}.gridss-assembly.bam.gridss.working/
rm -rf ${{output_pref}}.gridss.vcf.gridss.working/
rm -rf ${{output_pref}}.normal.temp.bam.gridss.working/
rm -rf ${{output_pref}}.tumor.temp.bam.gridss.working/

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
    for (tumor, normal) in sample_conf.gridss:
        output_dir = "%s/gridss/%s" % (run_conf.project_root, tumor)
        output_vcf = "%s/%s.gridss.vcf" % (output_dir, tumor)
        output_somatic_vcf = "%s/%s.gridss.somatic.vcf" % (output_dir, tumor)
        output_files[tumor] = []
        output_files[tumor].append(output_vcf + ".gz")
        output_files[tumor].append(output_vcf + ".gz.tbi")

        input_normal_cram = ""
        if normal != None:
            input_normal_cram = input_bams[normal]

        filter_option = ""
        if gcat_conf.safe_get(CONF_SECTION, "gridss_fulloutput_option", "False").lower() == "true":
            filter_option = "--fulloutput %s/%s%s" % (
                output_dir, tumor, gcat_conf.get(CONF_SECTION, "gridss_fulloutput_postfix")
            )

        arguments = {
            "SAMPLE": tumor,
            "INPUT_TUMOR_CRAM": input_bams[tumor],
            "INPUT_NORMAL_CRAM": input_normal_cram,
            "OUTPUT_VCF":  output_vcf,
            "OUTPUT_VCF_SOMATIC":  output_somatic_vcf,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "GRIDSS_OPTION": gcat_conf.get(CONF_SECTION, "gridss_option") + " " + gcat_conf.get(CONF_SECTION, "gridss_threads_option"),
            "GRIDSS_JAR": gcat_conf.get(CONF_SECTION, "gridss_jar"),
            "GRIDSS_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "gridss_java_option"),
            "SAMTOOLS_OPTION": gcat_conf.get(CONF_SECTION, "samtools_option") + " " + gcat_conf.get(CONF_SECTION, "samtools_threads_option"),
            "BGZIP_OPTION": gcat_conf.get(CONF_SECTION, "bgzip_option") + " " + gcat_conf.get(CONF_SECTION, "bgzip_threads_option"),
            "TABIX_OPTION": gcat_conf.get(CONF_SECTION, "tabix_option"),
            "FILTER_OPTION": filter_option,
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if tumor in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[tumor]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = tumor)
    
    return output_files
