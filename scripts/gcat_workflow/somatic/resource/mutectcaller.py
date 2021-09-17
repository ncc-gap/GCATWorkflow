#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

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

/usr/bin/java \\
  {MUTECT_JAVA_OPTION} \\
  -jar {GATK_JAR} Mutect2 \\
  -I={INPUT_TUMOR_CRAM} -tumor {SAMPLE} {INPUT_NORMAL_CRAM} \\
  -O={OUTPUT_VCF} \\
  -R={REFERENCE} {MUTECT_OPTION}

bgzip {BGZIP_OPTION} {OUTPUT_VCF}
tabix {TABIX_OPTION} {OUTPUT_VCF}.gz
rm -f {OUTPUT_VCF}.idx
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

mkdir -p $(dirname {OUTPUT_VCF})
{PBRUN} mutectcaller \
  --ref {REFERENCE} {MUTECT_OPTION} \\
  --in-tumor-bam {INPUT_TUMOR_CRAM} {INPUT_NORMAL_CRAM} \\
  --tumor-name {SAMPLE} \\
  --out-vcf {OUTPUT_VCF}

bgzip {BGZIP_OPTION} {OUTPUT_VCF}
tabix {TABIX_OPTION} {OUTPUT_VCF}.gz
rm -f {OUTPUT_VCF}.idx
"""

STAGE_NAME = "mutectcaller_parabricks"

def _compatible(input_bams, gcat_conf, run_conf, sample_conf):
    
    CONF_SECTION = "gatk_%s_compatible" % (STAGE_NAME)
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Compatible(params)
    
    output_files = {}
    for (tumor, normal) in sample_conf.mutect_call:
        output_vcf = "%s/mutectcaller/%s/%s.mutectcaller.vcf" % (run_conf.project_root, tumor, tumor)
        output_files[tumor] = [
            output_vcf + ".gz", output_vcf + ".gz.tbi"
        ]
        input_normal_cram = ""
        if normal != None:
            input_normal_cram = "-I %s -normal %s" % (input_bams[normal], normal)

        arguments = {
            "SAMPLE": tumor,
            "INPUT_TUMOR_CRAM": input_bams[tumor],
            "INPUT_NORMAL_CRAM": input_normal_cram,
            "OUTPUT_VCF":  output_vcf,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "GATK_JAR": gcat_conf.get(CONF_SECTION, "gatk_jar"),
            "MUTECT_OPTION": gcat_conf.get(CONF_SECTION, "mutect_option") + " " + gcat_conf.get(CONF_SECTION, "mutect_threads_option"),
            "MUTECT_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "mutect_java_option"),
            "BGZIP_OPTION": gcat_conf.get(CONF_SECTION, "bgzip_option") + " " + gcat_conf.get(CONF_SECTION, "bgzip_threads_option"),
            "TABIX_OPTION": gcat_conf.get(CONF_SECTION, "tabix_option"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if tumor in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[tumor]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = tumor)
    
    return output_files


# merge sorted bams into one and mark duplicate reads with biobambam
def _parabricks(input_bams, gcat_conf, run_conf, sample_conf):
    
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
    
    output_files = {}
    for (tumor, normal) in sample_conf.mutect_call:
        output_dir = "%s/mutectcaller/%s" % (run_conf.project_root, tumor) 
        os.makedirs(output_dir, exist_ok = True)
        output_vcf = "%s/%s.mutectcaller.vcf" % (output_dir, tumor)
        output_files[tumor] = [
            output_vcf + ".gz", output_vcf + ".gz.tbi"
        ]
        
        input_real_path = ""
        if not os.path.islink(input_bams[tumor]):
            input_real_path = input_bams[tumor]
        else:
            for path in sample_conf.bam_import_src[tumor]:
                if not os.path.islink(path) and path.split(".")[-1] in ["bam", "cram"]:
                    input_real_path = path

        input_normal_cram = ""
        if normal != None:
            input_real_path_normal = ""
            if not os.path.islink(input_bams[normal]):
                input_real_path_normal = input_bams[normal]
            else:
                for path in sample_conf.bam_import_src[normal]:
                    if not os.path.islink(path):
                        input_real_path_normal = path
            input_normal_cram = "--in-normal-bam %s --normal-name %s" % (input_real_path_normal, normal)

        arguments = {
            "SAMPLE": tumor,
            "INPUT_TUMOR_CRAM": input_real_path,
            "INPUT_NORMAL_CRAM": input_normal_cram,
            "OUTPUT_VCF": output_vcf,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "MUTECT_OPTION": gcat_conf.get(CONF_SECTION, "mutect_option") + " " + gcat_conf.get(CONF_SECTION, "mutect_threads_option"),
            "PBRUN": gcat_conf.get(CONF_SECTION, "pbrun"),
            "BGZIP_OPTION": gcat_conf.get(CONF_SECTION, "bgzip_option") + " " + gcat_conf.get(CONF_SECTION, "bgzip_threads_option"),
            "TABIX_OPTION": gcat_conf.get(CONF_SECTION, "tabix_option"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if tumor in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[tumor]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = tumor)
    
    return output_files

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    if gcat_conf.safe_get(STAGE_NAME, "gpu_support", "False").lower() == "true":
        return _parabricks(input_bams, gcat_conf, run_conf, sample_conf)
    return _compatible(input_bams, gcat_conf, run_conf, sample_conf)

