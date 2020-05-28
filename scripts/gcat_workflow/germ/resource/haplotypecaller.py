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
  {HAPLOTYPE_JAVA_OPTION} \\
  -jar {GATK_JAR} HaplotypeCaller \\
  -I={INPUT_CRAM} \\
  -O={OUTPUT_VCF} \\
  -R={REFERENCE} {HAPLOTYPE_OPTION}
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

{PBRUN} haplotypecaller \\
  --ref {REFERENCE} \\
  --in-bam {INPUT_CRAM} \\
  --out-variants {OUTPUT_VCF} {HAPLOTYPE_OPTION} \\
  --tmp-dir {OUTPUT_DIR}/tmp
"""

STAGE_NAME = "haplotypecaller-parabricks"

def _compatible(input_bams, gcat_conf, run_conf, sample_conf):
    
    CONF_SECTION = "gatk-%s-compatible" % (STAGE_NAME)
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Compatible(params)
    
    output_files = {}
    for sample in sample_conf.haplotype_call:
        output_vcf = "%s/haplotypecaller/%s/%s.haplotypecaller.vcf" % (run_conf.project_root, sample, sample)
        output_files[sample] = output_vcf
        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_VCF": output_vcf,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "GATK_JAR": gcat_conf.get(CONF_SECTION, "gatk_jar"),
            "HAPLOTYPE_OPTION": gcat_conf.get(CONF_SECTION, "haplotype_option"),
            "HAPLOTYPE_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "haplotype_java_option")
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    
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
    for sample in sample_conf.haplotype_call:
        output_dir = "%s/haplotypecaller/%s" % (run_conf.project_root, sample) 
        os.makedirs(output_dir, exist_ok = True)
        output_vcf = "%s/%s.haplotypecaller.vcf" % (output_dir, sample)
        output_files[sample] = output_vcf

        input_real_path = ""
        if not os.path.islink(input_bams[sample]):
            input_real_path = input_bams[sample]
        else:
            for path in sample_conf.bam_import_src[sample]:
                if not os.path.islink(path):
                    input_real_path = path

        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_real_path,
            "OUTPUT_VCF": output_vcf,
            "OUTPUT_DIR": output_dir,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "HAPLOTYPE_OPTION": gcat_conf.get(CONF_SECTION, "haplotype_option"),
            "PBRUN": gcat_conf.get(CONF_SECTION, "pbrun"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    
    return output_files

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    if gcat_conf.safe_get(STAGE_NAME, "gpu_support", "False").lower() == "true":
        return _parabricks(input_bams, gcat_conf, run_conf, sample_conf)
    return _compatible(input_bams, gcat_conf, run_conf, sample_conf)

