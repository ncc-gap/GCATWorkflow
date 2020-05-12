#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

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

{PBRUN} mutectcaller \
  --ref {REFERENCE} {HAPLOTYPE_OPTION}
  --in-tumor-bam {INPUT_CRAM} \
  --tumor-name {SAMPLE} \
  --out-vcf {OUTPUT_VCF}
  #--tmp-dir /scratch/tmp
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def _compatible(input_bams, gcat_conf, run_conf, sample_conf):
    
    STAGE_NAME = "gatk-haplotypecaller-parabricks"
    CONF_SECTION = "gatk-haplotypecaller-parabricks-compatible"
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Compatible(params)
    
    output_files = []
    for sample in sample_conf.haplotype_call:
        output_vcf = "haplotypecaller/%s/%s.gatk-hc.vcf" % (sample, sample)
        output_files.append(output_vcf)
        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_VCF":  "%s/%s" % (run_conf.project_root, output_vcf),
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
    
    STAGE_NAME = "gatk-haplotypecaller-parabricks"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": "",
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": ""
    }
    stage_class = Parabricks(params)
    
    output_files = []
    for sample in sample_conf.haplotype_call:
        output_vcf = "haplotypecaller/%s/%s.gatk-hc.vcf" % (sample, sample)
        output_files.append(output_vcf)
        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_VCF":  "%s/%s" % (run_conf.project_root, output_vcf),
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "HAPLOTYPE_OPTION": gcat_conf.get(CONF_SECTION, "haplotype_option"),
            "PBRUN": gcat_conf.path_get(CONF_SECTION, "pbrun"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    
    return output_files

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    if gcat_conf.safe_get("gatk-haplotypecaller-parabricks", "gpu_support", "False").lower() == "true":
        return _parabricks(input_bams, gcat_conf, run_conf, sample_conf)
    return _compatible(input_bams, gcat_conf, run_conf, sample_conf)

