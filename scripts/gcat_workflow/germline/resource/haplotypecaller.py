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
  -O={OUTPUT_AUTOSOME_GVCF} \\
  -L={INTERVAL_AUTOSOME} \\
  -R={REFERENCE} {HAPLOTYPE_OPTION_AUTOSOME} 

bgzip {BGZIP_OPTION} {OUTPUT_AUTOSOME_GVCF}
tabix {TABIX_OPTION} {OUTPUT_AUTOSOME_GVCF}.gz
rm -f {OUTPUT_AUTOSOME_GVCF}.idx

/usr/bin/java \\
  {HAPLOTYPE_JAVA_OPTION} \\
  -jar {GATK_JAR} HaplotypeCaller \\
  -I={INPUT_CRAM} \\
  -O={OUTPUT_PAR_GVCF} \\
  -L={INTERVAL_PAR} \\
  -R={REFERENCE} {HAPLOTYPE_OPTION_PAR} 

bgzip {BGZIP_OPTION} {OUTPUT_PAR_GVCF}
tabix {TABIX_OPTION} {OUTPUT_PAR_GVCF}.gz
rm -f {OUTPUT_PAR_GVCF}.idx

/usr/bin/java \\
  {HAPLOTYPE_JAVA_OPTION} \\
  -jar {GATK_JAR} HaplotypeCaller \\
  -I={INPUT_CRAM} \\
  -O={OUTPUT_CHRX_FEMALE_GVCF} \\
  -L={INTERVAL_CHRX} \\
  -R={REFERENCE} {HAPLOTYPE_OPTION_CHRX_FEMALE} 

bgzip {BGZIP_OPTION} {OUTPUT_CHRX_FEMALE_GVCF}
tabix {TABIX_OPTION} {OUTPUT_CHRX_FEMALE_GVCF}.gz
rm -f {OUTPUT_CHRX_FEMALE_GVCF}.idx

/usr/bin/java \\
  {HAPLOTYPE_JAVA_OPTION} \\
  -jar {GATK_JAR} HaplotypeCaller \\
  -I={INPUT_CRAM} \\
  -O={OUTPUT_CHRX_MALE_GVCF} \\
  -L={INTERVAL_CHRX} \\
  -R={REFERENCE} {HAPLOTYPE_OPTION_CHRX_MALE} 

bgzip {BGZIP_OPTION} {OUTPUT_CHRX_MALE_GVCF}
tabix {TABIX_OPTION} {OUTPUT_CHRX_MALE_GVCF}.gz
rm -f {OUTPUT_CHRX_MALE_GVCF}.idx

/usr/bin/java \\
  {HAPLOTYPE_JAVA_OPTION} \\
  -jar {GATK_JAR} HaplotypeCaller \\
  -I={INPUT_CRAM} \\
  -O={OUTPUT_CHRY_MALE_GVCF} \\
  -L={INTERVAL_CHRY} \\
  -R={REFERENCE} {HAPLOTYPE_OPTION_CHRY_MALE} 

bgzip {BGZIP_OPTION} {OUTPUT_CHRY_MALE_GVCF}
tabix {TABIX_OPTION} {OUTPUT_CHRY_MALE_GVCF}.gz
rm -f {OUTPUT_CHRY_MALE_GVCF}.idx

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

mkdir -p $(dirname {OUTPUT_AUTOSOME_GVCF})
{PBRUN} haplotypecaller \\
  --ref {REFERENCE} \\
  --in-bam {INPUT_CRAM} \\
  --out-variants {OUTPUT_AUTOSOME_GVCF} \\
  --interval-file {INTERVAL_AUTOSOME} {HAPLOTYPE_OPTION_AUTOSOME}

{BGZIP} {BGZIP_OPTION} {OUTPUT_AUTOSOME_GVCF}
{TABIX} {TABIX_OPTION} {OUTPUT_AUTOSOME_GVCF}.gz
rm -f {OUTPUT_AUTOSOME_GVCF}.idx

{PBRUN} haplotypecaller \\
  --ref {REFERENCE} \\
  --in-bam {INPUT_CRAM} \\
  --out-variants {OUTPUT_PAR_GVCF} \\
  --interval-file {INTERVAL_PAR} {HAPLOTYPE_OPTION_PAR}

{BGZIP} {BGZIP_OPTION} {OUTPUT_PAR_GVCF}
{TABIX} {TABIX_OPTION} {OUTPUT_PAR_GVCF}.gz
rm -f {OUTPUT_PAR_GVCF}.idx

{PBRUN} haplotypecaller \\
  --ref {REFERENCE} \\
  --in-bam {INPUT_CRAM} \\
  --out-variants {OUTPUT_CHRX_FEMALE_GVCF} \\
  --interval-file {INTERVAL_CHRX} {HAPLOTYPE_OPTION_CHRX_FEMALE}

{BGZIP} {BGZIP_OPTION} {OUTPUT_CHRX_FEMALE_GVCF}
{TABIX} {TABIX_OPTION} {OUTPUT_CHRX_FEMALE_GVCF}.gz
rm -f {OUTPUT_CHRX_FEMALE_GVCF}.idx

{PBRUN} haplotypecaller \\
  --ref {REFERENCE} \\
  --in-bam {INPUT_CRAM} \\
  --out-variants {OUTPUT_CHRX_MALE_GVCF} \\
  --interval-file {INTERVAL_CHRX} {HAPLOTYPE_OPTION_CHRX_MALE}

{BGZIP} {BGZIP_OPTION} {OUTPUT_CHRX_MALE_GVCF}
{TABIX} {TABIX_OPTION} {OUTPUT_CHRX_MALE_GVCF}.gz
rm -f {OUTPUT_CHRX_MALE_GVCF}.idx

{PBRUN} haplotypecaller \\
  --ref {REFERENCE} \\
  --in-bam {INPUT_CRAM} \\
  --out-variants {OUTPUT_CHRY_MALE_GVCF} \\
  --interval-file {INTERVAL_CHRY} {HAPLOTYPE_OPTION_CHRY_MALE}

{BGZIP} {BGZIP_OPTION} {OUTPUT_CHRY_MALE_GVCF}
{TABIX} {TABIX_OPTION} {OUTPUT_CHRY_MALE_GVCF}.gz
rm -f {OUTPUT_CHRY_MALE_GVCF}.idx
"""

STAGE_NAME = "haplotypecaller_parabricks"

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
    for sample in sample_conf.haplotype_call:
        output_dir = "%s/haplotypecaller/%s" % (run_conf.project_root, sample) 
        output_autosome_gvcf = "%s/%s.autosome.g.vcf" % (output_dir, sample)
        output_PAR_gvcf = "%s/%s.PAR.g.vcf" % (output_dir, sample)
        output_chrX_female_gvcf = "%s/%s.chrX.female.g.vcf" % (output_dir, sample)
        output_chrX_male_gvcf = "%s/%s.chrX.male.g.vcf" % (output_dir, sample)
        output_chrY_male_gvcf = "%s/%s.chrY.male.g.vcf" % (output_dir, sample)
        output_files[sample] = []
        output_files[sample].append(output_autosome_gvcf+".gz")
        output_files[sample].append(output_autosome_gvcf+".gz.tbi")
        output_files[sample].append(output_PAR_gvcf+".gz")
        output_files[sample].append(output_PAR_gvcf+".gz.tbi")
        output_files[sample].append(output_chrX_female_gvcf+".gz")
        output_files[sample].append(output_chrX_female_gvcf+".gz.tbi")
        output_files[sample].append(output_chrX_male_gvcf+".gz")
        output_files[sample].append(output_chrX_male_gvcf+".gz.tbi")
        output_files[sample].append(output_chrY_male_gvcf+".gz")
        output_files[sample].append(output_chrY_male_gvcf+".gz.tbi")
        haplotype_option = " " + gcat_conf.get(CONF_SECTION, "haplotype_threads_option") + " " + gcat_conf.get(CONF_SECTION, "haplotype_option")
        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_AUTOSOME_GVCF": output_autosome_gvcf,
            "OUTPUT_PAR_GVCF": output_PAR_gvcf,
            "OUTPUT_CHRX_FEMALE_GVCF": output_chrX_female_gvcf,
            "OUTPUT_CHRX_MALE_GVCF": output_chrX_male_gvcf,
            "OUTPUT_CHRY_MALE_GVCF": output_chrY_male_gvcf,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "GATK_JAR": gcat_conf.get(CONF_SECTION, "gatk_jar"),
            "INTERVAL_AUTOSOME": gcat_conf.path_get(CONF_SECTION, "interval_autosome"),
            "INTERVAL_PAR": gcat_conf.path_get(CONF_SECTION, "interval_par"),
            "INTERVAL_CHRX": gcat_conf.path_get(CONF_SECTION, "interval_chrx"),
            "INTERVAL_CHRY": gcat_conf.path_get(CONF_SECTION, "interval_chry"),
            "BGZIP_OPTION": gcat_conf.get(CONF_SECTION, "bgzip_option") + " " + gcat_conf.get(CONF_SECTION, "bgzip_threads_option"),
            "TABIX_OPTION": gcat_conf.get(CONF_SECTION, "tabix_option"),
            "HAPLOTYPE_OPTION_AUTOSOME": gcat_conf.get(CONF_SECTION, "haplotype_option_autosome") + haplotype_option,
            "HAPLOTYPE_OPTION_PAR": gcat_conf.get(CONF_SECTION, "haplotype_option_par") + haplotype_option,
            "HAPLOTYPE_OPTION_CHRX_FEMALE": gcat_conf.get(CONF_SECTION, "haplotype_option_chrx_female") + haplotype_option,
            "HAPLOTYPE_OPTION_CHRX_MALE": gcat_conf.get(CONF_SECTION, "haplotype_option_chrx_male") + haplotype_option,
            "HAPLOTYPE_OPTION_CHRY_MALE": gcat_conf.get(CONF_SECTION, "haplotype_option_chry_male") + haplotype_option,
            "HAPLOTYPE_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "haplotype_java_option")
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)
    
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
        output_autosome_gvcf = "%s/%s.autosome.g.vcf" % (output_dir, sample)
        output_PAR_gvcf = "%s/%s.PAR.g.vcf" % (output_dir, sample)
        output_chrX_female_gvcf = "%s/%s.chrX.female.g.vcf" % (output_dir, sample)
        output_chrX_male_gvcf = "%s/%s.chrX.male.g.vcf" % (output_dir, sample)
        output_chrY_male_gvcf = "%s/%s.chrY.male.g.vcf" % (output_dir, sample)
        output_files[sample] = []
        output_files[sample].append(output_autosome_gvcf+".gz")
        output_files[sample].append(output_autosome_gvcf+".gz.tbi")
        output_files[sample].append(output_PAR_gvcf+".gz")
        output_files[sample].append(output_PAR_gvcf+".gz.tbi")
        output_files[sample].append(output_chrX_female_gvcf+".gz")
        output_files[sample].append(output_chrX_female_gvcf+".gz.tbi")
        output_files[sample].append(output_chrX_male_gvcf+".gz")
        output_files[sample].append(output_chrX_male_gvcf+".gz.tbi")
        output_files[sample].append(output_chrY_male_gvcf+".gz")
        output_files[sample].append(output_chrY_male_gvcf+".gz.tbi")

        input_real_path = ""
        if not os.path.islink(input_bams[sample]):
            input_real_path = input_bams[sample]
        else:
            for path in sample_conf.bam_import_src[sample]:
                if not os.path.islink(path) and path.split(".")[-1] in ["bam", "cram"]:
                    input_real_path = path

        haplotype_option = " " + gcat_conf.get(CONF_SECTION, "haplotype_threads_option") + " " + gcat_conf.get(CONF_SECTION, "haplotype_option")
        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_real_path,
            "OUTPUT_AUTOSOME_GVCF": output_autosome_gvcf,
            "OUTPUT_PAR_GVCF": output_PAR_gvcf,
            "OUTPUT_CHRX_FEMALE_GVCF": output_chrX_female_gvcf,
            "OUTPUT_CHRX_MALE_GVCF": output_chrX_male_gvcf,
            "OUTPUT_CHRY_MALE_GVCF": output_chrY_male_gvcf,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "INTERVAL_AUTOSOME": gcat_conf.path_get(CONF_SECTION, "interval_autosome"),
            "INTERVAL_PAR": gcat_conf.path_get(CONF_SECTION, "interval_par"),
            "INTERVAL_CHRX": gcat_conf.path_get(CONF_SECTION, "interval_chrx"),
            "INTERVAL_CHRY": gcat_conf.path_get(CONF_SECTION, "interval_chry"),
            "BGZIP": gcat_conf.path_get(CONF_SECTION, "bgzip"),
            "BGZIP_OPTION": gcat_conf.get(CONF_SECTION, "bgzip_option") + " " + gcat_conf.get(CONF_SECTION, "bgzip_threads_option"),
            "TABIX": gcat_conf.path_get(CONF_SECTION, "tabix"),
            "TABIX_OPTION": gcat_conf.get(CONF_SECTION, "tabix_option"),
            "HAPLOTYPE_OPTION_AUTOSOME": gcat_conf.get(CONF_SECTION, "haplotype_option_autosome") + haplotype_option,
            "HAPLOTYPE_OPTION_PAR": gcat_conf.get(CONF_SECTION, "haplotype_option_par") + haplotype_option,
            "HAPLOTYPE_OPTION_CHRX_FEMALE": gcat_conf.get(CONF_SECTION, "haplotype_option_chrx_female") + haplotype_option,
            "HAPLOTYPE_OPTION_CHRX_MALE": gcat_conf.get(CONF_SECTION, "haplotype_option_chrx_male") + haplotype_option,
            "HAPLOTYPE_OPTION_CHRY_MALE": gcat_conf.get(CONF_SECTION, "haplotype_option_chry_male") + haplotype_option,
            "PBRUN": gcat_conf.get(CONF_SECTION, "pbrun"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)
    
    return output_files

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    if gcat_conf.safe_get(STAGE_NAME, "gpu_support", "False").lower() == "true":
        return _parabricks(input_bams, gcat_conf, run_conf, sample_conf)
    return _compatible(input_bams, gcat_conf, run_conf, sample_conf)

