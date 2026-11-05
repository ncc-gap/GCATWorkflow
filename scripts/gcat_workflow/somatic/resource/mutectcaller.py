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

if [ "_{USE_RECAL_FILE}" = "_T" ]; then
  /usr/bin/java \\
    {MUTECT_JAVA_OPTION} \\
    -jar {GATK_JAR} ApplyBQSR \\
    -R {REFERENCE} \\
    -I {INPUT_TUMOR_CRAM} \\
    --bqsr-recal-file {INPUT_TUMOR_RECAL} \\
    -O {OUTPUT_TUMOR_RECAL_CRAM}
  
  if [ "_{NORMAL_SAMPLE}" != "_" ]; then
    /usr/bin/java \\
      {MUTECT_JAVA_OPTION} \\
      -jar {GATK_JAR} ApplyBQSR \\
      -R {REFERENCE} \\
      -I {INPUT_NORMAL_CRAM} \\
      --bqsr-recal-file {INPUT_NORMAL_RECAL} \\
      -O {OUTPUT_NORMAL_RECAL_CRAM}
  fi
fi

/usr/bin/java \\
  {MUTECT_JAVA_OPTION} \\
  -jar {GATK_JAR} Mutect2 \\
  -I {MUTECT2_TUMOR_CRAM} -tumor {TUMOR_SAMPLE} {MUTECT2_NORMAL_CRAM} {NORMAL_SAMPLE} \\
  -O {OUTPUT_VCF} \\
  -R {REFERENCE} {MUTECT_OPTION}

bgzip {BGZIP_OPTION} {OUTPUT_VCF}
tabix {TABIX_OPTION} {OUTPUT_VCF}.gz
rm -f {OUTPUT_VCF}.idx

if [ "_{USE_RECAL_FILE}" = "_T" ]; then
  rm -f {OUTPUT_TUMOR_RECAL_CRAM}
  rm -f {OUTPUT_TUMOR_RECAL_CRAM}.crai
  rm -f {OUTPUT_NORMAL_RECAL_CRAM}
  rm -f {OUTPUT_NORMAL_RECAL_CRAM}.crai
fi
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

{BGZIP} {BGZIP_OPTION} {OUTPUT_VCF}
{TABIX} {TABIX_OPTION} {OUTPUT_VCF}.gz
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
        output_dir = "%s/mutectcaller/%s" % (run_conf.project_root, tumor) 
        os.makedirs(output_dir, exist_ok = True)
        
        output_vcf = "%s/%s.mutectcaller.vcf" % (output_dir, tumor)
        output_files[tumor] = [
            output_vcf + ".gz", output_vcf + ".gz.tbi"
        ]

        mutect_option = gcat_conf.get(CONF_SECTION, "mutect_option") + " " + gcat_conf.get(CONF_SECTION, "mutect_threads_option")
        
        use_recal_file = gcat_conf.safe_get(CONF_SECTION, "use_recal_file", "False").lower() == "true"
        recal_postfix = gcat_conf.get(CONF_SECTION, "recal_postfix")
        input_tumor_cram = input_bams[tumor]
        if use_recal_file:
            input_tumor_recal = "%s/%s%s" % (os.path.dirname(input_tumor_cram), tumor, recal_postfix)
            output_tumor_recal_cram = "%s/%s.recal.cram" % (output_dir, tumor)
            mutect2_tumor_cram = output_tumor_recal_cram
        else:
            input_tumor_recal = ""
            output_tumor_recal_cram = ""
            mutect2_tumor_cram = input_tumor_cram

        normal_sample = ""
        input_normal_cram = ""
        output_normal_recal_cram = ""
        mutect2_normal_cram = ""
        input_normal_recal = ""
        if normal != None:
            input_normal_cram = input_bams[normal]
            normal_sample = "-normal %s" % (normal)
            if use_recal_file:
                input_normal_recal = "%s/%s%s" % (os.path.dirname(input_normal_cram), normal, recal_postfix)
                output_normal_recal_cram = "%s/%s.recal.cram" % (output_dir, normal)
                mutect2_normal_cram = "-I %s" % (output_normal_recal_cram)
            else:
                mutect2_normal_cram = "-I %s" % (input_normal_cram)

        if gcat_conf.safe_get(CONF_SECTION, "out_f1r2_file", "False").lower() == "true":
            mutect_option += " --f1r2-tar-gz %s/%s.f1r2.tar.gz" % (output_dir, tumor)

        interval_files = gcat_conf.safe_get(CONF_SECTION, "interval_file", "")
        if interval_files != "":
            for interval_file in interval_files.split(","):
                if gcat_conf.exist_check and not os.path.exists(interval_file):
                    err_msg = "[%s] %s: path %s is not exists" % (section, option, interval_file)
                    raise ValueError(err_msg)
                mutect_option += " --intervals %s" % (interval_file)

        arguments = {
            "TUMOR_SAMPLE": tumor,
            "NORMAL_SAMPLE": normal_sample,
            "INPUT_TUMOR_CRAM": input_tumor_cram,
            "INPUT_NORMAL_CRAM": input_normal_cram,
            "INPUT_TUMOR_RECAL": input_tumor_recal,
            "INPUT_NORMAL_RECAL": input_normal_recal,
            "MUTECT2_TUMOR_CRAM": mutect2_tumor_cram,
            "MUTECT2_NORMAL_CRAM": mutect2_normal_cram,
            "USE_RECAL_FILE": "T" if use_recal_file else "F",
            "OUTPUT_TUMOR_RECAL_CRAM": output_tumor_recal_cram,
            "OUTPUT_NORMAL_RECAL_CRAM": output_normal_recal_cram,
            "OUTPUT_VCF":  output_vcf,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "GATK_JAR": gcat_conf.get(CONF_SECTION, "gatk_jar"),
            "MUTECT_OPTION": mutect_option,
            "MUTECT_JAVA_OPTION": gcat_conf.get(CONF_SECTION, "mutect_java_option"),
            "BGZIP_OPTION": gcat_conf.get(CONF_SECTION, "bgzip_option") + " " + gcat_conf.get(CONF_SECTION, "bgzip_threads_option"),
            "TABIX_OPTION": gcat_conf.get(CONF_SECTION, "tabix_option"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if tumor in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[tumor]
        if normal in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[normal]

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
        
        mutect_option = gcat_conf.get(CONF_SECTION, "mutect_option") + " " + gcat_conf.get(CONF_SECTION, "mutect_threads_option")

        input_real_path = ""
        if not os.path.islink(input_bams[tumor]):
            input_real_path = input_bams[tumor]
        else:
            for path in sample_conf.bam_import_src[tumor]:
                if not os.path.islink(path) and path.split(".")[-1] in ["bam", "cram"]:
                    input_real_path = path

        use_recal_file = gcat_conf.safe_get(CONF_SECTION, "use_recal_file", "False").lower() == "true"
        recal_postfix = gcat_conf.get(CONF_SECTION, "recal_postfix")
        if use_recal_file and input_real_path != "":
            mutect_option += " --in-tumor-recal-file %s/%s%s" % (os.path.dirname(input_real_path), tumor, recal_postfix)

        input_normal_cram = ""
        if normal != None:
            input_real_path_normal = ""
            if not os.path.islink(input_bams[normal]):
                input_real_path_normal = input_bams[normal]
            else:
                for path in sample_conf.bam_import_src[normal]:
                    if not os.path.islink(path) and path.split(".")[-1] in ["bam", "cram"]:
                        input_real_path_normal = path
            input_normal_cram = "--in-normal-bam %s --normal-name %s" % (input_real_path_normal, normal)
            if use_recal_file:
                mutect_option += " --in-normal-recal-file %s/%s%s" % (os.path.dirname(input_real_path_normal), normal, recal_postfix)

        if gcat_conf.safe_get(CONF_SECTION, "out_f1r2_file", "False").lower() == "true":
            mutect_option += " --mutect-f1r2-tar-gz %s/%s.f1r2.tar.gz" % (output_dir, tumor)

        interval_files = gcat_conf.safe_get(CONF_SECTION, "interval_file", "")
        if interval_files != "":
            for interval_file in interval_files.split(","):
                if gcat_conf.exist_check and not os.path.exists(interval_file):
                    err_msg = "[%s] %s: path %s is not exists" % (section, option, interval_file)
                    raise ValueError(err_msg)
                mutect_option += " --interval-file %s" % (interval_file)

        arguments = {
            "SAMPLE": tumor,
            "INPUT_TUMOR_CRAM": input_real_path,
            "INPUT_NORMAL_CRAM": input_normal_cram,
            "OUTPUT_VCF": output_vcf,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "MUTECT_OPTION": mutect_option,
            "PBRUN": gcat_conf.get(CONF_SECTION, "pbrun"),
            "BGZIP": gcat_conf.path_get(CONF_SECTION, "bgzip"),
            "BGZIP_OPTION": gcat_conf.get(CONF_SECTION, "bgzip_option") + " " + gcat_conf.get(CONF_SECTION, "bgzip_threads_option"),
            "TABIX": gcat_conf.path_get(CONF_SECTION, "tabix"),
            "TABIX_OPTION": gcat_conf.get(CONF_SECTION, "tabix_option"),
        }
        
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if tumor in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[tumor]
        if normal in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[normal]
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = tumor)
    
    return output_files

def configure(input_bams, gcat_conf, run_conf, sample_conf):

    if len(sample_conf.mutect_call) == 0:
        return {}

    if gcat_conf.safe_get(STAGE_NAME, "gpu_support", "False").lower() == "true":
        return _parabricks(input_bams, gcat_conf, run_conf, sample_conf)
    return _compatible(input_bams, gcat_conf, run_conf, sample_conf)

