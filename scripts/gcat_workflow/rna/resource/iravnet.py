#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

class Iravnet(stage_task.Stage_task):
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

iravnet get {INPUT_BAM} {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf.tmp1 {REF}

bcftools norm {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf.tmp1 -f {REF} > {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf.tmp2

iravnet annotate {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf.tmp2 {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf --gnomad_exome {GNOMAD_EXOME} --gnomad_genome {GNOMAD_GENOME} --clinvar {CLINVAR_DB}
bgzip -c {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf > {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf.gz
tabix -f -p vcf {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf.gz
bcftools filter -e "INFO/GNOMAD_GENOME[0] > 0.01 | INFO/GNOMAD_EXOME[0] > 0.01 | IR_MT / (SJ_WT + IR_MT) < 0.1 | IR_MT_MOH < 25"  ${OUTPUT_DIR}/{SAMPLE}.iravnet.vcf.gz > ${OUTPUT_DIR}/{SAMPLE}.iravnet.filt.vcf

iravnet filt_bam {OUTPUT_DIR}/{SAMPLE}.iravnet.filt.vcf {INPUT_BAM} {OUTPUT_DIR}/{SAMPLE}.iravnet.filt.bam

samtools index {OUTPUT_DIR}/{SAMPLE}.iravnet.filt.bam

bgzip -f {OUTPUT_DIR}/{SAMPLE}.iravnet.filt.vcf
tabix -f -p vcf {OUTPUT_DIR}/{SAMPLE}.iravnet.filt.vcf.gz

rm -rf {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf.tmp1
rm -rf {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf.tmp2
rm -rf {OUTPUT_DIR}/{SAMPLE}.iravnet.vcf
"""

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    import os
    import urllib
    
    STAGE_NAME = "iravnet"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": gcat_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Iravnet(params)
    
    output_files = {}
    dbs = [
        (SECTION_NAME, "reference"),
        (SECTION_NAME, "clinvar_db"),
        (SECTION_NAME, "gnomad_exome"),
        (SECTION_NAME, "gnomad_genome")
    ]
    local_dbs = []
    for (section, db_name) in dbs:
        parsed = urllib.parse.urlparse(gcat_conf.get(section, db_name))
        if parsed.scheme == "":
            local_dbs.append(gcat_conf.path_get(section, db_name))
            
    for sample in sample_conf.iravnet:
        output_dir = "%s/iravnet/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)    
        output_files[sample] = [
            "%s/%s.iravnet.filt.bam" % (output_dir, sample),
            "%s/%s.iravnet.filt.bam.bai" % (output_dir, sample)
        ]
        arguments = {
            "SAMPLE": sample,
            "INPUT_BAM": input_bams[sample],
            "OUTPUT_DIR": output_dir,
            "REF": gcat_conf.path_get(SECTION_NAME, "reference"),
            "CLINVAR_DB": gcat_conf.path_get(SECTION_NAME, "clinvar_db"),
            "GNOMAD_EXOME": gcat_conf.get(SECTION_NAME, "gnomad_exome"),
            "GNOMAD_GENOME": gcat_conf.get(SECTION_NAME, "gnomad_genome")
        }
       
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
        
        for db in local_dbs:
            singularity_bind.append(os.path.dirname(db))
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)

    return output_files
