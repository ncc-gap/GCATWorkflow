#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task
import glob

class Glimpse(stage_task.Stage_task):
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

function glimpse_split() {{
    CHR_NUM=${{1}}
    MAP=${{2}}

    REFBCF={GLIMPSE_REFERENCE_PREFIX}${{CHR_NUM}}.bcf
    REFVCF={GLIMPSE_REFERENCE_PREFIX}${{CHR_NUM}}.vcf.gz
    REFTSV={GLIMPSE_REFERENCE_PREFIX}${{CHR_NUM}}.tsv.gz
    REFCHUNK={GLIMPSE_REFERENCE_PREFIX}${{CHR_NUM}}.chunks.txt

    mkdir -p {OUTPUT_DIR}/temp

    bcftools mpileup -f {REFERENCE} -I -E -a 'FORMAT/DP' -T ${{REFVCF}} -r ${{CHR_NUM}} {INPUT_CRAM} -Ou | bcftools call -Aim -C alleles -T ${{REFTSV}} -Oz -o {OUTPUT_DIR}/temp/{SAMPLE}_mpileup.vcf.gz
    bcftools index -f {OUTPUT_DIR}/temp/{SAMPLE}_mpileup.vcf.gz

    while IFS="" read -r LINE || [ -n "$LINE" ];
    do
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        {GLIMPSE_PHASE_PATH} {GLIMPSE_PHASE_OPTION} \\
            --input {OUTPUT_DIR}/temp/{SAMPLE}_mpileup.vcf.gz \\
            --reference ${{REFBCF}} \\
            --map ${{MAP}} \\
            --input-region ${{IRG}} \\
            --output-region ${{ORG}} \\
            --output {OUTPUT_DIR}/temp/{SAMPLE}.${{CHR_NUM}}.${{ID}}.bcf
        bcftools index -f {OUTPUT_DIR}/temp/{SAMPLE}.${{CHR_NUM}}.${{ID}}.bcf
    done < ${{REFCHUNK}}

    ls {OUTPUT_DIR}/temp/{SAMPLE}.${{CHR_NUM}}.*.bcf > {OUTPUT_DIR}/temp/{SAMPLE}.list.${{CHR_NUM}}.txt
    {GLIMPSE_LIGATE_PATH} --input {OUTPUT_DIR}/temp/{SAMPLE}.list.${{CHR_NUM}}.txt --output {OUTPUT_DIR}/temp/{SAMPLE}.${{CHR_NUM}}.merged.bcf
    bcftools index -f {OUTPUT_DIR}/temp/{SAMPLE}.${{CHR_NUM}}.merged.bcf

    {GLIMPSE_SAMPLE_PATH} --input {OUTPUT_DIR}/temp/{SAMPLE}.${{CHR_NUM}}.merged.bcf --solve --output {OUTPUT_DIR}/{SAMPLE}.${{CHR_NUM}}.phased.bcf
    bcftools index -f {OUTPUT_DIR}/{SAMPLE}.${{CHR_NUM}}.phased.bcf

    rm -rf {OUTPUT_DIR}/temp
}}

glimpse_split chr1  {MAP_chr1}
glimpse_split chr2  {MAP_chr2}
glimpse_split chr3  {MAP_chr3}
glimpse_split chr4  {MAP_chr4}
glimpse_split chr5  {MAP_chr5}
glimpse_split chr6  {MAP_chr6}
glimpse_split chr7  {MAP_chr7}
glimpse_split chr8  {MAP_chr8}
glimpse_split chr9  {MAP_chr9}
glimpse_split chr10 {MAP_chr10}
glimpse_split chr11 {MAP_chr11}
glimpse_split chr12 {MAP_chr12}
glimpse_split chr13 {MAP_chr13}
glimpse_split chr14 {MAP_chr14}
glimpse_split chr15 {MAP_chr15}
glimpse_split chr16 {MAP_chr16}
glimpse_split chr17 {MAP_chr17}
glimpse_split chr18 {MAP_chr18}
glimpse_split chr19 {MAP_chr19}
glimpse_split chr20 {MAP_chr20}
glimpse_split chr21 {MAP_chr21}
glimpse_split chr22 {MAP_chr22}

bcftools concat -Oz \\
 {OUTPUT_DIR}/{SAMPLE}.chr1.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr2.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr3.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr4.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr5.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr6.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr7.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr8.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr9.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr10.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr11.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr12.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr13.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr14.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr15.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr16.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr17.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr18.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr19.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr20.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr21.phased.bcf \\
 {OUTPUT_DIR}/{SAMPLE}.chr22.phased.bcf \\
 > {OUTPUT_DIR}/{SAMPLE}.allchr.concat.vcf.gz

tabix -p vcf -f {OUTPUT_DIR}/{SAMPLE}.allchr.concat.vcf.gz
"""

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    
    output_files = {}
    if len(sample_conf.glimpse) == 0:
        return output_files

    STAGE_NAME = "glimpse"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Glimpse(params)

    for sample in sample_conf.glimpse:
        output_dir = "%s/glimpse/%s" % (run_conf.project_root, sample)

        output_files[sample] = [
            "%s/%s.allchr.concat.vcf.gz" % (output_dir, sample),
            "%s/%s.allchr.concat.vcf.gz.tbi" % (output_dir, sample)
        ]

        map_dir = gcat_conf.path_get(CONF_SECTION, "map_dir")
        map_files = {}
        for chr_num in range(1,23):
            map_files[chr_num] = sorted(glob.glob("%s/chr%d.*.gz" % (map_dir, chr_num)))[0]

        arguments = {
            "SAMPLE": sample,
            "INPUT_CRAM": input_bams[sample],
            "OUTPUT_DIR":  output_dir,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "GLIMPSE_PHASE_PATH": gcat_conf.get(CONF_SECTION, "glimpse_phase_path"),
            "GLIMPSE_LIGATE_PATH": gcat_conf.get(CONF_SECTION, "glimpse_ligate_path"),
            "GLIMPSE_SAMPLE_PATH": gcat_conf.get(CONF_SECTION, "glimpse_sample_path"),
            "GLIMPSE_PHASE_OPTION": gcat_conf.get(CONF_SECTION, "glimpse_phase_option") + " " + gcat_conf.get(CONF_SECTION, "glimpse_phase_threads_option"),
            "GLIMPSE_REFERENCE_PREFIX":  gcat_conf.get(CONF_SECTION, "glimpse_reference_prefix"),
            "MAP_chr1": map_files[1],
            "MAP_chr2": map_files[2],
            "MAP_chr3": map_files[3],
            "MAP_chr4": map_files[4],
            "MAP_chr5": map_files[5],
            "MAP_chr6": map_files[6],
            "MAP_chr7": map_files[7],
            "MAP_chr8": map_files[8],
            "MAP_chr9": map_files[9],
            "MAP_chr10": map_files[10],
            "MAP_chr11": map_files[11],
            "MAP_chr12": map_files[12],
            "MAP_chr13": map_files[13],
            "MAP_chr14": map_files[14],
            "MAP_chr15": map_files[15],
            "MAP_chr16": map_files[16],
            "MAP_chr17": map_files[17],
            "MAP_chr18": map_files[18],
            "MAP_chr19": map_files[19],
            "MAP_chr20": map_files[20],
            "MAP_chr21": map_files[21],
            "MAP_chr22": map_files[22],
        }
        
        singularity_bind = [
            run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference")), 
            os.path.dirname(gcat_conf.get(CONF_SECTION, "glimpse_reference_prefix")),
            gcat_conf.path_get(CONF_SECTION, "map_dir"),
        ]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)
    
    return output_files
