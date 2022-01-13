#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class Genomon_mutation_call(stage_task.Stage_task):
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

SAMTOOLS=samtools

# interval_list
interval_option=""
if [ "_{FISHER_INTERVAL_LIST}" != "_" ]; then
  interval_option="-L {FISHER_INTERVAL_LIST}"
fi

# print_header
mut_header=""

ext=""
if [ _{OUTPUT_FORMAT} = "_vcf" ]; then
   ext=vcf
fi

if [ _{SAMPLE2} = "_None" ]; then

  # Fisher's Exact Test
  time fisher single -o {OUTPUT_PREF}.fisher_mutations.${{ext}} -O {OUTPUT_FORMAT} --ref_fa {REFERENCE} -a {SAMPLE1} -1 {INPUT_BAM1} --samtools_path ${{SAMTOOLS}} {FISHER_SINGLE_OPTION} ${{interval_option}} --samtools_params "{FISHER_SAMTOOLS_OPTION}"

  # Local realignment using edlib. The candidate mutations are varidated.
  time mutfilter realignment --target_mutation_file {OUTPUT_PREF}.fisher_mutations.${{ext}} -O {OUTPUT_FORMAT} -A {SAMPLE1} -1 {INPUT_BAM1} --output {OUTPUT_PREF}.realignment_mutations.${{ext}} --ref_genome {REFERENCE}  {REALIGNMENT_OPTION}

  # Annotation if the candidate is on the simplerepeat. 
  time mutfilter simplerepeat --target_mutation_file {OUTPUT_PREF}.realignment_mutations.${{ext}} -O {OUTPUT_FORMAT} --output {OUTPUT_PREF}.simplerepeat_mutations.${{ext}} --simple_repeat_db {ANNOTATION_DB}/simpleRepeat.bed.gz {SIMPLE_REPEAT_OPTION}

  # header
  if [ _{OUTPUT_FORMAT} != "_vcf" ]; then
    mut_header="Chr,Start,End,Ref,Alt,depth,variantNum,bases,A_C_G_T,misRate,strandRatio,10%_posterior_quantile,posterior_mean,90%_posterior_quantile,readPairNum,variantPairNum,otherPairNum,10%_posterior_quantile(realignment),posterior_mean(realignment),90%_posterior_quantile(realignment),simple_repeat_pos,simple_repeat_seq"
    print_header=`echo ${{mut_header}} | tr "," "\t"`
    echo "${{print_header}}" > {OUTPUT_PREF}.genomon_mutation.result.${{ext}}
    cat {OUTPUT_PREF}.simplerepeat_mutations.${{ext}} >> {OUTPUT_PREF}.genomon_mutation.result.${{ext}}
  else
    cp {OUTPUT_PREF}.simplerepeat_mutations.${{ext}} {OUTPUT_PREF}.genomon_mutation.result.${{ext}}
  fi

  mutil filter -i {OUTPUT_PREF}.genomon_mutation.result.${{ext}} -O {OUTPUT_FORMAT} -1 {SAMPLE1} -o {OUTPUT_PREF}.genomon_mutation.result.filt.${{ext}} {FILTER_SINGLE_OPTION}

else
  # Fisher's Exact Test
  time fisher comparison -o {OUTPUT_PREF}.fisher_mutations.${{ext}} -O {OUTPUT_FORMAT} --ref_fa {REFERENCE} -a {SAMPLE1} -b {SAMPLE2} -1 {INPUT_BAM1} -2 {INPUT_BAM2} --samtools_path ${{SAMTOOLS}} {FISHER_PAIR_OPTION} ${{interval_option}} --samtools_params "{FISHER_SAMTOOLS_OPTION}"

  # Local realignment using blat. The candidate mutations are varidated.
  time mutfilter realignment --target_mutation_file {OUTPUT_PREF}.fisher_mutations.${{ext}} -O {OUTPUT_FORMAT} -A {SAMPLE1} -B {SAMPLE2} -1 {INPUT_BAM1} -2 {INPUT_BAM2} --output {OUTPUT_PREF}.realignment_mutations.${{ext}} --ref_genome {REFERENCE} {REALIGNMENT_OPTION}

  # Annotation if the candidate is near Indel. 
  time mutfilter indel --target_mutation_file {OUTPUT_PREF}.realignment_mutations.${{ext}} -O {OUTPUT_FORMAT} -A {SAMPLE1} -B {SAMPLE2} -2 {INPUT_BAM2} --output {OUTPUT_PREF}.indel_mutations.${{ext}} --ref_genome {REFERENCE} --samtools_path ${{SAMTOOLS}} {INDEL_OPTION} --samtools_params "{INDEL_SAMTOOLS_OPTION}"

  # Annotation if the candidate is near the breakpoint. 
  time mutfilter breakpoint --target_mutation_file {OUTPUT_PREF}.indel_mutations.${{ext}} -O {OUTPUT_FORMAT} -A {SAMPLE1} -B {SAMPLE2} -2 {INPUT_BAM2} --output {OUTPUT_PREF}.breakpoint_mutations.${{ext}} --ref_genome {REFERENCE} {BREAKPOINT_OPTION}

  # Annotation if the candidate is on the simplerepeat. 
  time mutfilter simplerepeat --target_mutation_file {OUTPUT_PREF}.breakpoint_mutations.${{ext}} -O {OUTPUT_FORMAT} --output {OUTPUT_PREF}.simplerepeat_mutations.${{ext}} --simple_repeat_db {ANNOTATION_DB}/simpleRepeat.bed.gz {SIMPLE_REPEAT_OPTION}

  # header
  if [ _{OUTPUT_FORMAT} != "_vcf" ]; then
    mut_header="Chr,Start,End,Ref,Alt,depth_tumor,variantNum_tumor,depth_normal,variantNum_normal,bases_tumor,bases_normal,A_C_G_T_tumor,A_C_G_T_normal,misRate_tumor,strandRatio_tumor,misRate_normal,strandRatio_normal,P-value(fisher),readPairNum_tumor,variantPairNum_tumor,otherPairNum_tumor,readPairNum_normal,variantPairNum_normal,otherPairNum_normal,P-value(fisher_realignment),indel_mismatch_count,indel_mismatch_rate,bp_mismatch_count,distance_from_breakpoint,simple_repeat_pos,simple_repeat_seq"
    print_header=`echo ${{mut_header}} | tr "," "\t"`
    echo "${{print_header}}" > {OUTPUT_PREF}.genomon_mutation.result.${{ext}}
    cat {OUTPUT_PREF}.simplerepeat_mutations.${{ext}} >> {OUTPUT_PREF}.genomon_mutation.result.${{ext}}
  else
    cp {OUTPUT_PREF}.simplerepeat_mutations.${{ext}} {OUTPUT_PREF}.genomon_mutation.result.${{ext}}
  fi

  mutil filter -i {OUTPUT_PREF}.genomon_mutation.result.${{ext}} -O {OUTPUT_FORMAT} -o {OUTPUT_PREF}.genomon_mutation.result.filt.${{ext}} {FILTER_PAIR_OPTION}
fi 

bgzip {BGZIP_OPTION} {OUTPUT_PREF}.genomon_mutation.result.${{ext}}
tabix {TABIX_OPTION} {OUTPUT_PREF}.genomon_mutation.result.${{ext}}.gz

bgzip {BGZIP_OPTION} {OUTPUT_PREF}.genomon_mutation.result.filt.${{ext}}
tabix {TABIX_OPTION} {OUTPUT_PREF}.genomon_mutation.result.filt.${{ext}}.gz

rm -f {OUTPUT_PREF}.fisher_mutations.${{ext}}
rm -f {OUTPUT_PREF}.realignment_mutations.${{ext}}
rm -f {OUTPUT_PREF}.indel_mutations.${{ext}}
rm -f {OUTPUT_PREF}.breakpoint_mutations.${{ext}}
rm -f {OUTPUT_PREF}.simplerepeat_mutations.${{ext}}

"""

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    
    STAGE_NAME = "genomon_mutation_call"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Genomon_mutation_call(params)
    
    output_files = {}
    for (tumor, normal) in sample_conf.genomon_mutation_call:
        output_dir = "%s/genomon_mutation/%s" % (run_conf.project_root, tumor)
        os.makedirs(output_dir, exist_ok=True)
        output_file = "%s/%s.genomon_mutation.result.vcf.gz" % (output_dir, tumor)
        output_files[tumor] = [
            output_file, output_file + ".tbi"
        ]
        input_normal_cram = ""
        if normal != None:
            input_normal_cram = input_bams[normal]

        arguments = {
            "INPUT_BAM1": input_bams[tumor],
            "INPUT_BAM2": input_normal_cram,
            "REFERENCE": gcat_conf.path_get(CONF_SECTION, "reference"),
            "FISHER_INTERVAL_LIST": gcat_conf.get(CONF_SECTION, "filter_interval_list"),
            "FISHER_PAIR_OPTION": gcat_conf.get(CONF_SECTION, "fisher_pair_option"),
            "FISHER_SINGLE_OPTION": gcat_conf.get(CONF_SECTION, "fisher_single_option"),
            "FISHER_SAMTOOLS_OPTION": gcat_conf.get(CONF_SECTION, "fisher_samtools_option"),
            "REALIGNMENT_OPTION": gcat_conf.get(CONF_SECTION, "realignment_option")+ " " + gcat_conf.get(CONF_SECTION, "realignment_threads_option"),
            "INDEL_OPTION": gcat_conf.get(CONF_SECTION, "indel_option") + " " + gcat_conf.get(CONF_SECTION, "indel_threads_option"),
            "INDEL_SAMTOOLS_OPTION": gcat_conf.get(CONF_SECTION, "indel_samtools_option"),
            "BREAKPOINT_OPTION": gcat_conf.get(CONF_SECTION, "breakpoint_option"),
            "SIMPLE_REPEAT_OPTION": gcat_conf.get(CONF_SECTION, "simple_repeat_option"),
            "OUTPUT_PREF": "%s/%s" % (output_dir, tumor),
            "SAMPLE1": tumor,
            "SAMPLE2": normal,
            "ANNOTATION_DB": gcat_conf.get(CONF_SECTION, "annotation_db"),
            "OUTPUT_FORMAT": gcat_conf.get(CONF_SECTION, "output_format"),
            "FILTER_PAIR_OPTION": gcat_conf.get(CONF_SECTION, "filter_pair_option"),
            "FILTER_SINGLE_OPTION": gcat_conf.get(CONF_SECTION, "filter_single_option"),
            "BGZIP_OPTION": gcat_conf.get(CONF_SECTION, "bgzip_option") + " " + gcat_conf.get(CONF_SECTION, "bgzip_threads_option"),
            "TABIX_OPTION": gcat_conf.get(CONF_SECTION, "tabix_option"),
        }
       
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if tumor in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[tumor]
        if normal in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[normal]
        if gcat_conf.get(CONF_SECTION, "filter_interval_list") != "":
            singularity_bind.append(gcat_conf.path_get(CONF_SECTION, "filter_interval_list"))
        singularity_bind.append(gcat_conf.path_get(CONF_SECTION, "annotation_db"))
        
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = tumor)
    
    return output_files
