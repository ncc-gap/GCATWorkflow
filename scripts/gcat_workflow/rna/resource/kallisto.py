#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

SCRIPT_HEADER = """#!/bin/bash
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
"""

class Kallisto(stage_task.Stage_task):
    def __init__(self, params):
        super().__init__(params)
        self.shell_script_template = """{bam_tofastq_command}
{cat_fastq_command}

/tools/kallisto-0.46.2/bin/kallisto quant \
    -i {REFERENCE_KALLISTO_INDEX} \
    --fusion \
    -o {OUTPUT_DIR} \
    {INPUT_FASTQ_1} {INPUT_FASTQ_2}

/tools/pizzly-0.37.3/bin/pizzly \
    --gtf {ANNOTATION_GTF} \
    --fasta {REFERENCE_FASTA} \
    --output {OUTPUT_DIR}/{SAMPLE_NAME}.pizzly \
    {PIZZLY_OPTION} \
    {OUTPUT_DIR}/fusion.txt

python /tools/pizzly-0.37.3/scripts/flatten_json.py \
    {OUTPUT_DIR}/{SAMPLE_NAME}.pizzly.unfiltered.json \
    > {OUTPUT_DIR}/{SAMPLE_NAME}.pizzly.unfiltered.table

python /tools/pizzly-0.37.3/scripts/flatten_json.py \
    {OUTPUT_DIR}/{SAMPLE_NAME}.pizzly.json \
    > {OUTPUT_DIR}/{SAMPLE_NAME}.pizzly.table

{remove_command}
"""

def configure(gcat_conf, run_conf, sample_conf):
    import os
    import gcat_workflow.rna.resource.bamtofastq as rs_bamtofastq
    
    STAGE_NAME = "kallisto"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": gcat_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Kallisto(params)
    bamtofastq_class = rs_bamtofastq.Bam_tofastq(params)
    
    output_files = {}
    for sample in sample_conf.kallisto:
        output_dir = "%s/kallisto/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)    
        output_files[sample] = "%s/%s.pizzly.table" % (output_dir, sample)
            
        f1_name = ""
        f2_name = ""
        bam_tofastq_command = SCRIPT_HEADER
        cat_fastq_command = ""
        remove_command = ""
        if sample in sample_conf.cram_import:
            raise Exception("cram to fastq is not yet supported.")
        elif sample in sample_conf.bam_import:
            f1_name = "%s/1_1_temp.fastq" % (output_dir)
            f2_name = "%s/2_1_temp.fastq" % (output_dir)
            
            bam_tofastq_command = bamtofastq_class.shell_script_template.format(
                param = gcat_conf.get("bam_tofastq", "params"),
                input_bam = sample_conf.bam_import[sample],
                f1_name = f1_name,
                f2_name = f2_name,
                o1_name = output_dir + '/unmatched_first_output.txt',
                o2_name = output_dir + '/unmatched_second_output.txt',
                t = output_dir + '/temp.txt',
                s = output_dir + '/single_end_output.txt',
                pass_file = output_dir + '/pass.txt'
            )
            remove_command = "rm %s %s %s" % (output_dir + '/pass.txt', f1_name, f2_name)
            
        else:
            paired = len(sample_conf.fastq[sample]) > 1
            if len(sample_conf.fastq[sample][0]) == 1:
                f1_name = sample_conf.fastq[sample][0][0]
                if paired:
                    f2_name = sample_conf.fastq[sample][1][0]
            else:
                cat_fastq_command += "cat {FASTQ1s} > {OUTPUT_DIR}/1_1_temp.fastq\n".format(
                    FASTQ1s = " ".join(sample_conf.fastq[sample][0]),
                    OUTPUT_DIR = output_dir
                )
                f1_name = "%s/1_1_temp.fastq" % (output_dir)
                if paired:
                    cat_fastq_command += "cat {FASTQ2s} > {OUTPUT_DIR}/2_1_temp.fastq\n".format(
                        FASTQ2s = " ".join(sample_conf.fastq[sample][1]),
                        OUTPUT_DIR = output_dir
                    )
                    f2_name = "%s/2_1_temp.fastq" % (output_dir)
                remove_command = "rm %s %s" % (f1_name, f2_name)
                
        arguments = {
            "SAMPLE": sample,
            "OUTPUT_DIR": output_dir,
            "bam_tofastq_command": bam_tofastq_command,
            "cat_fastq_command": cat_fastq_command,
            "remove_command": remove_command,
            "SAMPLE_NAME": sample,
            "INPUT_FASTQ_1": f1_name,
            "INPUT_FASTQ_2": f2_name,
            "REFERENCE_FASTA": gcat_conf.path_get(SECTION_NAME, "reference_fasta"),
            "REFERENCE_KALLISTO_INDEX": gcat_conf.path_get(SECTION_NAME, "reference_kallisto_index"),
            "ANNOTATION_GTF": gcat_conf.path_get(SECTION_NAME, "annotation_gtf"),
            "PIZZLY_OPTION": gcat_conf.get(SECTION_NAME, "pizzly_option"),
            "OUTPUT_DIR": output_dir
        }
        
        singularity_bind = [
            run_conf.project_root,
            os.path.dirname(gcat_conf.path_get(SECTION_NAME, "reference_fasta")),
            os.path.dirname(gcat_conf.path_get(SECTION_NAME, "reference_kallisto_index")),
            os.path.dirname(gcat_conf.path_get(SECTION_NAME, "annotation_gtf"))
        ]
        if sample in sample_conf.fastq_src:
            singularity_bind += sample_conf.fastq_src[sample][0] + sample_conf.fastq_src[sample][1]

        elif sample in sample_conf.bam_import:
            singularity_bind += sample_conf.bam_import_src[sample]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)

    return output_files
