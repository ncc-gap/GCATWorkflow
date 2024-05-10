#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

OUTPUT_FQ1_FORMAT = "sra_fastq_dump/{sample}/1_1.fasq"
OUTPUT_FQ2_FORMAT = "sra_fastq_dump/{sample}/1_2.fasq"

class SRA_fastq_dump(stage_task.Stage_task):
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

mkdir -p {OUTPUT_DIR}/temp
cd {OUTPUT_DIR}/temp

if [ "{SRA_PATH}" = "" ]; then
    # download
    prefetch --max-size 100000000 {PREFETCH_OPTION} {RUN_ID}
    fasterq-dump -v --split-files {DUMP_OPTION} {RUN_ID}/{RUN_ID}.sra
    rm -rf {RUN_ID}

    if [ -e {RUN_ID}.sra_1.fastq ]; then mv {RUN_ID}.sra_1.fastq 1_1.fastq; fi
    if [ -e {RUN_ID}.sra_2.fastq ]; then mv {RUN_ID}.sra_2.fastq 1_2.fastq; fi
    if [ -e {RUN_ID}.sra.fastq ]; then mv {RUN_ID}.sra.fastq 1_1.fastq; fi
    if [ -e {RUN_ID}_1.fastq ]; then mv {RUN_ID}_1.fastq 1_1.fastq; fi
    if [ -e {RUN_ID}_2.fastq ]; then mv {RUN_ID}_2.fastq 1_2.fastq; fi
    if [ -e {RUN_ID}.fastq ]; then mv {RUN_ID}.fastq 1_1.fastq; fi

else
    # for ddbj
    if [ "{WGET}" = "T" ]; then
        wget {SRA_PATH}
        fasterq-dump -v --split-files {DUMP_OPTION} {RUN_ID}.sra
        rm {RUN_ID}.sra
    else
        fasterq-dump -v --split-files {DUMP_OPTION} {SRA_PATH}
    fi
    if [ -e {RUN_ID}_1.fastq ]; then mv {RUN_ID}_1.fastq 1_1.fastq; fi
    if [ -e {RUN_ID}_2.fastq ]; then mv {RUN_ID}_2.fastq 1_2.fastq; fi
    if [ -e {RUN_ID}.fastq ]; then mv {RUN_ID}.fastq 1_1.fastq; fi
fi

# squeeze
seqtk squeeze {OUTPUT_DIR}/temp/1_1.fastq > {OUTPUT_DIR}/1_1.fastq

if [ -e {OUTPUT_DIR}/temp/1_2.fastq ]
then
    seqtk squeeze {OUTPUT_DIR}/temp/1_2.fastq > {OUTPUT_DIR}/1_2.fastq
fi

cd {OUTPUT_DIR}
rm -rf {OUTPUT_DIR}/temp

touch {PASS_FILE}
"""

def configure(gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "sra_fastq_dump"
    CONF_SECTION = "sra_fastq_dump"
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = SRA_fastq_dump(params)
    output_fastqs = {}
    for sample in sample_conf.sra_fastq_dump:
        output_dir = "%s/fastq/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)    
        f1_name = output_dir + "/1_1.fastq"
        f2_name = output_dir + "/1_2.fastq"
        output_fastqs[sample] = [[f1_name], [f2_name]]

        sra_path = sample_conf.sra_fastq_dump[sample][1]
        download = "F"
        singularity_bind = [
            run_conf.project_root,
        ]
        if sra_path == None:
            # prefetch from SRA
            sra_path = ""
        elif ":" in sra_path:
            # wget from ddbj
            download = "T"
        else:
            # direct access in ddbj
            singularity_bind.append(os.path.dirname(sra_path))
	
        arguments = {
            "RUN_ID": sample_conf.sra_fastq_dump[sample][0],
            "OUTPUT_DIR": output_dir,
            "PASS_FILE": output_dir + '/pass.txt',
            "SRA_PATH": sra_path,
            "WGET": download,
            "PREFETCH_OPTION": gcat_conf.get(CONF_SECTION, "prefetch_option"),
            "DUMP_OPTION": gcat_conf.get(CONF_SECTION, "dump_option"),
        }
        
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)
    return output_fastqs
