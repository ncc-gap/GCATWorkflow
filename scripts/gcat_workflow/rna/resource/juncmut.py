#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

class Juncmut(stage_task.Stage_task):
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

rm -rf {OUTPUT_DIR}/*
{DECOMPRESS_CMD}
juncmut get \
-input_file {INPUT_FILE} \
-output_file {OUTPUT_DIR}/{OUTPUT_FILE} \
-output_bam {OUTPUT_DIR}/{OUTPUT_BAM} \
-reference {REFERENCE} \
-rna_bam {INPUT_BAM} \
-control_file {CONTROL_FILE1} {CONTROL_FILE2} \
-genecode_gene_file {GENE_FILE}  \
{JUNCMUT_PRAM} -gnomad {GNOMAD}
{RM_CMD}
"""

def configure(input_bams, input_sj_tabs, gcat_conf, run_conf, sample_conf):
    import os
    import urllib
    
    STAGE_NAME = "juncmut"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": gcat_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Juncmut(params)
    
    output_files = {}
    dbs = [
        (SECTION_NAME, "reference"),
        (SECTION_NAME, "control_file1"),
        (SECTION_NAME, "control_file2"),
        (SECTION_NAME, "genecode_gene_file"),
        (SECTION_NAME, "gnomad")
    ]
    local_dbs = []
    for (section, db_name) in dbs:
        parsed = urllib.parse.urlparse(gcat_conf.get(section, db_name))
        if parsed.scheme == "":
            local_dbs.append(gcat_conf.path_get(section, db_name))
            
    for sample in sample_conf.juncmut:
        output_dir = "%s/juncmut/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)
        output_files[sample] = [
            "%s/%s.juncmut.filt.bam" % (output_dir, sample),
            "%s/%s.juncmut.filt.bam.bai" % (output_dir, sample)
        ]
        decomp = ""
        remove = ""
        sjtab = input_sj_tabs[sample]
        if input_sj_tabs[sample].endswith(".gz"):
            sjtab_comp = "%s/%s" % (output_dir, os.path.basename(input_sj_tabs[sample]))
            (sjtab, ext) = os.path.splitext(sjtab_comp)
            decomp = "cp {input} {sjtab}\ngunzip {sjtab}".format(
                input = input_sj_tabs[sample],
                sjtab = sjtab_comp
            )
            remove = "rm %s" % (sjtab)
        arguments = {
            "SAMPLE": sample,
            "INPUT_FILE": sjtab,
            "INPUT_BAM": input_bams[sample],
            "OUTPUT_DIR": output_dir,
            "OUTPUT_FILE": "%s.juncmut.txt" % (sample),
            "OUTPUT_BAM": "%s.juncmut.filt.bam" % (sample),
            "REFERENCE": gcat_conf.path_get(SECTION_NAME, "reference"),
            "CONTROL_FILE1": gcat_conf.get(SECTION_NAME, "control_file1"),
            "CONTROL_FILE2": gcat_conf.get(SECTION_NAME, "control_file2"),
            "GENE_FILE": gcat_conf.get(SECTION_NAME, "genecode_gene_file"),
            "GNOMAD": gcat_conf.get(SECTION_NAME, "gnomad"),
            "JUNCMUT_PRAM": gcat_conf.get(SECTION_NAME, "juncmut_pram"),
            "DECOMPRESS_CMD": decomp,
            "RM_CMD": remove
        }
       
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
        
        for db in local_dbs:
            singularity_bind.append(os.path.dirname(db))
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)

    return output_files
