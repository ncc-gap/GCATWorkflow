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

juncmut detect \
  {SJ_OUTTAB} \
  {BAM} \
  {OUTPUT_PREFIX}.juncmut.txt \
  {REFERENCE} \
  {GENCODE} \
  --control_file {CONTROL_FILE1} {CONTROL_FILE2} {JUNCMUT_DETECT_PARAM}

juncmut filt_bam \
  {OUTPUT_PREFIX}.juncmut.txt \
  {BAM} \
  {OUTPUT_PREFIX}.juncmut.filt.bam \
  {GENCODE}

juncmut sjclass \
  {OUTPUT_PREFIX}.juncmut.txt \
  {OUTPUT_PREFIX}.juncmut.sjclass.txt \
  {BAM} \
  {SJ_OUTTAB} \
  {REFERENCE} \
  {GENCODE} {JUNCMUT_SJCLASS_PARAM}

juncmut alu \
  {OUTPUT_PREFIX}.juncmut.sjclass.txt \
  {OUTPUT_PREFIX}.juncmut.sjclass.alu.txt \
  {RMSK_BED} \
  {REFERENCE} {JUNCMUT_ALU_PARAM}

juncmut annot \
  {OUTPUT_PREFIX}.juncmut.sjclass.alu.txt \
  {OUTPUT_PREFIX}.juncmut.sjclass.alu.annot.txt \
  {REFERENCE} {JUNCMUT_ANNOT_PARAM}

juncmut filt \
  {OUTPUT_PREFIX}.juncmut.sjclass.alu.annot.txt \
  {OUTPUT_PREFIX}.juncmut.sjclass.alu.annot.filt.tmp.txt

mv {OUTPUT_PREFIX}.juncmut.sjclass.alu.annot.filt.tmp.txt {OUTPUT_PREFIX}.juncmut.sjclass.alu.annot.filt.txt

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
        (SECTION_NAME, "gnomad"),
    ]
    optional_dbs = [
        (SECTION_NAME, "cgc_file"),
        (SECTION_NAME, "clinvar_file"),
        (SECTION_NAME, "acmg_file"),
        (SECTION_NAME, "clinvar_star234_file"),
        (SECTION_NAME, "pancan_file"),
        (SECTION_NAME, "dosage_sensitivity_file"),
        (SECTION_NAME, "cgd_file"),
    ]
    local_dbs = []
    for (section, db_name) in dbs:
        parsed = urllib.parse.urlparse(gcat_conf.get(section, db_name))
        if parsed.scheme == "":
            local_dbs.append(gcat_conf.path_get(section, db_name))

    juncmut_annot_param = ["--gnomad %s" % (gcat_conf.path_get(SECTION_NAME, "gnomad"))]
    for (section, db_name) in optional_dbs:
        value = gcat_conf.safe_get(section, db_name, "")
        if value != "":
            local_dbs.append(gcat_conf.path_get(section, db_name))
            juncmut_annot_param += ["--%s %s" % (db_name, gcat_conf.path_get(section, db_name))]

    for sample in sample_conf.juncmut:
        output_dir = "%s/juncmut/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)
        output_files[sample] = [
            "%s/%s.juncmut.sjclass.alu.annot.filt.txt" % (output_dir, sample),
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
            "SJ_OUTTAB": sjtab,
            "BAM": input_bams[sample],
            "OUTPUT_DIR": output_dir,
            "OUTPUT_PREFIX": "%s/%s" % (output_dir, sample),
            "JUNCMUT_DETECT_PARAM": gcat_conf.get(SECTION_NAME, "juncmut_detect_param"),
            "JUNCMUT_SJCLASS_PARAM": gcat_conf.get(SECTION_NAME, "juncmut_sjclass_param"),
            "JUNCMUT_ALU_PARAM": gcat_conf.get(SECTION_NAME, "juncmut_alu_param"),
            "JUNCMUT_ANNOT_PARAM": " ".join(juncmut_annot_param),
            "REFERENCE": gcat_conf.path_get(SECTION_NAME, "reference"),
            "CONTROL_FILE1": gcat_conf.path_get(SECTION_NAME, "control_file1"),
            "CONTROL_FILE2": gcat_conf.path_get(SECTION_NAME, "control_file2"),
            "GENCODE": gcat_conf.path_get(SECTION_NAME, "genecode_gene_file"),
            "RMSK_BED": gcat_conf.path_get(SECTION_NAME, "rmsk_bed"),
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
