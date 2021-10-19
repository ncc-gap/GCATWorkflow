#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

BAM_POSTFIX = ".Aligned.sortedByCoord.out.bam"
BAI_POSTFIX = ".Aligned.sortedByCoord.out.bam.bai"
CRAM_POSTFIX = ".Aligned.sortedByCoord.out.cram"
CHIMERIC_JUNCTION_POSTFIX = ".Chimeric.out.junction"
CHIMERIC_SAM_POSTFIX = ".Chimeric.out.sam"
SJ_TAB_POSTFIX = ".SJ.out.tab.gz"

OUTPUT_BAM_FORMAT = "star/{sample}/{sample}" + BAM_POSTFIX
OUTPUT_CHIMERIC_JUNCTION_FORMAT = "star/{sample}/{sample}" + CHIMERIC_JUNCTION_POSTFIX
OUTPUT_CHIMERIC_SAM_FORMAT = "star/{sample}/{sample}" + CHIMERIC_SAM_POSTFIX
OUTPUT_SJ_TAB_FORMAT = "star/{sample}/{sample}" + SJ_TAB_POSTFIX

class Cram_tobam(stage_task.Stage_task):
    def __init__(self, params):
        super().__init__(params)
        self.shell_script_template = """#! /bin/bash
set -eux

rm -rf {OUTPUT_DIR}/*

samtools view -b {OPTION} -T {REFERENCE} {INPUT_CRAM} -o {OUTPUT_BAM}
samtools index {OUTPUT_BAM}

ln -s {INPUT_CHIMERIC_SAM} {OUTPUT_CHIMERIC_SAM}
ln -s {INPUT_SJ_TAB} {OUTPUT_SJ_TAB}
"""

def configure(gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "cram_tobam"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": gcat_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Cram_tobam(params)
    
    output_bams = {}
    for sample in sample_conf.cram_import:
        input_dir = os.path.dirname(sample_conf.cram_import[sample])
        input_chimeric_sam = sample_conf.cram_import[sample].replace(CRAM_POSTFIX, CHIMERIC_SAM_POSTFIX)
        if not os.path.exists(input_chimeric_sam):
            raise ValueError("Not exist junction file: %s" % (input_chimeric_sam))
        input_sj_tab =  sample_conf.cram_import[sample].replace(CRAM_POSTFIX, SJ_TAB_POSTFIX)
        if not os.path.exists(input_sj_tab):
            raise ValueError("Not exist junction file: %s" % (input_sj_tab))
        
        output_dir = "%s/star/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)

        output_bam = "%s/%s%s" % (output_dir, sample, BAM_POSTFIX)
        output_sam = "%s/%s%s" % (output_dir, sample, CHIMERIC_SAM_POSTFIX)
        output_sjtab = "%s/%s%s" % (output_dir, sample, SJ_TAB_POSTFIX)

        output_bams[sample] = output_bam
 
        arguments = {
            "INPUT_CRAM": sample_conf.cram_import[sample],
            "INPUT_CHIMERIC_SAM": input_chimeric_sam,
            "INPUT_SJ_TAB": input_sj_tab,
            "OUTPUT_DIR": output_dir,
            "OUTPUT_BAM": output_bam,
            "OUTPUT_CHIMERIC_SAM": output_sam,
            "OUTPUT_SJ_TAB": output_sjtab,
            "REFERENCE": gcat_conf.path_get(SECTION_NAME, "reference"),
            "OPTION": " ".join([
                gcat_conf.get(SECTION_NAME, "samtools_option"),
                gcat_conf.get(SECTION_NAME, "samtools_threads_option"),
            ]),
        }
        
        singularity_bind = [
            run_conf.project_root,
            gcat_conf.get(SECTION_NAME, "reference"),
        ] + sample_conf.cram_import_src[sample]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)
    return output_bams
