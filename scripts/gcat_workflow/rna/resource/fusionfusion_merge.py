#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

OUTPUT_FORMAT = "fusionfusion/control_panel/{sample}.merged.Chimeric.count"

class Fusionfusion_merge(stage_task.Stage_task):
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

echo -n > {MERGED_LIST}
for TMP_SAMPLE in {PANEL_SAMPLES}; do
    echo {OUTPUT_DIR}/${{TMP_SAMPLE}}/${{TMP_SAMPLE}}.Chimeric.count >> {MERGED_LIST}
done

chimera_utils merge_control {OPTION} {MERGED_LIST} {OUTPUT}
"""

def configure(gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "fusionfusion_merge"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": gcat_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Fusionfusion_merge(params)
    
    output_files = {}
    for (sample_dummy, panel) in sample_conf.fusionfusion:
        if panel == None:
            continue 
        if panel in output_files:
            continue
        
        output_file = OUTPUT_FORMAT.format(sample = panel)
        output_files[panel] = run_conf.project_root + "/" + output_file
        
        output_dir = "%s/fusionfusion/control_panel" % (run_conf.project_root)
        os.makedirs(output_dir, exist_ok=True)

        arguments = {
            "MERGED_LIST": "%s/%s.Chimeric.count.list" % (output_dir, panel),
            "OUTPUT": "%s/%s.merged.Chimeric.count" % (output_dir, panel),
            "PANEL_SAMPLES": " ".join(sample_conf.control_panel[panel]),
            "OUTPUT_DIR": "%s/fusionfusion" % (run_conf.project_root),
            "OPTION": gcat_conf.get(SECTION_NAME, "chimera_utils_merge_option")
        }
        
        singularity_bind = [run_conf.project_root]
        for sample in sample_conf.control_panel[panel]:
            if sample in sample_conf.bam_import_src:
                singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = panel)

    return output_files
