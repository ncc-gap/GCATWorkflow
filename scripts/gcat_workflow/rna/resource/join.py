#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

class Join(stage_task.Stage_task):
    def __init__(self, params):
        super().__init__(params)
        self.shell_script_template = """
#!/bin/bash
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

{RM_FASTQS}
touch {JOIN_FILE}
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(remove_files, gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "join"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": "",
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": ""
    }
    stage_class = Join(params)
    
    output_dir = "%s/join" % (run_conf.project_root)
    os.makedirs(output_dir, exist_ok=True)
    fastqs = []
    for sample in remove_files:
        for li in remove_files[sample]:
            fastqs.extend(li)
    
    rm_fastqs = ""
    if bool(gcat_conf.get(SECTION_NAME, "qsub_option")) and len(fastqs) > 0:
        rm_fastqs = "rm -f " +  " ".join(fastqs)

    arguments = {
        "RM_FASTQS": rm_fastqs,
        "JOIN_FILE": "%s/all.txt" % (output_dir)
    }
    
    singularity_bind = []
    
    stage_class.write_script(arguments, singularity_bind, run_conf)

    return []
