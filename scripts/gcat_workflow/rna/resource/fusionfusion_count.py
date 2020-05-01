#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

OUTPUT_FORMAT = "fusionfusion/{sample}/{sample}.Chimeric.count"

class Fusionfusion_count(stage_task.Stage_task):
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

chimera_utils count {OPTION} {INPUT} {OUTPUT}
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(input_files, gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "fusionfusion_count"
    SECTION_NAME = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(SECTION_NAME, "image"),
        "qsub_option": gcat_conf.get(SECTION_NAME, "qsub_option"),
        "singularity_option": gcat_conf.get(SECTION_NAME, "singularity_option")
    }
    stage_class = Fusionfusion_count(params)
    
    samples = []
    for (sample, panel) in sample_conf.fusionfusion:
        samples.append(sample)
        if panel == None:
            continue
        for i in sample_conf.control_panel[panel]:
            if not i in samples:
                samples.append(i)
                
    output_files = {}
    for sample in samples:
        output_dir = "%s/fusionfusion/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)
        output_files[sample] = OUTPUT_FORMAT.format(sample = sample)
        
        arguments = {
            "INPUT": input_files[sample],
            "OUTPUT": "%s/%s.Chimeric.count" % (output_dir, sample),
            "OPTION": gcat_conf.get(SECTION_NAME, "chimera_utils_count_option"),
        }
       
        singularity_bind = [run_conf.project_root]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)

    return output_files
