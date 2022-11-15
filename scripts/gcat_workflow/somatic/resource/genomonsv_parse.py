#! /usr/bin/env python

import os
import gcat_workflow.core.stage_task_abc as stage_task

class genomonsv_parse(stage_task.Stage_task):
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

INPUT_BAM={input_cram}

if [ "{use_bam}" = "T" ]
then
  samtools view \\
    {samtools_view_option} \\
    -b \\
    -T {reference} \\
    -o {temp_bam} \\
    {input_cram}

  samtools index \\
    {samtools_index_option} \\
    {temp_bam}

  INPUT_BAM={temp_bam}
fi

GenomonSV parse ${{INPUT_BAM}} {output_prefix} --reference {reference} {param}

if [ "{use_bam}" = "T" ]
then
  rm {temp_bam}.bai
  rm {temp_bam}
fi
"""

def configure(input_bams, gcat_conf, run_conf, sample_conf):
    
    STAGE_NAME = "genomonsv_parse"
    CONF_SECTION = STAGE_NAME
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = genomonsv_parse(params)

    samples = []
    for (tumor, normal, panel) in sample_conf.genomon_sv:
        samples.append(tumor)
        if normal != None:
            samples.append(normal)
        if panel != None:
            samples.extend(sample_conf.control_panel[panel])
    samples = list(set(samples))

    output_files = {}
    for sample in samples:
        output_dir = '%s/genomonsv/%s' % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)
        output_files[sample] = "%s/%s.junction.clustered.bedpe.gz.tbi"% (output_dir, sample)

        use_bam = "F"
        temp_bam = input_bams[sample]
        if gcat_conf.getboolean(CONF_SECTION, "use_bam"):
            use_bam = "T"
            temp_bam = "%s/%s.bam" % (output_dir, sample)

        arguments = {
            "input_cram": input_bams[sample],
            "output_prefix": "%s/%s"% (output_dir, sample),
            "param": gcat_conf.get(CONF_SECTION, "params"),
            "reference": gcat_conf.path_get(CONF_SECTION, "reference"),
            "use_bam": use_bam,
            "temp_bam": temp_bam,
            "samtools_view_option": gcat_conf.get(CONF_SECTION, "samtools_view_option") + " " + gcat_conf.get(CONF_SECTION, "samtools_view_threads_option"),
            "samtools_index_option": gcat_conf.get(CONF_SECTION, "samtools_index_option") + " " + gcat_conf.get(CONF_SECTION, "samtools_index_threads_option"),
        }
        
        singularity_bind = [run_conf.project_root, os.path.dirname(gcat_conf.path_get(CONF_SECTION, "reference"))]
        if sample in sample_conf.bam_import_src:
            singularity_bind += sample_conf.bam_import_src[sample]
            
        stage_class.write_script(arguments, singularity_bind, run_conf, gcat_conf, sample = sample)
    
    return output_files
