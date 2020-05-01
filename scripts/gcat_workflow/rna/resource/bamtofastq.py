#! /usr/bin/env python

import gcat_workflow.core.stage_task_abc as stage_task

class Bam_tofastq(stage_task.Stage_task):
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

export LD_LIBRARY_PATH=/usr/local/lib
bams=( `echo "{input_bam}" | tr -s ';' ' '`)
if [ ${{#bams[@]}} -eq 1 ]; then
    bam=${{bams[0]}}
    echo $bam
    /usr/local/bin/bamtofastq {param} filename=${{bam}} F={f1_name}.tmp F2={f2_name}.tmp T={t} S={s}.tmp O={o1_name} O2={o2_name}
    if [ -s {f1_name}.tmp ]; then
        mv {f1_name}.tmp {f1_name}
        mv {f2_name}.tmp {f2_name}
        mv {s}.tmp {s}
    elif [ -s {s}.tmp ]; then
        mv {s}.tmp {f1_name}
        rm {f1_name}.tmp {f2_name}.tmp
    fi
else
    touch {f1_name} {f2_name} {t} {s} {o1_name} {o2_name}
    
    for bam in ${{bams[@]}}; do
        echo $bam
        /usr/local/bin/bamtofastq {param} filename=${{bam}} F={f1_name}.tmp F2={f2_name}.tmp T={t}.tmp S={s}.tmp O={o1_name}.tmp O2={o2_name}.tmp
        if [ -s {f1_name}.tmp ]; then
            cat {f1_name}.tmp >> {f1_name}
            cat {f2_name}.tmp >> {f2_name}
            cat {s}.tmp >> {s}
        elif [ -s {s}.tmp ]; then
            cat {s}.tmp >> {f1_name}
        fi
            
        cat {t}.tmp >> {t}
        cat {o1_name}.tmp >> {o1_name}
        cat {o2_name}.tmp >> {o2_name}
        
        rm {f1_name}.tmp {f2_name}.tmp {t}.tmp {s}.tmp {o1_name}.tmp {o2_name}.tmp
    done
fi
touch {pass_file}
"""

# merge sorted bams into one and mark duplicate reads with biobambam
def configure(gcat_conf, run_conf, sample_conf):
    import os
    
    STAGE_NAME = "bam_tofastq"
    CONF_SECTION = "bam_tofastq"
    params = {
        "work_dir": run_conf.project_root,
        "stage_name": STAGE_NAME,
        "image": gcat_conf.path_get(CONF_SECTION, "image"),
        "qsub_option": gcat_conf.get(CONF_SECTION, "qsub_option"),
        "singularity_option": gcat_conf.get(CONF_SECTION, "singularity_option")
    }
    stage_class = Bam_tofastq(params)
    output_fastqs = {}
    for sample in sample_conf.bam_tofastq:
        output_dir = "%s/fastq/%s" % (run_conf.project_root, sample)
        os.makedirs(output_dir, exist_ok=True)    
        f1_name = output_dir + "/1_1.fastq"
        f2_name = output_dir + "/1_2.fastq"
        output_fastqs[sample] = [[f1_name], [f2_name]]

        arguments = {
            "param": gcat_conf.get(CONF_SECTION, "params"),
            "input_bam": sample_conf.bam_tofastq[sample],
            "f1_name": f1_name,
            "f2_name": f2_name,
            "o1_name": output_dir + '/unmatched_first_output.txt',
            "o2_name": output_dir + '/unmatched_second_output.txt',
            "t": output_dir + '/temp.txt',
            "s": output_dir + '/single_end_output.txt',
            "pass_file": output_dir + '/pass.txt'
        }
        
        singularity_bind = [
            run_conf.project_root,
        ] + sample_conf.bam_tofastq_src[sample]
        
        stage_class.write_script(arguments, singularity_bind, run_conf, sample = sample)
    return output_fastqs
