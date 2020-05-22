#! /usr/bin/env python

import os

# link files attached bam
def link_import_attached_files(run_conf, bam_import_stage, bam_postfix, junction_postfix, subdir = "star"):
    for sample in bam_import_stage:
        source_file = bam_import_stage[sample].replace(bam_postfix, junction_postfix)
        if not os.path.exists(source_file):
            raise ValueError("Not exist junction file: %s" % (source_file))
        
        link_dir = "%s/%s/%s" % (run_conf.project_root, subdir, sample)
        os.makedirs(link_dir, exist_ok=True)
        
        link_file = link_dir +'/'+ sample + junction_postfix
        if not os.path.exists(link_file):
            os.symlink(source_file, link_file)

def main(gcat_conf, run_conf, sample_conf):
    
    # preparation
    import gcat_workflow.core.setup_common as setup
    input_stages = (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq)
    setup.create_directories(gcat_conf, run_conf, input_stages, 'rna/data/snakefile.txt')
    bam_tofastq_stages = (sample_conf.bam_tofastq, )
    setup.touch_bam_tofastq(run_conf, bam_tofastq_stages)
    
    # dump conf.yaml
    import gcat_workflow.rna.resource.star_align as rs_align
    y = setup.dump_yaml_input_section(
        run_conf,
        (sample_conf.bam_tofastq, ),
        sample_conf.fastq,
        sample_conf.bam_import, 
        rs_align.OUTPUT_BAM_FORMAT
    )

    # link fastq
    linked_fastqs = setup.link_input_fastq(run_conf, sample_conf.fastq, sample_conf.fastq_src)
    
    # link import bam
    output_bams = setup.link_import_bam(
        run_conf, sample_conf.bam_import, 
        rs_align.BAM_POSTFIX, 
        rs_align.BAI_POSTFIX,
        rs_align.OUTPUT_BAM_FORMAT.split("/")[0]
    )
    
    # link files attached bam
    link_import_attached_files(
        run_conf, sample_conf.bam_import,
        rs_align.BAM_POSTFIX, 
        rs_align.CHIMERIC_JUNCTION_POSTFIX, 
        rs_align.OUTPUT_BAM_FORMAT.split("/")[0]
    )
    link_import_attached_files(
        run_conf, sample_conf.bam_import,
        rs_align.BAM_POSTFIX, 
        rs_align.CHIMERIC_SAM_POSTFIX, 
        rs_align.OUTPUT_BAM_FORMAT.split("/")[0]
    )
    
    # ######################
    # create scripts
    # ######################
    # bam to fastq
    import gcat_workflow.rna.resource.bamtofastq as rs_bamtofastq
    output_fastqs = rs_bamtofastq.configure(gcat_conf, run_conf, sample_conf)
    
    # star
    for sample in output_fastqs:
        sample_conf.fastq[sample] = output_fastqs[sample]
        sample_conf.fastq_src[sample] = [[], []]

    for sample in linked_fastqs:
        sample_conf.fastq[sample] = linked_fastqs[sample]["fastq"]
        sample_conf.fastq_src[sample] = linked_fastqs[sample]["src"]

    align_bams = rs_align.configure(gcat_conf, run_conf, sample_conf)
    output_bams.update(align_bams)

    # fusion-fusion
    output_bam_sams = {}
    for sample in output_bams:
        output_bam_sams[sample] = output_bams[sample].replace(rs_align.BAM_POSTFIX, rs_align.CHIMERIC_SAM_POSTFIX)

    import gcat_workflow.rna.resource.fusionfusion_count as rs_fusionfusion_count
    output_fusionfusion_counts = rs_fusionfusion_count.configure(output_bam_sams, gcat_conf, run_conf, sample_conf)
    
    import gcat_workflow.rna.resource.fusionfusion_merge as rs_fusionfusion_merge
    output_fusionfusion_merges = rs_fusionfusion_merge.configure(output_fusionfusion_counts, gcat_conf, run_conf, sample_conf)
    
    import gcat_workflow.rna.resource.fusionfusion as rs_fusionfusion
    output_fusionfusions = rs_fusionfusion.configure(output_bam_sams, output_fusionfusion_merges, gcat_conf, run_conf, sample_conf)
    
    # STAR-fusion
    output_bam_junctions = {}
    for sample in output_bams:
        output_bam_junctions[sample] = output_bams[sample].replace(rs_align.BAM_POSTFIX, rs_align.CHIMERIC_JUNCTION_POSTFIX)

    import gcat_workflow.rna.resource.star_fusion as rs_star_fusion
    output_star_fusions = rs_star_fusion.configure(output_bam_junctions, gcat_conf, run_conf, sample_conf)
    
    # ir_count
    import gcat_workflow.rna.resource.ir_count as rs_ir_count
    output_ir_counts = rs_ir_count.configure(output_bams, gcat_conf, run_conf, sample_conf)
    
    # iravnet
    import gcat_workflow.rna.resource.iravnet as rs_iravnet
    output_iravnets = rs_iravnet.configure(output_bams, gcat_conf, run_conf, sample_conf)
    
    # expression
    import gcat_workflow.rna.resource.expression as rs_expression
    output_expressions = rs_expression.configure(output_bams, gcat_conf, run_conf, sample_conf)
    
    # kallisto
    import gcat_workflow.rna.resource.kallisto as rs_kallisto
    output_kallistos = rs_kallisto.configure(gcat_conf, run_conf, sample_conf)
    
    # join
    import gcat_workflow.rna.resource.join as rs_join
    rs_join.configure(output_fastqs, gcat_conf, run_conf, sample_conf)

    # ######################
    # dump conf.yaml
    # ######################
    def __dic_values(dic):
        values = []
        for key in dic:
            if type(dic[key]) == list:
                values.extend(dic[key])
            else:
                values.append(dic[key])
        return values

    #y["output_files"].extend(__dic_values(output_fusionfusion_counts))
    #y["output_files"].extend(__dic_values(output_fusionfusion_merges))
    y["output_files"].extend(output_fusionfusions)
    y["output_files"].extend(output_star_fusions)
    y["output_files"].extend(output_ir_counts)
    y["output_files"].extend(output_iravnets)
    y["output_files"].extend(output_expressions)
    y["output_files"].extend(output_kallistos)
    
    y["fusionfusion_count_samples"] = {}
    for [sample, panel] in sample_conf.fusionfusion:
        y["fusionfusion_count_samples"][sample] = rs_align.OUTPUT_CHIMERIC_SAM_FORMAT.format(sample=sample)
        if panel == None:
            continue
        for i in sample_conf.control_panel[panel]:
            y["fusionfusion_count_samples"][i] = rs_align.OUTPUT_CHIMERIC_SAM_FORMAT.format(sample=i)

    y["fusionfusion_merge_samples"] = {}
    for [sample, panel] in sample_conf.fusionfusion:
        if panel == None:
            continue
        y["fusionfusion_merge_samples"][panel] = []
        for i in sample_conf.control_panel[panel]:
            y["fusionfusion_merge_samples"][panel].append(rs_fusionfusion_count.OUTPUT_FORMAT.format(sample=i))
    
    y["fusionfusion_samples"] = {}
    for [sample, panel] in sample_conf.fusionfusion:
        y["fusionfusion_samples"][sample] = [
            rs_align.OUTPUT_CHIMERIC_SAM_FORMAT.format(sample=sample)
        ]
        if panel != None:
            y["fusionfusion_samples"][sample].append(rs_fusionfusion_merge.OUTPUT_FORMAT.format(sample=panel))
        
    
    y["star_fusion_samples"] = {}
    for sample in sample_conf.star_fusion:
        y["star_fusion_samples"][sample] = rs_align.OUTPUT_CHIMERIC_JUNCTION_FORMAT.format(sample=sample)
    
    y["expression_samples"] = {}
    for sample in sample_conf.expression:
        y["expression_samples"][sample] = rs_align.OUTPUT_BAM_FORMAT.format(sample=sample)
    
    y["ir_count_samples"] = {}
    for sample in sample_conf.ir_count:
        y["ir_count_samples"][sample] = rs_align.OUTPUT_BAM_FORMAT.format(sample=sample)

    y["iravnet_samples"] = {}
    for sample in sample_conf.iravnet:
        y["iravnet_samples"][sample] = rs_align.OUTPUT_BAM_FORMAT.format(sample=sample)

    y["kallisto_samples"] = {}
    for sample in sample_conf.kallisto:
        if sample in y["aln_samples"]:
            y["kallisto_samples"][sample] = y["aln_samples"][sample]
        else:
            y["kallisto_samples"][sample] = y["bam_import"][sample]
            
    import yaml
    open(run_conf.project_root + "/config.yml", "w").write(yaml.dump(y))
    
