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

def touch_sra_fastq_dump(run_conf, sra_fastq_dump_stage):
    for sample in sra_fastq_dump_stage:
        wdir = run_conf.project_root + '/sra_fastq_dump/' + sample
        os.makedirs(wdir, exist_ok=True)
        open(wdir + '/' + sample + ".txt", "w").close()

def main(gcat_conf, run_conf, sample_conf):
    
    # preparation
    import gcat_workflow.core.setup_common as setup
    input_stages = (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq, sample_conf.sra_fastq_dump)
    setup.create_directories(gcat_conf, run_conf, input_stages, 'rna/data/snakefile.txt')
    setup.touch_bam_tofastq(run_conf, (sample_conf.bam_tofastq,))
    touch_sra_fastq_dump(run_conf,  sample_conf.sra_fastq_dump)

    # dump conf.yaml
    import gcat_workflow.rna.resource.star_align as rs_align
    y = setup.dump_yaml_input_section(
        run_conf,
        (sample_conf.bam_tofastq, ),
        sample_conf.fastq,
        sample_conf.bam_import, 
        rs_align.OUTPUT_BAM_FORMAT,
        gcat_conf.get("join", "remove_bam").lower() == "true"
    )

    #import json
    #print(json.dumps(y, indent=4, sort_keys=True, separators=(',', ': ')))

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
    #link_import_attached_files(
    #    run_conf, 
    #    sample_conf.bam_import, 
    #    rs_align.BAM_POSTFIX, 
    #    rs_align.CHIMERIC_JUNCTION_POSTFIX, 
    #    rs_align.OUTPUT_BAM_FORMAT.split("/")[0]
    #)
    link_import_attached_files(
        run_conf, sample_conf.bam_import,
        rs_align.BAM_POSTFIX, 
        rs_align.CHIMERIC_SAM_POSTFIX, 
        rs_align.OUTPUT_BAM_FORMAT.split("/")[0]
    )
    link_import_attached_files(
        run_conf, sample_conf.bam_import,
        rs_align.BAM_POSTFIX, 
        rs_align.SJ_TAB_POSTFIX, 
        rs_align.OUTPUT_BAM_FORMAT.split("/")[0]
    )
    
    # ######################
    # create scripts
    # ######################
    # bam to fastq
    import gcat_workflow.rna.resource.bamtofastq as rs_bamtofastq
    output_fastqs = rs_bamtofastq.configure(gcat_conf, run_conf, sample_conf)
    
    # SRA fastq-dump
    import gcat_workflow.rna.resource.sra_fastq_dump as rs_sra_fastq_dump
    output_sra_fastq_dump = rs_sra_fastq_dump.configure(gcat_conf, run_conf, sample_conf)
    output_fastqs.update(output_sra_fastq_dump)
    
    # cram to bam
    import gcat_workflow.rna.resource.cram_tobam as rs_cramtobam
    cramto_bams = rs_cramtobam.configure(gcat_conf, run_conf, sample_conf)
    output_bams.update(cramto_bams)

    # star
    org_fastq_samples = list(sample_conf.fastq.keys())
    for sample in output_fastqs:
        sample_conf.fastq[sample] = output_fastqs[sample]
        sample_conf.fastq_src[sample] = [[], []]

    for sample in linked_fastqs:
        sample_conf.fastq[sample] = linked_fastqs[sample]["fastq"]
        sample_conf.fastq_src[sample] = linked_fastqs[sample]["src"]

    align_bams = rs_align.configure(gcat_conf, run_conf, sample_conf, org_fastq_samples)
    output_bams.update(align_bams)

    # fusion-fusion
    output_bam_sams = {}
    for sample in output_bams:
        output_bam_sams[sample] = output_bams[sample].replace(rs_align.BAM_POSTFIX, rs_align.CHIMERIC_SAM_POSTFIX)

    import gcat_workflow.rna.resource.fusionfusion_count as rs_fusionfusion_count
    rs_fusionfusion_count.configure(output_bam_sams, gcat_conf, run_conf, sample_conf)
    
    import gcat_workflow.rna.resource.fusionfusion_merge as rs_fusionfusion_merge
    output_fusionfusion_merges = rs_fusionfusion_merge.configure(gcat_conf, run_conf, sample_conf)

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
    
    # juncmut
    output_sj_tabs = {}
    for sample in output_bams:
        output_sj_tabs[sample] = output_bams[sample].replace(rs_align.BAM_POSTFIX, rs_align.SJ_TAB_POSTFIX)

    import gcat_workflow.rna.resource.juncmut as rs_juncmut
    output_juncmuts = rs_juncmut.configure(output_bams, output_sj_tabs, gcat_conf, run_conf, sample_conf)
    
    # expression
    import gcat_workflow.rna.resource.expression as rs_expression
    output_expressions = rs_expression.configure(output_bams, gcat_conf, run_conf, sample_conf)
    
    # kallisto
    import gcat_workflow.rna.resource.kallisto as rs_kallisto
    output_kallistos = rs_kallisto.configure(gcat_conf, run_conf, sample_conf)
    
    # join
    import gcat_workflow.rna.resource.join as rs_join
    bams = {}
    bams.update(cramto_bams)
    bams.update(align_bams)
    output_joins = rs_join.configure(
        output_bams.keys(), 
        list(sample_conf.bam_import.keys()) + list(sample_conf.cram_import.keys()), 
        output_fastqs, 
        bams, 
        gcat_conf, run_conf, sample_conf
    )

    # ######################
    # dump conf.yaml
    # ######################
    def __to_relpath(fullpath):
        return fullpath.replace(run_conf.project_root + "/", "", 1)
  
    def __dic_values(dic):
        values = []
        for key in dic:
            if type(dic[key]) == list:
                for path in dic[key]:
                    values.append(__to_relpath(path))
            else:
                values.append(__to_relpath(dic[key]))
        return values

    def __update_dic(target, src):
        for key in src:
            if not key in target:
                target[key] = []
            if type(src[key]) == list:
                for path in src[key]:
                    target[key].append(__to_relpath(path))
            else:
                target[key].append(__to_relpath(src[key]))

    y["sra_fastq_dump"] = {}
    output_dumps = {}
    for sample in sample_conf.sra_fastq_dump:
        y["sra_fastq_dump"][sample] = "sra_fastq_dump/%s/%s.txt" % (sample, sample)
        y["aln_samples"][sample] = "fastq/%s/pass.txt" % (sample)
        output_dumps[sample] = rs_align.OUTPUT_BAM_FORMAT.format(sample=sample)

    y["cram_tobam"] = {}
    for sample in sample_conf.cram_import:
        y["cram_tobam"][sample] = sample_conf.cram_import[sample]

    output_files = y.pop("output_files")
    dic_output_files = {}
    for path in output_files:
        sample = path.split("/")[-2]
        if not sample in dic_output_files:
            dic_output_files[sample] = []
        dic_output_files[sample].append(path)

    output_stars = {}
    output_cramto_bams = {}
    for sample in output_bams:
        outputs = [
            output_bams[sample],
            output_bams[sample].replace(rs_align.BAM_POSTFIX, rs_align.BAI_POSTFIX),
            #output_bams[sample].replace(rs_align.BAM_POSTFIX, rs_align.CHIMERIC_JUNCTION_POSTFIX),
            output_bams[sample].replace(rs_align.BAM_POSTFIX, rs_align.CHIMERIC_SAM_POSTFIX),
            output_bams[sample].replace(rs_align.BAM_POSTFIX, rs_align.SJ_TAB_POSTFIX),
        ]
        if sample in sample_conf.cram_import:
            output_cramto_bams[sample] = outputs
        else:
            output_stars[sample] = outputs

    output_fusionfusion_merges_sample = {}
    for panel in output_fusionfusion_merges:
        for sample in sample_conf.control_panel[panel]:
            output_fusionfusion_merges_sample[sample] = output_fusionfusion_merges[panel]

    __update_dic(dic_output_files, output_stars)
    __update_dic(dic_output_files, output_cramto_bams)
    __update_dic(dic_output_files, output_dumps)
    __update_dic(dic_output_files, output_fusionfusion_merges_sample)
    __update_dic(dic_output_files, output_fusionfusions)
    __update_dic(dic_output_files, output_star_fusions)
    __update_dic(dic_output_files, output_ir_counts)
    __update_dic(dic_output_files, output_iravnets)
    __update_dic(dic_output_files, output_juncmuts)
    __update_dic(dic_output_files, output_expressions)
    __update_dic(dic_output_files, output_kallistos)
    y["output_files"] = dic_output_files
    y["join"] = __dic_values(output_joins)

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

    y["juncmut_samples"] = {}
    for sample in sample_conf.juncmut:
        y["juncmut_samples"][sample] = rs_align.OUTPUT_SJ_TAB_FORMAT.format(sample=sample)

    y["kallisto_samples"] = {}
    for sample in sample_conf.kallisto:
        if sample in y["aln_samples"]:
            y["kallisto_samples"][sample] = y["aln_samples"][sample]
        else:
            y["kallisto_samples"][sample] = y["bam_import"][sample]

    import yaml
    open(run_conf.project_root + "/config.yml", "w").write(yaml.dump(y))
    
