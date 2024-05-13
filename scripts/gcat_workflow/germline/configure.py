#! /usr/bin/env python

def main(gcat_conf, run_conf, sample_conf):
    
    # preparation
    import gcat_workflow.core.setup_common as setup
    input_stages = (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq)
    setup.create_directories(gcat_conf, run_conf, input_stages, 'germline/data/snakefile.txt')
    setup.touch_bam_tofastq(run_conf, (sample_conf.bam_tofastq,))
    
    # dump conf.yaml
    import gcat_workflow.germline.resource.post_align as rs_post_align
    y = setup.dump_yaml_input_section(
        run_conf, 
        (sample_conf.bam_tofastq, ),
        sample_conf.fastq,
        sample_conf.bam_import, 
        rs_post_align.OUTPUT_FORMAT
    )

    # link fastq
    linked_fastqs = setup.link_input_fastq(run_conf, sample_conf.fastq, sample_conf.fastq_src)
    
    # link import bam
    output_bams = setup.link_import_bam(
        run_conf,
        sample_conf.bam_import, 
        rs_post_align.BAM_POSTFIX, 
        rs_post_align.BAI_POSTFIX,
        rs_post_align.OUTPUT_FORMAT.split("/")[0]
    )
    
    # ######################
    # create scripts
    # ######################
    # bam to fastq
    import gcat_workflow.germline.resource.bamtofastq as rs_bamtofastq
    output_fastqs = rs_bamtofastq.configure(gcat_conf, run_conf, sample_conf)
    
    # bwa
    for sample in output_fastqs:
        sample_conf.fastq[sample] = output_fastqs[sample]
        sample_conf.fastq_src[sample] = [[], []]

    for sample in linked_fastqs:
        sample_conf.fastq[sample] = linked_fastqs[sample]["fastq"]
        sample_conf.fastq_src[sample] = linked_fastqs[sample]["src"]
    
    import gcat_workflow.germline.resource.parabricks_align as rs_align
    align_bams = rs_align.configure(gcat_conf, run_conf, sample_conf)
    post_align_bams = rs_post_align.configure(align_bams, gcat_conf, run_conf, sample_conf)
    output_bams.update(post_align_bams)

    # mutation
    import gcat_workflow.germline.resource.haplotypecaller as rs_mutation
    output_mutations = rs_mutation.configure(output_bams, gcat_conf, run_conf, sample_conf)
    
    # wgs summary
    import gcat_workflow.germline.resource.collectwgsmetrics as rs_wgs_summary
    output_wgs_metrics = rs_wgs_summary.configure(output_bams, gcat_conf, run_conf, sample_conf)

    # hs summary
    import gcat_workflow.somatic.resource.collecthsmetrics as rs_hs_summary
    output_hs_metrics = rs_hs_summary.configure(output_bams, gcat_conf, run_conf, sample_conf)

    # multiple summary
    import gcat_workflow.germline.resource.collectmultiplemetrics as rs_multiple_summary
    output_multiple_metrics = rs_multiple_summary.configure(output_bams, gcat_conf, run_conf, sample_conf)
    
    # gridss
    import gcat_workflow.germline.resource.gridss as rs_gridss
    output_gridss = rs_gridss.configure(output_bams, gcat_conf, run_conf, sample_conf)

    # manta 
    import gcat_workflow.germline.resource.manta as rs_manta
    output_manta = rs_manta.configure(output_bams, gcat_conf, run_conf, sample_conf)

    # melt
    import gcat_workflow.germline.resource.melt as rs_melt
    output_melt = rs_melt.configure(output_bams, gcat_conf, run_conf, sample_conf)

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

    y["output_files"].extend(__dic_values(output_mutations))
    y["output_files"].extend(__dic_values(output_wgs_metrics))
    y["output_files"].extend(__dic_values(output_hs_metrics))
    y["output_files"].extend(__dic_values(output_multiple_metrics))
    y["output_files"].extend(__dic_values(output_gridss))
    y["output_files"].extend(__dic_values(output_manta))
    y["output_files"].extend(__dic_values(output_melt))

    y["post_aln_samples"] = {}
    for sample in y["aln_samples"]:
        y["post_aln_samples"][sample] = rs_align.OUTPUT_FORMAT.format(sample=sample)
    
    y["haplotypecaller_samples"] = {}
    for sample in sample_conf.haplotype_call:
        y["haplotypecaller_samples"][sample] = rs_post_align.OUTPUT_FORMAT.format(sample=sample)
        
    y["collect_wgs_metrics_samples"] = {}
    for sample in sample_conf.wgs_metrics:
        y["collect_wgs_metrics_samples"][sample] = rs_post_align.OUTPUT_FORMAT.format(sample=sample)

    y["collect_multiple_metrics_samples"] = {}
    for sample in sample_conf.multiple_metrics:
        y["collect_multiple_metrics_samples"][sample] = rs_post_align.OUTPUT_FORMAT.format(sample=sample)
        
    y["collect_hs_metrics_samples"] = {}
    for sample in sample_conf.hs_metrics:
        y["collect_hs_metrics_samples"][sample] = rs_post_align.OUTPUT_FORMAT.format(sample=sample)
    
    y["gridss_samples"] = {}
    for sample in sample_conf.gridss:
        y["gridss_samples"][sample] = rs_post_align.OUTPUT_FORMAT.format(sample=sample)
        
    y["manta_samples"] = {}
    for sample in sample_conf.manta:
        y["manta_samples"][sample] = rs_post_align.OUTPUT_FORMAT.format(sample=sample)
        
    y["melt_samples"] = {}
    for sample in sample_conf.melt:
        y["melt_samples"][sample] = rs_post_align.OUTPUT_FORMAT.format(sample=sample)
        
    import yaml
    open(run_conf.project_root + "/config.yml", "w").write(yaml.dump(y))

