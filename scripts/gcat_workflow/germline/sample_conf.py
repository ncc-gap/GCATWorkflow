#! /usr/bin/env python
import gcat_workflow.core.sample_conf_abc as abc

class Sample_conf(abc.Sample_conf_abc):
    SECTION_FASTQ = "fastq"
    SECTION_BAM_IMPORT = "bam_import"
    SECTION_BAM_TOFASTQ = "bam_tofastq"
    SECTION_HTCALL = "haplotypecaller_parabricks"
    SECTION_WGS_METRICS = "collect_wgs_metrics"
    SECTION_HS_METRICS = "collect_hs_metrics"
    SECTION_MULTIPLE_METRICS = "collect_multiple_metrics"
    SECTION_GRIDSS = "gridss"
    SECTION_MANTA = "manta"
    SECTION_MELT = "melt"
    SECTION_READGROUP = "readgroup"
    
    def __init__(self, sample_conf_file, exist_check = True):

        self.fastq = {}
        self.fastq_src = {}
        self.bam_tofastq = {}
        self.bam_tofastq_src = {}
        self.bam_import = {}
        self.bam_import_src = {}
        self.haplotype_call = []
        self.wgs_metrics = []
        self.hs_metrics = []
        self.multiple_metrics = []
        self.gridss = []
        self.manta = []
        self.melt = []
        self.readgroup = {}
        self.readgroup_src = {}
        self.exist_check = exist_check
        
        self.parse_file(sample_conf_file)

    def parse_data(self, _data):
        
        input_sections = [self.SECTION_FASTQ, self.SECTION_BAM_IMPORT, self.SECTION_BAM_TOFASTQ]
        analysis_sections = [self.SECTION_HTCALL, self.SECTION_WGS_METRICS, self.SECTION_HS_METRICS, self.SECTION_MULTIPLE_METRICS, self.SECTION_GRIDSS, self.SECTION_MANTA, self.SECTION_MELT]
        controlpanel_sections = []
        extend_sections = [self.SECTION_READGROUP]
        splited = self.split_section_data(_data, input_sections, analysis_sections, controlpanel_sections, extend_sections)
        
        bwa_samples = []
        if self.SECTION_FASTQ in splited:
            parsed_fastq = self.parse_data_fastq_pair(splited[self.SECTION_FASTQ])
            self.fastq.update(parsed_fastq["fastq"])
            self.fastq_src.update(parsed_fastq["fastq_src"])
            bwa_samples += self.fastq.keys()
            
        if self.SECTION_BAM_TOFASTQ in splited:
            parsed_bam_tofastq = self.parse_data_bam_tofastq(splited[self.SECTION_BAM_TOFASTQ])
            self.bam_tofastq.update(parsed_bam_tofastq["bam_tofastq"])
            self.bam_tofastq_src.update(parsed_bam_tofastq["bam_tofastq_src"])
            bwa_samples += self.bam_tofastq.keys()
            
        if self.SECTION_BAM_IMPORT in splited:
            parsed_bam_import = self.parse_data_bam_import(splited[self.SECTION_BAM_IMPORT])
            self.bam_import.update(parsed_bam_import["bam_import"])
            self.bam_import_src.update(parsed_bam_import["bam_import_src"])
            
        if self.SECTION_HTCALL in splited:
            self.haplotype_call += self.parse_data_general(splited[self.SECTION_HTCALL])

        if self.SECTION_WGS_METRICS in splited:
            self.wgs_metrics += self.parse_data_general(splited[self.SECTION_WGS_METRICS])

        if self.SECTION_HS_METRICS in splited:
            self.hs_metrics += self.parse_data_general(splited[self.SECTION_HS_METRICS])

        if self.SECTION_MULTIPLE_METRICS in splited:
            self.multiple_metrics += self.parse_data_general(splited[self.SECTION_MULTIPLE_METRICS])

        if self.SECTION_GRIDSS in splited:
            self.gridss += self.parse_data_general(splited[self.SECTION_GRIDSS])

        if self.SECTION_MANTA in splited:
            self.manta += self.parse_data_general(splited[self.SECTION_MANTA])
        
        if self.SECTION_MELT in splited:
            self.melt += self.parse_data_general(splited[self.SECTION_MELT])
        
        if len(bwa_samples) > 0:
            if not self.SECTION_READGROUP in splited:
                err_msg = "[%s] section is not set" % (self.SECTION_READGROUP)
                raise ValueError(err_msg)
            parsed_bam_tofastq = self.parse_data_readgroup(splited[self.SECTION_READGROUP], bwa_samples, self.SECTION_READGROUP)
            self.readgroup.update(parsed_bam_tofastq["metadata"])
            self.readgroup_src.update(parsed_bam_tofastq["metadata_src"])
                
