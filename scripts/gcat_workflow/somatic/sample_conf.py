#! /usr/bin/env python
import gcat_workflow.core.sample_conf_abc as abc

class Sample_conf(abc.Sample_conf_abc):
    SECTION_FASTQ = "fastq"
    SECTION_BAM_IMPORT = "bam-import"
    SECTION_BAM_TOFASTQ = "bam-tofastq"
    SECTION_MTCALL = "mutectcaller-parabricks"
    SECTION_WGS_METRICS = "collect-wgs-metrics"
    SECTION_MULTIPLE_METRICS = "collect-multiple-metrics"
    SECTION_GRIDSS = "gridss"
    SECTION_MANTA = "manta"
    SECTION_GENOMON_SV = "genomon-sv"
    SECTION_CONTROL_PANEL = "controlpanel"
    
    def __init__(self, sample_conf_file, exist_check = True):

        self.fastq = {}
        self.fastq_src = {}
        self.bam_tofastq = {}
        self.bam_tofastq_src = {}
        self.bam_import = {}
        self.bam_import_src = {}
        self.control_panel = {}
        self.mutect_call = []
        self.wgs_metrics = []
        self.multiple_metrics = []
        self.gridss = []
        self.manta = []
        self.genomon_sv = []
        self.exist_check = exist_check
        
        self.parse_file(sample_conf_file)

    def parse_data(self, _data):
        
        input_sections = [self.SECTION_FASTQ, self.SECTION_BAM_IMPORT, self.SECTION_BAM_TOFASTQ]
        analysis_sections = [
            self.SECTION_MTCALL, 
            self.SECTION_WGS_METRICS, 
            self.SECTION_MULTIPLE_METRICS, 
            self.SECTION_GRIDSS, 
            self.SECTION_MANTA, 
            self.SECTION_GENOMON_SV
        ]
        controlpanel_sections = [self.SECTION_CONTROL_PANEL]
        splited = self.split_section_data(_data, input_sections, analysis_sections, controlpanel_sections)
        
        sample_ids = []
        if self.SECTION_FASTQ in splited:
            parsed_fastq = self.parse_data_fastq_pair(splited[self.SECTION_FASTQ])
            self.fastq.update(parsed_fastq["fastq"])
            self.fastq_src.update(parsed_fastq["fastq_src"])
            sample_ids += parsed_fastq["fastq"].keys()
            
        if self.SECTION_BAM_TOFASTQ in splited:
            parsed_bam_tofastq = self.parse_data_bam_tofastq(splited[self.SECTION_BAM_TOFASTQ])
            self.bam_tofastq.update(parsed_bam_tofastq["bam_tofastq"])
            self.bam_tofastq_src.update(parsed_bam_tofastq["bam_tofastq_src"])
            sample_ids += parsed_bam_tofastq["bam_tofastq"].keys()
            
        if self.SECTION_BAM_IMPORT in splited:
            parsed_bam_import = self.parse_data_bam_import(splited[self.SECTION_BAM_IMPORT])
            self.bam_import.update(parsed_bam_import["bam_import"])
            self.bam_import_src.update(parsed_bam_import["bam_import_src"])
            sample_ids += parsed_bam_import["bam_import"].keys()
            
        if self.SECTION_MTCALL in splited:
            self.mutect_call += self.parse_data_tumor_normal(splited[self.SECTION_MTCALL], sample_ids, self.SECTION_MTCALL)
        
        if self.SECTION_WGS_METRICS in splited:
            self.wgs_metrics += self.parse_data_general(splited[self.SECTION_WGS_METRICS])

        if self.SECTION_MULTIPLE_METRICS in splited:
            self.multiple_metrics += self.parse_data_general(splited[self.SECTION_MULTIPLE_METRICS])

        if self.SECTION_GRIDSS in splited:
            self.gridss += self.parse_data_general(splited[self.SECTION_GRIDSS])

        if self.SECTION_MANTA in splited:
            self.manta += self.parse_data_general(splited[self.SECTION_MANTA])
        
        if self.SECTION_CONTROL_PANEL in splited:
            self.control_panel.update(self.parse_data_controlpanel(
                splited[self.SECTION_CONTROL_PANEL], 
                self.SECTION_CONTROL_PANEL
            ))
        
        if self.SECTION_GENOMON_SV in splited:
            self.genomon_sv += self.parse_data_tumor_normal_controlpanel(
                splited[self.SECTION_GENOMON_SV],
                sample_ids,
                self.control_panel.keys(),
                self.SECTION_GENOMON_SV
            )
