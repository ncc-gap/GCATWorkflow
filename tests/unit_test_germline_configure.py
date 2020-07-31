# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:18:34 2016

@author: okada

"""

import os
import sys
import shutil
import unittest
import subprocess
import snakemake

def func_path (root, name):
    wdir = root + "/" + name
    ss_path = root + "/" + name + ".csv"
    return (wdir, ss_path)

BAM_IMP = "bam_import"
BAM_2FQ = "bam_tofastq"
ALN = "bwa_alignment_parabricks"
HT_CALL = "haplotypecaller_parabricks"
SUMMARY1 = "collect_wgs_metrics"
SUMMARY2 = "collect_multiple_metrics"

class ConfigureTest(unittest.TestCase):
    
    DATA_DIR = "/tmp/temp-test/gcat_test_germline_configure"
    SAMPLE_DIR = DATA_DIR + "/samples"
    REMOVE = False
    SS_NAME = "/test.csv"
    GC_NAME = "/gcat.cfg"
    GC_NAME_P = "/gcat_parabricks.cfg"

    # init class
    @classmethod
    def setUpClass(self):
        if os.path.exists(self.DATA_DIR):
            shutil.rmtree(self.DATA_DIR)
        os.makedirs(self.SAMPLE_DIR, exist_ok = True)
        os.makedirs(self.DATA_DIR + "/reference", exist_ok = True)
        os.makedirs(self.DATA_DIR + "/image", exist_ok = True)
        os.makedirs(self.DATA_DIR + "/parabricks", exist_ok = True)
        touch_files = [
            "/samples/A1.fastq",
            "/samples/A2.fastq",
            "/samples/B1.fq",
            "/samples/B2.fq",
            "/samples/C1_1.fq",
            "/samples/C1_2.fq",
            "/samples/C2_1.fq",
            "/samples/C2_2.fq",
            "/samples/A.markdup.cram",
            "/samples/A.markdup.cram.crai",
            "/samples/B.markdup.cram",
            "/samples/B.markdup.crai",
            "/reference/XXX.fa",
            "/image/YYY.simg",
            "/parabricks/pbrun",
        ]
        for p in touch_files:
            open(self.DATA_DIR + p, "w").close()
        
        for p in ["/samples/A.metadata.txt", "/samples/B.metadata.txt", "/samples/C.metadata.txt"]:
            f = open(self.DATA_DIR + p, "w")
            f.write("@RG:xxx\n@RG:xxx\n@RG:xxx")
            f.close()

        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq
pool2,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq

[{bam2fq}]
A_control,{sample_dir}/A.markdup.cram

[{bamimp}]
pool3,{sample_dir}/B.markdup.cram

[{ht_call}]
A_tumor
A_control
pool1
pool2
pool3

[{summary1}]
A_tumor

[{summary2}]
A_tumor

[manta]
A_tumor

[melt]
A_tumor

[gridss]
A_tumor

[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
pool1,{sample_dir}/B.metadata.txt
pool2,{sample_dir}/C.metadata.txt
A_control,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ, bamimp = BAM_IMP, ht_call = HT_CALL, summary1 = SUMMARY1, summary2 = SUMMARY2)
        
        f = open(self.DATA_DIR + self.SS_NAME, "w")
        f.write(data_sample)
        f.close()

        conf_template = """[{bam2fq}]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg

[gatk_{aln}_compatible]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[{aln}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
reference = {sample_dir}/reference/XXX.fa

[post_{aln}]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[gatk_{ht_call}_compatible]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[{ht_call}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
reference = {sample_dir}/reference/XXX.fa

[gatk_{summary1}_compatible]
qsub_option = -l s_vmem=32G,mem_req=32G
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[{summary1}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
qsub_option = -l s_vmem=32G,mem_req=32G
reference = {sample_dir}/reference/XXX.fa

[gatk_{summary2}_compatible]
qsub_option = -l s_vmem=32G,mem_req=32G
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[{summary2}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
qsub_option = -l s_vmem=32G,mem_req=32G
reference = {sample_dir}/reference/XXX.fa

[gridss]
qsub_option = -l s_vmem=4G,mem_req=4G -pe def_slot 8
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[manta]
qsub_option = -l s_vmem=2G,mem_req=2G -pe def_slot 8
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[melt]
qsub_option = -l s_vmem=32G,mem_req=32G
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
"""
        # Not parabricks
        data_conf = conf_template.format(
            sample_dir = self.DATA_DIR, 
            bam2fq = BAM_2FQ, aln = ALN, ht_call = HT_CALL, summary1 = SUMMARY1, summary2 = SUMMARY2,
            gpu_support = "False"
        )
        
        f = open(self.DATA_DIR + self.GC_NAME, "w")
        f.write(data_conf)
        f.close()
        
        # parabricks
        data_conf2 = conf_template.format(
            sample_dir = self.DATA_DIR, 
            bam2fq = BAM_2FQ, aln = ALN, ht_call = HT_CALL, summary1 = SUMMARY1, summary2 = SUMMARY2,
            gpu_support = "True"
        )
        
        f = open(self.DATA_DIR + self.GC_NAME_P, "w")
        f.write(data_conf2)
        f.close()
        
    # terminated class
    @classmethod
    def tearDownClass(self):
        if self.REMOVE:
            shutil.rmtree(self.DATA_DIR)

    # init method
    def setUp(self):
        pass

    # terminated method
    def tearDown(self):
        pass
    
    def test1_01_version(self):
        subprocess.check_call(['python', 'gcat_workflow', '--version'])
    
    def test1_02_version(self):
        subprocess.check_call(['python', 'gcat_runner', '--version'])
    
    def test2_01_configure_drmaa_nogpu(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        options = [
            "germline",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test2_02_configure_drmaa_gpu(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        options = [
            "germline",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME_P,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test2_03_configure_qsub_nogpu(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        options = [
            "germline",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME,
            "--runner", "drmaa",
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)
    
    def test2_04_configure_qsub_gpu(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        options = [
            "germline",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME_P,
            "--runner", "drmaa",
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test2_05_configure_slurm_nogpu(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        options = [
            "germline",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME,
            "--runner", "slurm",
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)
    
    def test2_06_configure_slurm_gpu(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        options = [
            "germline",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME_P,
            "--runner", "slurm",
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_01_bwa_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "germline",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_02_1_bwa_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bam2fq}]
A_tumor,{sample_dir}/A.markdup.cram
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "germline",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_03_bwa_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "germline",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_01_htc_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[{ht_call}]
A_tumor
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, ht_call = HT_CALL)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "germline",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_02_htc_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bam2fq}]
A_tumor,{sample_dir}/A.markdup.cram
[{ht_call}]
A_tumor
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ, ht_call = HT_CALL)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "germline",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_03_htc_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[{ht_call}]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP, ht_call = HT_CALL)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "germline",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


if __name__ == '__main__':
    unittest.main()
