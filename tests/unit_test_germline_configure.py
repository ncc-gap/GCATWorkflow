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
from tests import snakemake_compat as snakemake

def func_path (root, name):
    wdir = root + "/" + name
    ss_path = root + "/" + name + ".csv"
    conf_path = root + "/" + name + ".cfg"
    return (wdir, ss_path, conf_path)

BAM_IMP = "bam_import"
BAM_2FQ = "bam_tofastq"
ALN = "bwa_alignment_parabricks"
HT_CALL = "haplotypecaller_parabricks"
SUMMARY1 = "collect_wgs_metrics"
SUMMARY2 = "collect_multiple_metrics"
SUMMARY3 = "collect_hs_metrics"

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
        os.makedirs(self.DATA_DIR + "/tools", exist_ok = True)
        os.makedirs(self.DATA_DIR + "/glimpse", exist_ok = True)
        touch_files = [
            "/samples/A1.fastq",
            "/samples/A2.fastq",
            "/samples/B1.fq",
            "/samples/B2.fq",
            "/samples/C1_1.fq",
            "/samples/C1_2.fq",
            "/samples/C2_1.fq",
            "/samples/C2_2.fq",
            "/samples/D1_1.fq",
            "/samples/D1_2.fq",
            "/samples/D2_1.fq",
            "/samples/D2_2.fq",
            "/samples/D3_1.fq",
            "/samples/D3_2.fq",
            "/samples/D4_1.fq",
            "/samples/D4_2.fq",
            "/samples/D5_1.fq",
            "/samples/D5_2.fq",
            "/samples/D6_1.fq",
            "/samples/D6_2.fq",
            "/samples/D7_1.fq",
            "/samples/D7_2.fq",
            "/samples/D8_1.fq",
            "/samples/D8_2.fq",
            "/samples/D9_1.fq",
            "/samples/D9_2.fq",
            "/samples/D10_1.fq",
            "/samples/D10_2.fq",
            "/samples/D11_1.fq",
            "/samples/D11_2.fq",
            "/samples/D12_1.fq",
            "/samples/D12_2.fq",
            "/samples/D13_1.fq",
            "/samples/D13_2.fq",
            "/samples/D14_1.fq",
            "/samples/D14_2.fq",
            "/samples/A.markdup.cram",
            "/samples/A.markdup.cram.crai",
            "/samples/B.markdup.cram",
            "/samples/B.markdup.crai",
            "/reference/XXX.fa",
            "/reference/XXX.bed",
            "/image/YYY.simg",
            "/parabricks/pbrun",
            "/tools/bgzip",
            "/tools/tabix",
            "/glimpse/chr1.map.gz",
            "/glimpse/chr2.map.gz",
            "/glimpse/chr3.map.gz",
            "/glimpse/chr4.map.gz",
            "/glimpse/chr5.map.gz",
            "/glimpse/chr6.map.gz",
            "/glimpse/chr7.map.gz",
            "/glimpse/chr8.map.gz",
            "/glimpse/chr9.map.gz",
            "/glimpse/chr10.map.gz",
            "/glimpse/chr11.map.gz",
            "/glimpse/chr12.map.gz",
            "/glimpse/chr13.map.gz",
            "/glimpse/chr14.map.gz",
            "/glimpse/chr15.map.gz",
            "/glimpse/chr16.map.gz",
            "/glimpse/chr17.map.gz",
            "/glimpse/chr18.map.gz",
            "/glimpse/chr19.map.gz",
            "/glimpse/chr20.map.gz",
            "/glimpse/chr21.map.gz",
            "/glimpse/chr22.map.gz",
        ]
        for p in touch_files:
            open(self.DATA_DIR + p, "w").close()
        
        link_files = [
            ("/C1_1.fq", "/link_C1_1.fq"),
            ("/C1_2.fq", "/link_C1_2.fq"),
            ("/A.markdup.cram", "/link_A.markdup.cram"),
            ("/A.markdup.crai", "/link_A.markdup.cram.crai"),
            ("/B.markdup.cram", "/link_B.markdup.cram"),
            ("/B.markdup.crai", "/link_B.markdup.crai"),
        ]
        for s in link_files:
            os.symlink(self.SAMPLE_DIR + s[0], self.SAMPLE_DIR + s[1])
        
        for p in [("/samples/A.metadata.txt", "A", 1), ("/samples/B.metadata.txt", "B", 1), ("/samples/C.metadata.txt", "C", 2), ("/samples/D.metadata.txt", "D", 14)]:
            f = open(self.DATA_DIR + p[0], "w")
            for i in range(p[2]):
                f.write("@RG:%s_%d\n" % (p[1], i))
            f.close()

        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq
pool2,{sample_dir}/link_C1_1.fq;{sample_dir}/link_C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq
D_tumor,{sample_dir}/D1_1.fq;{sample_dir}/D2_1.fq;{sample_dir}/D3_1.fq;{sample_dir}/D4_1.fq;{sample_dir}/D5_1.fq;{sample_dir}/D6_1.fq;{sample_dir}/D7_1.fq;{sample_dir}/D8_1.fq;{sample_dir}/D9_1.fq;{sample_dir}/D10_1.fq;{sample_dir}/D11_1.fq;{sample_dir}/D12_1.fq;{sample_dir}/D13_1.fq;{sample_dir}/D14_1.fq,{sample_dir}/D1_2.fq;{sample_dir}/D2_2.fq;{sample_dir}/D3_2.fq;{sample_dir}/D4_2.fq;{sample_dir}/D5_2.fq;{sample_dir}/D6_2.fq;{sample_dir}/D7_2.fq;{sample_dir}/D8_2.fq;{sample_dir}/D9_2.fq;{sample_dir}/D10_2.fq;{sample_dir}/D11_2.fq;{sample_dir}/D12_2.fq;{sample_dir}/D13_2.fq;{sample_dir}/D14_2.fq

[{bam2fq}]
A_control,{sample_dir}/A.markdup.cram
A_control2,{sample_dir}/link_A.markdup.cram

[{bamimp}]
pool3,{sample_dir}/A.markdup.cram
pool4,{sample_dir}/link_B.markdup.cram

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

[{summary3}]
A_tumor

[manta]
A_tumor

[melt]
A_tumor

[gridss]
A_tumor

[glimpse]
A_tumor

[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
pool1,{sample_dir}/B.metadata.txt
pool2,{sample_dir}/C.metadata.txt
A_control,{sample_dir}/A.metadata.txt
A_control2,{sample_dir}/A.metadata.txt
D_tumor,{sample_dir}/D.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ, bamimp = BAM_IMP, ht_call = HT_CALL, summary1 = SUMMARY1, summary2 = SUMMARY2, summary3 = SUMMARY3)
        
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
gatk_recal = True

[{aln}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
reference = {sample_dir}/reference/XXX.fa
fq2bam_markdup_metrics = True
fq2bam_recal = True

[post_{aln}]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[gatk_{ht_call}_compatible]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
interval_autosome = {sample_dir}/reference/XXX.bed
interval_par = {sample_dir}/reference/XXX.bed
interval_chrx = {sample_dir}/reference/XXX.bed
interval_chry = {sample_dir}/reference/XXX.bed

[{ht_call}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
reference = {sample_dir}/reference/XXX.fa
interval_autosome = {sample_dir}/reference/XXX.bed
interval_par = {sample_dir}/reference/XXX.bed
interval_chrx = {sample_dir}/reference/XXX.bed
interval_chry = {sample_dir}/reference/XXX.bed
bgzip = {sample_dir}/tools/bgzip
tabix = {sample_dir}/tools/tabix

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

[{summary3}]
qsub_option = -l s_vmem=32G,mem_req=32G
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
bait_intervals = {sample_dir}/reference/XXX.bed
target_intervals = {sample_dir}/reference/XXX.bed

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

[glimpse]
qsub_option = -l s_vmem=4G,mem_req=4G -pe def_slot 8
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
map_dir = {sample_dir}/glimpse
glimpse_reference_prefix ={sample_dir}
"""
        # Not parabricks
        data_conf = conf_template.format(
            sample_dir = self.DATA_DIR, 
            bam2fq = BAM_2FQ, aln = ALN, ht_call = HT_CALL, summary1 = SUMMARY1, summary2 = SUMMARY2, summary3 = SUMMARY3,
            gpu_support = "False"
        )
        
        f = open(self.DATA_DIR + self.GC_NAME, "w")
        f.write(data_conf)
        f.close()
        
        # parabricks
        data_conf2 = conf_template.format(
            sample_dir = self.DATA_DIR, 
            bam2fq = BAM_2FQ, aln = ALN, ht_call = HT_CALL, summary1 = SUMMARY1, summary2 = SUMMARY2, summary3 = SUMMARY3,
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
        (wdir, _, _) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
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
        (wdir, _, _) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
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
        (wdir, _, _) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
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
        (wdir, _, _) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
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
        (wdir, _, _) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
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
        (wdir, _, _) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
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

    def test2_21_dag(self):
        (wdir, ss_path, _) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)

        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq

[{ht_call}]
A_tumor

[{summary1}]
A_tumor

[{summary2}]
A_tumor

[manta]
A_tumor

[gridss]
A_tumor

[glimpse]
A_tumor

[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ, bamimp = BAM_IMP, ht_call = HT_CALL, summary1 = SUMMARY1, summary2 = SUMMARY2)
        
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
        sys.stdout = open(wdir + '/germline.dot', 'w')
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True, printdag = True)
        sys.stdout.close()
        sys.stdout = sys.__stdout__
        self.assertTrue(success)
        if subprocess.call('which dot', shell=True) == 0:
            subprocess.check_call('dot -Tpng {wdir}/germline.dot > {wdir}/dag_germline.png'.format(wdir=wdir), shell=True)

    def test3_01_bwa_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR)
        
        data_conf = """[gatk_{aln}_compatible]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
gatk_recal = True

[post_{aln}]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR, aln = ALN)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_02_bwa_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR)
        
        data_conf = """[{aln}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
reference = {sample_dir}/reference/XXX.fa
fq2bam_markdup_metrics = True
fq2bam_recal = True

[post_{aln}]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR, bam2fq = BAM_2FQ, aln = ALN, gpu_support = "True")

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_03_bwa_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bam2fq}]
A_tumor,{sample_dir}/A.markdup.cram
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ)

        data_conf = """[{bam2fq}]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg

[gatk_{aln}_compatible]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
gatk_recal = True

[post_{aln}]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR, bam2fq = BAM_2FQ, aln = ALN)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test3_04_bwa_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bam2fq}]
A_tumor,{sample_dir}/A.markdup.cram
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ)

        data_conf = """[{bam2fq}]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg

[{aln}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
reference = {sample_dir}/reference/XXX.fa
fq2bam_markdup_metrics = True
fq2bam_recal = True

[post_{aln}]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR, bam2fq = BAM_2FQ, aln = ALN, gpu_support = "True")

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_05_bwa_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP)

        data_conf = ""

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_01_htc_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[{ht_call}]
A_tumor
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, ht_call = HT_CALL)
        
        data_conf = """[gatk_{aln}_compatible]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
gatk_recal = True

[post_{aln}]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[gatk_{ht_call}_compatible]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
interval_autosome = {sample_dir}/reference/XXX.bed
interval_par = {sample_dir}/reference/XXX.bed
interval_chrx = {sample_dir}/reference/XXX.bed
interval_chry = {sample_dir}/reference/XXX.bed
""".format(sample_dir = self.DATA_DIR, bam2fq = BAM_2FQ, aln = ALN, ht_call = HT_CALL)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_02_htc_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[{ht_call}]
A_tumor
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, ht_call = HT_CALL)
        
        data_conf = """[{aln}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
reference = {sample_dir}/reference/XXX.fa
fq2bam_markdup_metrics = True
fq2bam_recal = True

[post_{aln}]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[{ht_call}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
reference = {sample_dir}/reference/XXX.fa
interval_autosome = {sample_dir}/reference/XXX.bed
interval_par = {sample_dir}/reference/XXX.bed
interval_chrx = {sample_dir}/reference/XXX.bed
interval_chry = {sample_dir}/reference/XXX.bed
bgzip = {sample_dir}/tools/bgzip
tabix = {sample_dir}/tools/tabix
""".format(sample_dir = self.DATA_DIR, bam2fq = BAM_2FQ, aln = ALN, ht_call = HT_CALL, gpu_support = "True")

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_03_htc_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bam2fq}]
A_tumor,{sample_dir}/A.markdup.cram
[{ht_call}]
A_tumor
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ, ht_call = HT_CALL)
        
        data_conf = """[{bam2fq}]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg

[gatk_{aln}_compatible]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
gatk_recal = True

[post_{aln}]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[gatk_{ht_call}_compatible]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
interval_autosome = {sample_dir}/reference/XXX.bed
interval_par = {sample_dir}/reference/XXX.bed
interval_chrx = {sample_dir}/reference/XXX.bed
interval_chry = {sample_dir}/reference/XXX.bed
""".format(sample_dir = self.DATA_DIR, bam2fq = BAM_2FQ, aln = ALN, ht_call = HT_CALL)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test4_04_htc_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bam2fq}]
A_tumor,{sample_dir}/A.markdup.cram
[{ht_call}]
A_tumor
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ, ht_call = HT_CALL)
        
        data_conf = """[{bam2fq}]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg

[{aln}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
reference = {sample_dir}/reference/XXX.fa
fq2bam_markdup_metrics = True
fq2bam_recal = True

[post_{aln}]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[{ht_call}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
reference = {sample_dir}/reference/XXX.fa
interval_autosome = {sample_dir}/reference/XXX.bed
interval_par = {sample_dir}/reference/XXX.bed
interval_chrx = {sample_dir}/reference/XXX.bed
interval_chry = {sample_dir}/reference/XXX.bed
bgzip = {sample_dir}/tools/bgzip
tabix = {sample_dir}/tools/tabix
""".format(sample_dir = self.DATA_DIR, bam2fq = BAM_2FQ, aln = ALN, ht_call = HT_CALL, gpu_support = "True")

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_05_htc_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[{ht_call}]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP, ht_call = HT_CALL)
        
        data_conf = """[gatk_{ht_call}_compatible]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
interval_autosome = {sample_dir}/reference/XXX.bed
interval_par = {sample_dir}/reference/XXX.bed
interval_chrx = {sample_dir}/reference/XXX.bed
interval_chry = {sample_dir}/reference/XXX.bed
""".format(sample_dir = self.DATA_DIR, ht_call = HT_CALL)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test4_06_htc_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[{ht_call}]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP, ht_call = HT_CALL)
        
        data_conf = """[{ht_call}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
reference = {sample_dir}/reference/XXX.fa
interval_autosome = {sample_dir}/reference/XXX.bed
interval_par = {sample_dir}/reference/XXX.bed
interval_chrx = {sample_dir}/reference/XXX.bed
interval_chry = {sample_dir}/reference/XXX.bed
bgzip = {sample_dir}/tools/bgzip
tabix = {sample_dir}/tools/tabix
""".format(sample_dir = self.DATA_DIR, ht_call = HT_CALL, gpu_support = "True")

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_01_metrics_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[{summary1}]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP, summary1 = SUMMARY1)
        
        data_conf = """[gatk_{summary1}_compatible]
qsub_option = -l s_vmem=32G,mem_req=32G
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR, summary1 = SUMMARY1)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_02_metrics_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[{summary1}]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP, summary1 = SUMMARY1)
        
        data_conf = """[{summary1}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
qsub_option = -l s_vmem=32G,mem_req=32G
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR, summary1 = SUMMARY1, gpu_support = "True")

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_03_metrics_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[{summary2}]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP, summary2 = SUMMARY2)
        
        data_conf = """[gatk_{summary2}_compatible]
qsub_option = -l s_vmem=32G,mem_req=32G
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR, summary2 = SUMMARY2)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test5_04_metrics_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[{summary2}]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP, summary2 = SUMMARY2)
        
        data_conf = """[{summary2}]
gpu_support = {gpu_support}
pbrun = {sample_dir}/parabricks/pbrun
qsub_option = -l s_vmem=32G,mem_req=32G
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR, summary2 = SUMMARY2, gpu_support = "True")

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_05_metrics_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[{summary3}]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP, summary3 = SUMMARY3)
        
        data_conf = """[{summary3}]
qsub_option = -l s_vmem=32G,mem_req=32G
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
bait_intervals = {sample_dir}/reference/XXX.bed
target_intervals = {sample_dir}/reference/XXX.bed
""".format(sample_dir = self.DATA_DIR, summary3 = SUMMARY3)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test6_01_manta_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[manta]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP)
        
        data_conf = """[manta]
qsub_option = -l s_vmem=2G,mem_req=2G -pe def_slot 8
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test7_01_gridss_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[gridss]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP)
        
        data_conf = """[gridss]
qsub_option = -l s_vmem=2G,mem_req=2G -pe def_slot 8
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test8_01_melt_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[melt]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP)
        
        data_conf = """[melt]
qsub_option = -l s_vmem=2G,mem_req=2G -pe def_slot 8
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test9_01_glimpse_limited(self):
        (wdir, ss_path, conf_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram
[glimpse]
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP)
        
        data_conf = """[glimpse]
qsub_option = -l s_vmem=4G,mem_req=4G -pe def_slot 8
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
map_dir = {sample_dir}/glimpse
glimpse_reference_prefix ={sample_dir}/
""".format(sample_dir = self.DATA_DIR)

        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()

        f = open(conf_path, "w")
        f.write(data_conf)
        f.close()

        options = [
            "germline",
            ss_path,
            wdir,
            conf_path,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

if __name__ == '__main__':
    unittest.main()
