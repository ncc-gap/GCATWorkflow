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

class ConfigureTest(unittest.TestCase):
    
    DATA_DIR = "/tmp/temp-test/gcat_test_rna_configure"
    SAMPLE_DIR = DATA_DIR + "/samples"
    REMOVE = False
    SS_NAME = "/test.csv"
    GC_NAME = "/gcat.cfg"

    # init class
    @classmethod
    def setUpClass(self):
        shutil.rmtree(self.DATA_DIR)
        os.makedirs(self.SAMPLE_DIR, exist_ok = True)
        os.makedirs(self.DATA_DIR + "/reference", exist_ok = True)
        os.makedirs(self.DATA_DIR + "/image", exist_ok = True)
        touch_files = [
            "/samples/A1.fastq",
            "/samples/A2.fastq",
            "/samples/B1.fq",
            "/samples/B2.fq",
            "/samples/C1_1.fq",
            "/samples/C1_2.fq",
            "/samples/C2_1.fq",
            "/samples/C2_2.fq",
            "/samples/D1.fastq",
            "/samples/E1.fq",
            "/samples/F1.fq",
            "/samples/F2.fq",
            "/samples/G.Aligned.sortedByCoord.out.bam",
            "/samples/G.Aligned.sortedByCoord.out.bam.bai",
            "/samples/H.Aligned.sortedByCoord.out.bam",
            "/samples/H.Aligned.sortedByCoord.out.bai",
            "/samples/I.Aligned.sortedByCoord.out.bam",
            "/samples/I.Aligned.sortedByCoord.out.bam.bai",
            "/samples/I.Aligned.sortedByCoord.out.bam.bai",
            "/samples/I.Chimeric.out.junction",
            "/samples/I.Chimeric.out.sam",
            "/reference/XXX.fa",
            "/image/YYY.simg",
            "/reference/ZZZ.gtf"
        ]
        for p in touch_files:
            open(self.DATA_DIR + p, "w").close()
        
        data_sample = """[fastq],,,,
sampleA,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
sampleB,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,
sampleC,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq,,
sampleD,{sample_dir}/D1.fastq
sampleE,{sample_dir}/E1.fq
sampleF,{sample_dir}/F1.fq;{sample_dir}/F2.fq,,,,

[bam_tofastq_pair],,,,
sampleG,{sample_dir}/G.Aligned.sortedByCoord.out.bam,,,

[bam_tofastq_single],,,,
sampleH,{sample_dir}/H.Aligned.sortedByCoord.out.bam,,,

[bam_import],,,,
sampleI,{sample_dir}/I.Aligned.sortedByCoord.out.bam,,,
,,,,
[fusionfusion],,,,
sampleA,list1
sampleD,list1
sampleG,None
sampleH,,,,
sampleI,,,,
,,,,
[expression],,,,
sampleA,,,,
sampleD,,,,
,,,,
[qc],,,,
sampleA,,,,
sampleB,
sampleC,
sampleD,,,,
sampleE,
sampleF
sampleG,,,,
sampleH,,,,
sampleI,,,,
,,,,
[star_fusion],,,,
sampleA,,,,
sampleD,,,,
sampleG,,,,
sampleH,,,,
sampleI,,,,
,,,,
[iravnet],,,,
sampleA,,,,
sampleD,,,,
sampleG,,,,
sampleH,,,,
sampleI,,,,
,,,,
[intron_retention],,,,
sampleA,,,,
sampleD,,,,
sampleG,,,,
sampleH,,,,
sampleI,,,,
,,,,
[kalisto],,,,
sampleA,,,,
sampleD,,,,
sampleG,,,,
sampleH,,,,
sampleI,,,,
,,,,
[controlpanel]
list1,sampleB,sampleC,sampleE,sampleF
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(self.DATA_DIR + self.SS_NAME, "w")
        f.write(data_sample)
        f.close()

        data_conf = """[bam_tofastq]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 
params = collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY tryoq=0

[star_alignment]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 
star_option = --runThreadN 6
star_genome = {sample_dir}/reference
remove_fastq = True

[star_fusion]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 
star_fusion_option = --runThreadN 6
star_genome = {sample_dir}/reference

[fusionfusion]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 

[expression]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 
gtf = {sample_dir}/reference/ZZZ.gtf

[intron_retention]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 

[iravnet]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 

[kalisto]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
singularity_option = 

""".format(sample_dir = self.DATA_DIR)
        
        f = open(self.DATA_DIR + self.GC_NAME, "w")
        f.write(data_conf)
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
    
    def test2_01_configure(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        options = [
            "rna",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test2_02_configure(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)

        options = [
            "rna",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME,
            "--disable_drmaa",
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_01_star_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_02_1_star_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_single]
sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_02_2_star_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_pair]
sample,{sample_dir}/G.Aligned.sortedByCoord.out.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test3_03_star_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        

        data_sample = """[bam_import]
sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test4_01_fusionfusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[fusionfusion]
sample,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_02_fusionfusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_single]
sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
[fusionfusion]
sample,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_03_fusionfusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[fusionfusion]
sample,None
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test5_01_expression_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[expression]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_02_expression_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_single]
sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
[expression]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_03_expression_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[expression]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test6_01_star_fusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[star_fusion]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test6_02_star_fusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_single]
sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
[star_fusion]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test6_03_star_fusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[star_fusion]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test7_01_iravnet_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[iravnet]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test7_02_iravnet_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_single]
sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
[iravnet]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test7_03_iravnet_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[iravnet]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test8_01_ir_count_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[intron_retention]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test8_02_ir_count_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_single]
sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
[intron_retention]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test8_03_ir_count_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[intron_retention]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test9_01_kalisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[kalisto]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test9_02_kalisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq_single]
sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
[kalisto]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test9_03_kalisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[kalisto]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

if __name__ == '__main__':
    unittest.main()
