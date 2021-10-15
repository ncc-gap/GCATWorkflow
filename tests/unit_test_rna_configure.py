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
    GC2_NAME = "/gcat2.cfg"

    # init class
    @classmethod
    def setUpClass(self):
        if os.path.exists(self.DATA_DIR):
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
            "/samples/I.Chimeric.out.sam",
            "/samples/I.SJ.out.tab.gz",
            "/samples/L.Aligned.sortedByCoord.out.bam",
            "/samples/L.Aligned.sortedByCoord.out.bai",
            "/samples/L.Chimeric.out.sam",
            "/samples/L.SJ.out.tab.gz",
            "/samples/M.Aligned.sortedByCoord.out.cram",
            "/samples/M.Aligned.sortedByCoord.out.cram.crai",
            "/samples/M.Chimeric.out.sam",
            "/samples/M.SJ.out.tab.gz",
            "/samples/N.Aligned.sortedByCoord.out.cram",
            "/samples/N.Aligned.sortedByCoord.out.crai",
            "/samples/N.Chimeric.out.sam",
            "/samples/N.SJ.out.tab.gz",
            "/samples/run.sra",
            "/reference/XXX.fa",
            "/image/YYY.simg",
            "/reference/ZZZ.gtf",
            "/reference/ZZZ.vcf.gz",
            "/reference/ZZZ.bed",
            "/reference/ZZZ.idx",
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

[bam_tofastq],,,,
sampleG,{sample_dir}/G.Aligned.sortedByCoord.out.bam,,,
sampleH,{sample_dir}/H.Aligned.sortedByCoord.out.bam,,,

[bam_import],,,,
sampleI,{sample_dir}/I.Aligned.sortedByCoord.out.bam,,,
sampleL,{sample_dir}/L.Aligned.sortedByCoord.out.bam,,,
sampleM,{sample_dir}/M.Aligned.sortedByCoord.out.cram,,,
sampleN,{sample_dir}/N.Aligned.sortedByCoord.out.cram,,,

[sra_fastq_dump],,,,
sampleJ,RUNID123456
sampleK,RUNID123456,http://dummy.com/data/run.sra
sampleO,RUNID123456,{sample_dir}/run.sra

[fusionfusion],,,,
sampleA,list1
sampleD,list1
sampleG,None
sampleH,,,,
sampleI,,,,
sampleJ
sampleK
sampleL
sampleM
sampleN
sampleO
,,,,,,,,,,,,
[expression],,,,
sampleA,,,,
sampleD,,,,
sampleJ
sampleK
sampleL
sampleM
sampleN
sampleO
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
sampleJ
sampleK
sampleL
sampleM
sampleN
sampleO
[iravnet],,,,
sampleA,,,,
sampleB,,,,
sampleC,,,,
sampleD,,,,
sampleE,,,,
sampleF,,,,
sampleG,,,,
sampleH,,,,
sampleI,,,,
sampleJ
sampleK
[juncmut],,,,
sampleA,,,,
sampleB,,,,
sampleC,,,,
sampleD,,,,
sampleE,,,,
sampleF,,,,
sampleG,,,,
sampleH,,,,
sampleI,,,,
sampleJ
sampleK
sampleL
sampleM
sampleN
sampleO
[intron_retention],,,,
sampleA,,,,
sampleD,,,,
sampleG,,,,
sampleH,,,,
sampleI,,,,
sampleJ
[kallisto],,,,
sampleA,,,,
sampleD,,,,
sampleG,,,,
sampleH,,,,
sampleI,,,,
sampleJ
sampleK
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

[sra_fastq_dump]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg

[cram_tobam]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[star_alignment]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
star_genome = {sample_dir}/reference
remove_fastq = True
[star_fusion]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
star_genome = {sample_dir}/reference

[fusionfusion_count]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg

[fusionfusion_merge]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg

[fusionfusion]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[expression]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
gtf = {sample_dir}/reference/ZZZ.gtf

[intron_retention]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg

[iravnet]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
clinvar_db = {sample_dir}/reference/ZZZ.vcf.gz

[juncmut]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
control_file1 = {sample_dir}/reference/ZZZ.vcf.gz
control_file2 = {sample_dir}/reference/ZZZ.vcf.gz
genecode_gene_file = {sample_dir}/reference/ZZZ.vcf.gz

[kallisto]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference_fasta = {sample_dir}/reference/XXX.fa
reference_kallisto_index = {sample_dir}/reference/ZZZ.idx
annotation_gtf = {sample_dir}/reference/ZZZ.gtf

[join]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
remove_bam = True
bam_tocram = True
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR)
        
        f = open(self.DATA_DIR + self.GC_NAME, "w")
        f.write(data_conf)
        f.close()

        data_conf2 = """[bam_tofastq]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg

[sra_fastq_dump]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg

[cram_tobam]
qsub_option = -l s_vmem=2G,mem_req=2G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[star_alignment]
qsub_option = -l s_vmem=10.6G,mem_req=10.6G -l os7
image = {sample_dir}/image/YYY.simg
star_genome = {sample_dir}/reference

[star_fusion]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
star_genome = {sample_dir}/reference

[fusionfusion_count]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg

[fusionfusion_merge]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg

[fusionfusion]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa

[expression]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
gtf = {sample_dir}/reference/ZZZ.gtf

[intron_retention]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg

[iravnet]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
clinvar_db = {sample_dir}/reference/ZZZ.vcf.gz
target_file = {sample_dir}/reference/ZZZ.bed

[juncmut]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference = {sample_dir}/reference/XXX.fa
control_file1 = {sample_dir}/reference/ZZZ.vcf.gz
control_file2 = {sample_dir}/reference/ZZZ.vcf.gz
genecode_gene_file = {sample_dir}/reference/ZZZ.vcf.gz

[kallisto]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
reference_fasta = {sample_dir}/reference/XXX.fa
reference_kallisto_index = {sample_dir}/reference/ZZZ.idx
annotation_gtf = {sample_dir}/reference/ZZZ.gtf

[join]
qsub_option = -l s_vmem=5.3G,mem_req=5.3G -l os7
image = {sample_dir}/image/YYY.simg
remove_bam = True
bam_tocram = True
reference = {sample_dir}/reference/XXX.fa
""".format(sample_dir = self.DATA_DIR)
        
        f = open(self.DATA_DIR + self.GC2_NAME, "w")
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
            "--runner", "drmaa",
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test2_03_configure(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)

        options = [
            "rna",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC_NAME,
            "--runner", "slurm",
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test2_11_configure(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        options = [
            "rna",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test2_12_configure(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)

        options = [
            "rna",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
            "--runner", "drmaa",
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test2_13_configure(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)

        options = [
            "rna",
            self.DATA_DIR + self.SS_NAME,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
            "--runner", "slurm",
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)
    
    def test2_21_dag(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)

        data_sample = """[fastq]
sampleA,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq
[fusionfusion]
sampleA,list1
[expression]
sampleA
[star_fusion]
#sampleA
[iravnet]
sampleA
[juncmut]
sampleA
[intron_retention]
sampleA
[kallisto]
sampleA
[controlpanel]
list1,pool1
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
        sys.stdout = open(wdir + '/rna_star.dot', 'w')
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True, printdag = True)
        sys.stdout.close()
        sys.stdout = sys.__stdout__
        self.assertTrue(success)
        if subprocess.call('which dot', shell=True) == 0:
            subprocess.check_call('dot -Tpng {wdir}/rna_star.dot > {wdir}/dag_rna.png'.format(wdir=wdir), shell=True)

    def test2_22_dag(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)

        data_sample = """[sra_fastq_dump]
sampleR,RUNID123456
[expression]
sampleR
[iravnet]
sampleR
[juncmut]
sampleR
[intron_retention]
sampleR
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
        sys.stdout = open(wdir + '/rna_sra.dot', 'w')
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True, printdag = True)
        sys.stdout.close()
        sys.stdout = sys.__stdout__
        self.assertTrue(success)
        if subprocess.call('which dot', shell=True) == 0:
            subprocess.check_call('dot -Tpng {wdir}/rna_sra.dot > {wdir}/dag_rna_sradump.png'.format(wdir=wdir), shell=True)


    def test2_23_dag(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)

        data_sample = """[bam_import]
sampleA,{sample_dir}/M.Aligned.sortedByCoord.out.cram
pool1,{sample_dir}/I.Aligned.sortedByCoord.out.bam
pool2,{sample_dir}/N.Aligned.sortedByCoord.out.cram
[fusionfusion]
sampleA,list1
[expression]
sampleA
[star_fusion]
#sampleA
[iravnet]
sampleA
[juncmut]
sampleA
[intron_retention]
sampleA
[controlpanel]
list1,pool1,pool2
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
        sys.stdout = open(wdir + '/rna_cram_tobam.dot', 'w')
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True, printdag = True)
        sys.stdout.close()
        sys.stdout = sys.__stdout__
        self.assertTrue(success)
        if subprocess.call('which dot', shell=True) == 0:
            subprocess.check_call('dot -Tpng {wdir}/rna_cram_tobam.dot > {wdir}/dag_rna_sradump.png'.format(wdir=wdir), shell=True)


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

    def test3_02_star_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
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

    def test3_04_star_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        

        data_sample = """[sra_fastq_dump],,,,
sample,RUNID123456
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
sampleA,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
sampleB,{sample_dir}/B1.fq,{sample_dir}/B2.fq
pool1,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq
[fusionfusion]
sampleA,list1
sampleB,None
[controlpanel]
list1,pool1
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
        
        data_sample = """[bam_tofastq]
sampleA,{sample_dir}/H.Aligned.sortedByCoord.out.bam
sampleB,{sample_dir}/H.Aligned.sortedByCoord.out.bam
pool1,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[fusionfusion]
sampleA,list1
sampleB,None
[controlpanel]
list1,pool1
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

    def test4_04_fusionfusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[sra_fastq_dump]
sample,RUNID123456
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

    def test4_11_fusionfusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sampleA,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
sampleB,{sample_dir}/B1.fq,{sample_dir}/B2.fq
pool1,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq
[fusionfusion]
sampleA,list1
sampleB,None
[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_12_fusionfusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
sampleA,{sample_dir}/H.Aligned.sortedByCoord.out.bam
sampleB,{sample_dir}/H.Aligned.sortedByCoord.out.bam
pool1,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[fusionfusion]
sampleA,list1
sampleB,None
[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_13_fusionfusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sampleA,{sample_dir}/M.Aligned.sortedByCoord.out.cram
sampleB,{sample_dir}/I.Aligned.sortedByCoord.out.bam
pool1,{sample_dir}/N.Aligned.sortedByCoord.out.cram
[fusionfusion]
sampleA,list1
sampleB,None
[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test4_14_fusionfusion_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[sra_fastq_dump]
sampleA,RUNID123456
sampleB,RUNID123789
pool1,RUNID123012
[fusionfusion]
sampleA,list1
sampleB,None
[controlpanel]
list1,pool1
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
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
        
        data_sample = """[bam_tofastq]
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

    def test5_04_expression_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[sra_fastq_dump],,,,
sample,RUNID123456
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

    def test5_11_expression_limited(self):
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_12_expression_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_13_expression_limited(self):
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test5_14_expression_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[sra_fastq_dump],,,,
sample,RUNID123456
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

#    def test6_01_star_fusion_limited(self):
#        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
#        
#        data_sample = """[fastq]
#sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
#[star_fusion]
#sample
#""".format(sample_dir = self.SAMPLE_DIR)
#        
#        f = open(ss_path, "w")
#        f.write(data_sample)
#        f.close()
#        options = [
#            "rna",
#            ss_path,
#            wdir,
#            self.DATA_DIR + self.GC_NAME,
#        ]
#        subprocess.check_call(['python', 'gcat_workflow'] + options)
#        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
#        self.assertTrue(success)
#
#    def test6_02_star_fusion_limited(self):
#        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
#        
#        data_sample = """[bam_tofastq]
#sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
#[star_fusion]
#sample
#""".format(sample_dir = self.SAMPLE_DIR)
#        
#        f = open(ss_path, "w")
#        f.write(data_sample)
#        f.close()
#        options = [
#            "rna",
#            ss_path,
#            wdir,
#            self.DATA_DIR + self.GC_NAME,
#        ]
#        subprocess.check_call(['python', 'gcat_workflow'] + options)
#        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
#        self.assertTrue(success)
#
#    def test6_03_star_fusion_limited(self):
#        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
#        
#        data_sample = """[bam_import]
#sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
#[star_fusion]
#sample
#""".format(sample_dir = self.SAMPLE_DIR)
#        
#        f = open(ss_path, "w")
#        f.write(data_sample)
#        f.close()
#        options = [
#            "rna",
#            ss_path,
#            wdir,
#            self.DATA_DIR + self.GC_NAME,
#        ]
#        subprocess.check_call(['python', 'gcat_workflow'] + options)
#        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
#        self.assertTrue(success)
#
#    def test6_04_star_fusion_limited(self):
#        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
#        
#        data_sample = """[sra_fastq_dump],,,,
#sample,RUNID123456
#[star_fusion]
#sample
#""".format(sample_dir = self.SAMPLE_DIR)
#        
#        f = open(ss_path, "w")
#        f.write(data_sample)
#        f.close()
#        options = [
#            "rna",
#            ss_path,
#            wdir,
#            self.DATA_DIR + self.GC_NAME,
#        ]
#        subprocess.check_call(['python', 'gcat_workflow'] + options)
#        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
#        self.assertTrue(success)
#
#    def test6_11_star_fusion_limited(self):
#        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
#        
#        data_sample = """[fastq]
#sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
#[star_fusion]
#sample
#""".format(sample_dir = self.SAMPLE_DIR)
#        
#        f = open(ss_path, "w")
#        f.write(data_sample)
#        f.close()
#        options = [
#            "rna",
#            ss_path,
#            wdir,
#            self.DATA_DIR + self.GC2_NAME,
#        ]
#        subprocess.check_call(['python', 'gcat_workflow'] + options)
#        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
#        self.assertTrue(success)
#
#    def test6_12_star_fusion_limited(self):
#        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
#        
#        data_sample = """[bam_tofastq]
#sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
#[star_fusion]
#sample
#""".format(sample_dir = self.SAMPLE_DIR)
#        
#        f = open(ss_path, "w")
#        f.write(data_sample)
#        f.close()
#        options = [
#            "rna",
#            ss_path,
#            wdir,
#            self.DATA_DIR + self.GC2_NAME,
#        ]
#        subprocess.check_call(['python', 'gcat_workflow'] + options)
#        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
#        self.assertTrue(success)
#
#    def test6_13_star_fusion_limited(self):
#        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
#        
#        data_sample = """[bam_import]
#sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
#[star_fusion]
#sample
#""".format(sample_dir = self.SAMPLE_DIR)
#        
#        f = open(ss_path, "w")
#        f.write(data_sample)
#        f.close()
#        options = [
#            "rna",
#            ss_path,
#            wdir,
#            self.DATA_DIR + self.GC2_NAME,
#        ]
#        subprocess.check_call(['python', 'gcat_workflow'] + options)
#        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
#        self.assertTrue(success)
#
#    def test6_14_star_fusion_limited(self):
#        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
#        
#        data_sample = """[sra_fastq_dump],,,,
#sample,RUNID123456
#[star_fusion]
#sample
#""".format(sample_dir = self.SAMPLE_DIR)
#        
#        f = open(ss_path, "w")
#        f.write(data_sample)
#        f.close()
#        options = [
#            "rna",
#            ss_path,
#            wdir,
#            self.DATA_DIR + self.GC2_NAME,
#        ]
#        subprocess.check_call(['python', 'gcat_workflow'] + options)
#        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
#        self.assertTrue(success)
#
    def test7_11_iravnet_limited(self):
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test7_12_iravnet_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test7_13_iravnet_limited(self):
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test7_14_iravnet_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[sra_fastq_dump],,,,
sample,RUNID123456
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
            self.DATA_DIR + self.GC2_NAME,
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
        
        data_sample = """[bam_tofastq]
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

    def test8_04_ir_count_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[sra_fastq_dump],,,,
sample,RUNID123456
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

    def test8_11_ir_count_limited(self):
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test8_12_ir_count_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test8_13_ir_count_limited(self):
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test8_14_ir_count_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[sra_fastq_dump],,,,
sample,RUNID123456
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
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test9_01_kallisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[kallisto]
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

    def test9_02_kallisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
[kallisto]
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

    def test9_03_kallisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[kallisto]
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

    def test9_04_kallisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[sra_fastq_dump],,,,
sample,RUNID123456
[kallisto]
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

    def test9_11_kallisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[kallisto]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test9_12_kallisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
[kallisto]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test9_13_kallisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[kallisto]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test9_14_kallisto_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[sra_fastq_dump],,,,
sample,RUNID123456
[kallisto]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)


    def test10_11_juncmut_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[fastq]
sample,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[juncmut]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test10_12_juncmut_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_tofastq]
sample,{sample_dir}/H.Aligned.sortedByCoord.out.bam
[juncmut]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test10_13_juncmut_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sample,{sample_dir}/I.Aligned.sortedByCoord.out.bam
[juncmut]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test10_14_juncmut_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[sra_fastq_dump],,,,
sample,RUNID123456
[juncmut]
sample
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

    def test10_15_juncmut_limited(self):
        (wdir, ss_path) = func_path (self.DATA_DIR, sys._getframe().f_code.co_name)
        
        data_sample = """[bam_import]
sampleI,{sample_dir}/I.Aligned.sortedByCoord.out.bam
sampleL,{sample_dir}/L.Aligned.sortedByCoord.out.bam
sampleM,{sample_dir}/M.Aligned.sortedByCoord.out.cram
sampleN,{sample_dir}/N.Aligned.sortedByCoord.out.cram
[juncmut]
sampleI
sampleL
sampleM
sampleN
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data_sample)
        f.close()
        options = [
            "rna",
            ss_path,
            wdir,
            self.DATA_DIR + self.GC2_NAME,
        ]
        subprocess.check_call(['python', 'gcat_workflow'] + options)
        success = snakemake.snakemake(wdir + '/snakefile', workdir = wdir, dryrun = True)
        self.assertTrue(success)

if __name__ == '__main__':
    unittest.main()
