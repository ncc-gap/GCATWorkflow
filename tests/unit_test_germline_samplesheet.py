# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:18:34 2016

@author: okada

"""

import os
import sys
import shutil
import unittest
import gcat_workflow.germline.sample_conf as sc

BAM_IMP = "bam_import"
BAM_2FQ = "bam_tofastq"
ALN = "bwa_alignment_parabricks"
HT_CALL = "haplotypecaller_parabricks"
SUMMARY1 = "collect_wgs_metrics"
SUMMARY2 = "collect_multiple_metrics"

class ConfigureTest(unittest.TestCase):
    
    SAMPLE_DIR = "/tmp/temp-test/gcat_test_germline_samplesheet"
    REMOVE = True
    
    # init class
    @classmethod
    def setUpClass(self):
        os.makedirs(self.SAMPLE_DIR, exist_ok = True)
        touch_files = [
            "/A1.fastq",
            "/A2.fastq",
            "/B1.fq",
            "/B2.fq",
            "/C1_1.fq",
            "/C1_2.fq",
            "/C2_1.fq",
            "/C2_2.fq",
            "/A.markdup.cram",
            "/A.markdup.cram.crai",
            "/B.markdup.cram",
            "/B.markdup.crai",
            "/A.metadata.txt",
            "/B.metadata.txt",
            "/C.metadata.txt",
        ]
        for p in touch_files:
            open(self.SAMPLE_DIR + p, "w").close()
        
        link_files = [
            ("/C1_1.fq", "/link_C1_1.fq"),
            ("/C1_2.fq", "/link_C1_2.fq"),
            ("/A.markdup.cram", "/link_A.markdup.cram"),
            ("/A.markdup.crai", "/link_A.markdup.cram.crai"),
            ("/B.markdup.cram", "/link_B.markdup.cram"),
            ("/B.markdup.crai", "/link_B.markdup.crai"),
            ("/A.metadata.txt", "/link_A.metadata.txt"),
        ]
        for s in link_files:
            os.symlink(self.SAMPLE_DIR + s[0], self.SAMPLE_DIR + s[1])
        
    # terminated class
    @classmethod
    def tearDownClass(self):
        if self.REMOVE:
            shutil.rmtree(self.SAMPLE_DIR)

    # init method
    def setUp(self):
        pass

    # terminated method
    def tearDown(self):
        pass
    
    def test1_01_allinone(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq
pool2,{sample_dir}/link_C1_1.fq;{sample_dir}/link_C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq

[{bam2fq}]
A_control,{sample_dir}/A.markdup.cram
A_control2,{sample_dir}/link_A.markdup.cram

[{bamimp}]
pool3,{sample_dir}/B.markdup.cram
pool4,{sample_dir}/link_B.markdup.cram

[{htcall}]
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
A_control2,{sample_dir}/link_A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ, bamimp = BAM_IMP, htcall = HT_CALL, summary1 = SUMMARY1, summary2 = SUMMARY2)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        sample_conf = sc.Sample_conf(ss_path)
        
        self.assertEqual(sample_conf.fastq, {
            'A_tumor': [[self.SAMPLE_DIR + '/A1.fastq'], [self.SAMPLE_DIR + '/A2.fastq']], 
            'pool1': [[self.SAMPLE_DIR + '/B1.fq'], [self.SAMPLE_DIR + '/B2.fq']], 
            'pool2': [[self.SAMPLE_DIR + '/link_C1_1.fq', self.SAMPLE_DIR + '/link_C1_2.fq'], [self.SAMPLE_DIR + '/C2_1.fq', self.SAMPLE_DIR + '/C2_2.fq']],
        })
        
        self.assertEqual(sample_conf.fastq_src, {
            'A_tumor': [[self.SAMPLE_DIR + '/A1.fastq'], [self.SAMPLE_DIR + '/A2.fastq']], 
            'pool1': [[self.SAMPLE_DIR + '/B1.fq'], [self.SAMPLE_DIR + '/B2.fq']], 
            'pool2': [[self.SAMPLE_DIR + '/link_C1_1.fq', self.SAMPLE_DIR + '/C1_1.fq', self.SAMPLE_DIR + '/link_C1_2.fq', self.SAMPLE_DIR + '/C1_2.fq'], [self.SAMPLE_DIR + '/C2_1.fq', self.SAMPLE_DIR + '/C2_2.fq']],
        })

        self.assertEqual(sample_conf.bam_tofastq, {'A_control': self.SAMPLE_DIR + '/A.markdup.cram', 'A_control2': self.SAMPLE_DIR + '/link_A.markdup.cram'})
        self.assertEqual(sample_conf.bam_tofastq_src, {
            'A_control': [self.SAMPLE_DIR + '/A.markdup.cram'],
            'A_control2': [self.SAMPLE_DIR + '/link_A.markdup.cram', self.SAMPLE_DIR + '/A.markdup.cram']
        })
        self.assertEqual(sample_conf.bam_import, {'pool3': self.SAMPLE_DIR + '/B.markdup.cram', 'pool4': self.SAMPLE_DIR + '/link_B.markdup.cram'})
        self.assertEqual(sample_conf.bam_import_src, {
            'pool3': [self.SAMPLE_DIR + '/B.markdup.cram', self.SAMPLE_DIR + '/B.markdup.crai'],
            'pool4': [self.SAMPLE_DIR + '/link_B.markdup.cram', self.SAMPLE_DIR + '/link_B.markdup.crai', self.SAMPLE_DIR + '/B.markdup.cram', self.SAMPLE_DIR + '/B.markdup.crai']
        })
        self.assertEqual(sample_conf.haplotype_call, ['A_tumor','A_control','pool1','pool2','pool3'])
        self.assertEqual(sample_conf.manta, ['A_tumor'])
        self.assertEqual(sample_conf.melt, ['A_tumor'])
        self.assertEqual(sample_conf.gridss, ['A_tumor'])
        self.assertEqual(sample_conf.wgs_metrics, ['A_tumor'])
        self.assertEqual(sample_conf.multiple_metrics, ['A_tumor'])
        self.assertEqual(sample_conf.readgroup, {
            'A_tumor': self.SAMPLE_DIR + '/A.metadata.txt',
            'pool1': self.SAMPLE_DIR + '/B.metadata.txt',
            'pool2': self.SAMPLE_DIR + '/C.metadata.txt',
            'A_control': self.SAMPLE_DIR + '/A.metadata.txt',
            'A_control2': self.SAMPLE_DIR + '/link_A.metadata.txt'
        })
        self.assertEqual(sample_conf.readgroup_src, {
            'A_tumor': [self.SAMPLE_DIR + '/A.metadata.txt'],
            'pool1': [self.SAMPLE_DIR + '/B.metadata.txt'],
            'pool2': [self.SAMPLE_DIR + '/C.metadata.txt'],
            'A_control': [self.SAMPLE_DIR + '/A.metadata.txt'],
            'A_control2': [self.SAMPLE_DIR + '/link_A.metadata.txt', self.SAMPLE_DIR + '/A.metadata.txt']
        })

    # --------------------------------------------------------------------
    # Not Exist
    # --------------------------------------------------------------------
    def test2_01_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A3.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,
pool2,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test2_02_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B3.fq,,
pool2,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_2.fq,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test2_03_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq,,
pool2,{sample_dir}/C1_1.fq;{sample_dir}/C1_2.fq,{sample_dir}/C2_1.fq;{sample_dir}/C2_3.fq,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test2_04_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bam2fq}],,,,
A_control,{sample_dir}/X.markdup.cram,,,
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test2_05_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bamimp}],,,,
pool3,{sample_dir}/X.markdup.cram,,,
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test2_06_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
[readgroup]
A_tumor,{sample_dir}/A1.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)
    
    def test2_07_not_exists(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bam2fq}]
A_control,{sample_dir}/A.markdup.cram
[readgroup]
A_control,{sample_dir}/A1.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)
        
    # --------------------------------------------------------------------
    # Not defined
    # --------------------------------------------------------------------
    def test3_01_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq

[{htcall}]
B_tumor
""".format(sample_dir = self.SAMPLE_DIR, htcall = HT_CALL)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test3_02_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test3_03_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq
pool1,{sample_dir}/B1.fq,{sample_dir}/B2.fq
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)
        
    def test3_04_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bam2fq}]
A_control,{sample_dir}/A.markdup.cram
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test3_05_undefine(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bam2fq}]
A_tumor,{sample_dir}/A.markdup.cram
A_control,{sample_dir}/A.markdup.cram
[readgroup]
A_tumor,{sample_dir}/A.metadata.txt
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)
        
    # --------------------------------------------------------------------
    # Duplicate
    # --------------------------------------------------------------------
    def test4_01_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_02_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bam2fq}]
A_control,{sample_dir}/A.markdup.cram
A_control,{sample_dir}/A.markdup.cram
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_03_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bamimp}],,,,
pool3,{sample_dir}/B.markdup.cram,,,
pool3,{sample_dir}/B.markdup.cram,,,
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_04_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
[{bam2fq}],,,,
A_tumor,{sample_dir}/A.markdup.cram,,,
[{bamimp}],,,,
pool1,{sample_dir}/B.markdup.cram,,,
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ, bamimp = BAM_IMP)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_05_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
[{bam2fq}],,,,
pool3,{sample_dir}/A.markdup.cram,,,
[{bamimp}],,,,
pool3,{sample_dir}/B.markdup.cram,,,
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ, bamimp = BAM_IMP)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_06_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
pool3,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,,
[{bam2fq}],,,,
A_tumor,{sample_dir}/A.markdup.cram,,,
[{bamimp}],,,,
pool3,{sample_dir}/B.markdup.cram,,,
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ, bamimp = BAM_IMP)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_08_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq

[{htcall}]
A_tumor
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, htcall = HT_CALL)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_09_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq

[manta]
A_tumor
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, htcall = HT_CALL)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test4_10_duplicate(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq]
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq

[melt]
A_tumor
A_tumor
""".format(sample_dir = self.SAMPLE_DIR, htcall = HT_CALL)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    # --------------------------------------------------------------------
    # Format error
    # --------------------------------------------------------------------
    def test5_01_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_02_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[fastq],,,,
A_tumor,{sample_dir}/A1.fastq,{sample_dir}/A2.fastq,{sample_dir}/A2.fastq
""".format(sample_dir = self.SAMPLE_DIR)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_03_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bam2fq}]
{sample_dir}/A.markdup.cram
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_04_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bam2fq}]
A_tumor,{sample_dir}/A.markdup.cram,{sample_dir}/B.markdup.cram
""".format(sample_dir = self.SAMPLE_DIR, bam2fq = BAM_2FQ)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_05_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bamimp}]
{sample_dir}/B.markdup.cram
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)

    def test5_06_unformat(self):
        ss_path = self.SAMPLE_DIR + "/" + sys._getframe().f_code.co_name + ".csv"
        data = """[{bamimp}]
A_tumor,{sample_dir}/A.markdup.cram,{sample_dir}/B.markdup.cram
""".format(sample_dir = self.SAMPLE_DIR, bamimp = BAM_IMP)
        
        f = open(ss_path, "w")
        f.write(data)
        f.close()
        try:
            fail = False
            sc.Sample_conf(ss_path)
        except Exception as e:
            print(e)
            fail = True

        self.assertTrue(fail)


if __name__ == '__main__':
    unittest.main()

