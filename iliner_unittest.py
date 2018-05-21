B#!/bin/usr/python

import numpy as np
import sys
import unittest

from iliner_load_parse import LoadFiles, FileParser
from iliner_make_files import FileConvert
from iliner_parse_vcf import VcfReader, VcfParser, VcfMakeFiles 

class TestFileParser(unittest.TestCase):
  
    def setUp(self):
        """Make instances with dummy good files and dummy bad files
        """
        good_sample = LoadFiles.load_sample("good_chr22_50p.sample")
        good_haps = LoadFiles.load_haps("good_chr22_50p.haps")
        bad_sample = LoadFiles.load_sample("bad_chr22_50p.sample")
        bad_haps = LoadFiles.load_haps("bad_chr22_50p.haps")

        self.good_set = LoadFiles(good_sample, good_haps)
        self.bad_set = LoadFiles(bad_sample, bad_haps)
    
    def test_good_check_sample_file(self):
        noerror = True
        try:
            FileParser.check_sample_file(self.good_set) 
        except ValueError:
            noerror = False
        self.assertEqual(noerror, True)

    def test_bad_check_sample_file(self):
        didraiseerror = False
        try:
            FileParser.check_sample_file(self.bad_set)
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True) 

    def test_good_same_sample(self):
        noerror = True
        try:
            FileParser.check_same_sample(self.good_set)
        except ValueError:
            noerror = False
        self.assertEqual(noerror, True)

    def test_bad_same_sample(self):
        didraiseerror = False
        try:
            FileParser.check_same_sample(self.bad_set)
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True)

    def test_good_check_alleles(self):
        noerror = True
        try:
            FileParser.check_alleles(self.good_set.haps_file.iloc[:,[3]].values, self.good_set.haps_file.iloc[:,[4]].values)
        except ValueError:
            noerror = False
        self.assertEqual(noerror, True)

    def test_bad_check_alleles(self):
        didraiseerror = False
        try:
            FileParser.check_alleles(self.bad_set.haps_file.iloc[:,[3]].values, self.bad_set.haps_file.iloc[:,[4]])
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True)

    def test_good_validate_nucleotides(self):
        nucfunc = np.vectorize(FileParser.nucleotide_test)
        good_refs = nucfunc(self.good_set.haps_file[[3]])
        good_alts = nucfunc(self.good_set.haps_file[[4]])
        noerror = True
        try:
            FileParser.validate_nucleotides(good_refs, "Ref All")
        except ValueError:
            noerror = False
        self.assertEqual(noerror, True)
        noerrors = True
        try:
            FileParser.validate_nucleotides(good_alts, 'Alt All')
        except ValueError:
            noerrors = False
        self.assertEqual(noerrors, True)


    def test_bad_validate_nucleotides(self):
        nucfunc = np.vectorize(FileParser.nucleotide_test)
        bad_refs = nucfunc(self.bad_set.haps_file[[3]])
        bad_alts = nucfunc(self.bad_set.haps_file[[4]])
        didraiseerror = False
        try:
            FileParser.validate_nucleotides(bad_refs, "Ref All")
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True)
        didraiseerrors = False
        try:
            FileParser.validate_nucleotides(bad_alts, 'Alt, All')
        except ValueError:
            didraiseerrors = True
        self.assertEqual(didraiseerrors, True)

    def test_good_validate_binary(self):
        noerror = True
        try:
            FileParser.validate_binary(self.good_set)
        except ValueError:
            noerror = False
        self.assertEqual(noerror, True)

    def test_bad_validate_binary(self):
        didraiseerror = False
        try:
            FileParser.validate_binary(self.bad_set)
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True)

class TestFileMaker(unittest.TestCase):

    def setUp(self):
        """Again make instances with dummy good files and dummy bad files
        """
        good_sample = LoadFiles.load_sample("good_chr22_50p.sample")
        good_haps = LoadFiles.load_haps("good_chr22_50p.haps")
        bad_sample = LoadFiles.load_sample("bad_chr22_50p.sample")
        bad_haps = LoadFiles.load_haps("bad_chr22_50p.haps")

        self.good_set = LoadFiles(good_sample, good_haps)
        self.bad_set = LoadFiles(bad_sample, bad_haps)

    def test_start_pedfile(self):
        FileParser.check_sample_file(self.good_set)
        ped_p = FileConvert.start_pedfile(self.good_set)
        self.assertEqual(len(ped_p[0]), 6)

    def test_bad_start_pedfile(self):
        try:
            FileParser.check_sample_file(self.bad_set)
        except ValueError:
            ped_p = FileConvert.start_pedfile(self.bad_set)
        self.assertNotEqual(len(ped_p[0]), 6)
    
    def test_convert_haps(self):
        hapslist = FileConvert.convert_haps(self.good_set)
        ref_alt = self.good_set.haps_file.loc[0,3:4].values.tolist()
        hapsalls = hapslist[0][5:]
        hapsall_l = []
        for a in hapsalls:
            if a not in hapsall_l:
                hapsall_l.append(a)
        self.assertEqual(sorted(ref_alt), sorted(hapsall_l))

class TestVCFCodes(unittest.TestCase):
    
    def setUp(self):
        
        good_vcf, good_comments = VcfReader.read_vcf("25p_chr22.phased.vcf")
        bad_vcf, bad_comments = VcfReader.read_vcf("bad_25p_chr22.phased.vcf")

        self.good_vcf = VcfReader(good_vcf, good_comments)
        self.bad_vcf = VcfReader(bad_vcf, bad_comments)
    
    def test_no_comments(self):
        didraiseerror = False
        try:
            VcfReader.read_vcf("no_comments_25p_chr22.phased.vcf")
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True)
                                  
    def test_good_validate_header(self):

        good_header = self.good_vcf.vcf.columns.values.tolist()
        noerror = True
        try:
            VcfParser.validate_header(good_header)
        except ValueError:
            noerror = False
        self.assertEqual(noerror, True)

    def test_bad_validate_header(self):
        bad_header = self.bad_vcf.vcf.columns.values.tolist()
        didraiseerror = False
        try:
            VcfParser.validate_header(bad_header)
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True)

    def test_good_validate_chrom(self):
        noerror = True
        try:
            VcfParser.validate_chrom(self.good_vcf)
        except ValueError:
            noerror = False
        self.assertEqual(noerror, True)

    def test_bad_validate_chrom(self):
        didraiseerror = False
        try:
            VcfParser.validate_chrom(self.bad_vcf)
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True)
    
    def test_good_diff_alt_ref(self):
        noerror = True
        try: 
            VcfParser.diff_alt_ref(self.good_vcf.vcf['REF'], self.good_vcf.vcf['ALT'])
        except ValueError:
            noerror = False
        self.assertEqual(noerror, True)

    def test_bad_diff_alt_ref(self):
        didraiseerror = False
        try:
            """
            Header had error used incorrect header names for this test
            """
            VcfParser.diff_alt_ref(self.bad_vcf.vcf['RF'], self.bad_vcf.vcf['ALM'])
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True)

    def test_good_validate_nucleotides(self):
        noerror = True
        nucfunc = np.vectorize(VcfParser.nucleotide_test)
        good_ref = nucfunc(self.good_vcf.vcf['REF'])
        good_alt = nucfunc(self.good_vcf.vcf['ALT'])
        try:
            VcfParser.validate_nucleotides(self.good_vcf, good_ref, 'Reference Allele')
        except ValueError:
            noerror = False
        self.assertEqual(noerror, True)
        noerrors = True
        try:
            VcfParser.validate_nucleotides(self.good_vcf, good_alt, 'Alternative Allele')
        except ValueError:
            noerrors = False
        self.assertEqual(noerrors, True)

    def test_bad_validate_nucleotides(self):
        didraiseerror = False
        nucfunc = np.vectorize(VcfParser.nucleotide_test)
        bad_ref = nucfunc(self.bad_vcf.vcf['RF'])
        bad_alt = nucfunc(self.bad_vcf.vcf['ALM'])
        try:
            VcfParser.validate_nucleotides(self.bad_vcf, bad_ref, 'Reference Allele')
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True)
        didraiseerrors = False
        try:
            VcfParser.validate_nucleotides(self.bad_vcf, bad_alt, 'Alternative Allele')
        except ValueError:
            didraiseerrors = True
        self.assertEqual(didraiseerrors, True)

    def test_good_validate_binary(self):
        noerror = True
        good_binary = np.isin(self.good_vcf.vcf.iloc[:,9:], ['0|0', '0|1', '1|0', '1|1'])
        try:
            VcfParser.validate_binary(self.good_vcf, good_binary)
        except ValueError:
            noerror = False
        self.assertEqual(noerror, True)

    def test_bad_validate_binary(self):
        didraiseerror = False
        bad_binary = np.isin(self.bad_vcf.vcf.iloc[:,9:], ['0|0', '0|1', '1|0', '1|1'])
        try:
            VcfParser.validate_binary(self.bad_vcf, bad_binary)
        except ValueError:
            didraiseerror = True
        self.assertEqual(didraiseerror, True)

    
        

if __name__ == '__main__':
    unittest.main()
