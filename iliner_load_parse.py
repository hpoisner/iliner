#!/usr/bin/python3

from datetime import datetime
import logging
import numpy as np
import pandas as pd
import os
import sys

class LoadFiles(object):
    """LoadFiles loads sample files and haplotype files
       Files must be phased
    """

    def __init__(self, sample, haps_file):
        self.sample = sample
        self.haps_file = haps_file
    
    @staticmethod
    def load_sample(path):
        """
        Args:
            path(:obj:`filepath`): file path for sample file

        Raises:
            :obj:`IOError`: if file provided does not exist

        Returns:
            :obj:`pandas dataframe`: file as a dataframe
            Format:
            [ID_1] [ID_2] [missing] [covariates]
        """
        if os.path.isfile(path):
            samp_file = pd.read_table(path, dtype=str, sep=' ', skiprows=[1])
            return samp_file
        else:
            raise IOError("Sample file {} does not exist.".format(path))

    @staticmethod
    def load_haps(path):
        """

        Args:
            path(:obj:`filepath`): file path for haplotypefile
            
        Raises:
            :obj:`IOError`: if file provided does not exist

        Returns:
            :obj:`pandas dataframe`: file as a dataframe
            Format:
            [Chromosome #] [RSID] [Physical Position] [Reference Allele] [Alternative Allele] [Binary Representations of Alleles]
        """
        if os.path.isfile(path):
            haps_list = []
            with open(path) as f:
                for line in f:
                    line = line.strip('\n').split(' ')
                    haps_list.append(line)            
            haps_file = pd.DataFrame(haps_list)
            return haps_file
        else:
            raise IOError("Haplotype file {} does not exist.".format(path))



class FileParser(LoadFiles):
    """FileParser parses the loaded sample and haplotype files
    """

    
    def check_sample_file(self):
        """Ensures that the sample file includes the proper information.
        
        Args:
            self.sample(:obj:`pandas datatrame`): the sample file as a dataframe

        Raises:
            :obj:`ValueError`: if file is missing [ID_1, ID_2, or sex]

        Returns:
            :obj:`pandas dataframe`: sample file without headers and validated 
        """
        sample_header = self.sample.columns.values.tolist()
        key_columns = ['ID_1', 'ID_2','sex']
        desired_columns = np.isin(sample_header, key_columns)
        good_locs = np.where(desired_columns)
        actual_locs = good_locs[0].tolist()
        if len(actual_locs) != 3:
            raise ValueError("Your sample file should contain columns called ID_1, ID_2, and sex.")
        else:
            self.sample = self.sample[['ID_1', 'ID_2', 'sex']]
            

    def remove_indels(self):
       """
       Args:
           self

       Returns:
           self.haps_file without indels
       """
       mask = (self.haps_file[[3]].str.len() == 1) & (self.haps_file[[4]].str.len() == 1)
       self.haps_file = self.haps_file.loc[mask]
                         
    def check_same_sample(self):
        """Ensures that for each person in the sample file, they have two haplotypes.
        
        Args:
            sample_file(:obj:`pandas dataframe`): the sample file as a dataframe
            haps_file(:obj:`pandas dataframe`): the haplotype file as a dataframe

        Raises:
            :obj:`ValueError`: Each person does not have exactly 2 haplotypes
        """
        sample_shape = self.sample.shape
        haps_shape = self.haps_file.shape
        if sample_shape[0] == (haps_shape[1] - 5)/2:
            pass
        else:
            raise ValueError("Each person does not have two haplotypes.")

    def check_alleles(ref, alt):
        """Ensures that the reference and alternative allles are different

        Args:
            haps_file(:obj:`pandas dataframe`): the haplotype filea as a dataframe

        Raises:
            :obj:`ValueError`: References and alternative alleles should be different.
            Also provides location list
        """
        diff_alleles = np.array(ref == alt)
        true_values = np.where(diff_alleles == True)[0]
        true_locs = true_values.tolist()
        true_locs = [int(i)+1 for i in true_locs]
        if len(true_locs) > 0:
            raise ValueError("In these rows your reference and alternative alleles are recorded as the same allele: {}.".format(' '.join(map(str, true_locs))))
        else:
            pass

    @staticmethod
    def nucleotide_test(nucleotide):
        """Tests if nucleotides are A,C,T, or G
        
        Args:
            nucleotide(:obj:`string`): The string of the nucleotides

        Returns:
            :obj:`boolean value`: True or False
        """
        allowed_chars = set('ACTG')
        if set(nucleotide.upper()).issubset(allowed_chars):
            return True
        else:
            return False

    @staticmethod
    def validate_nucleotides(array, allele_type):
        """Validates nucleotides at all locations

        Args:
            array(:obj:`boolean array`): Array of True and False
            allele_type(:obj:`string`): String describng alleles in array as reference or alternative

        Raises:
            :obj:`ValueError`: A message explaining where the false values are in the file
        """
        false_values = np.where(array == False)[0]
        false_locs = false_values.tolist()
        false_locs = [int(i)+1 for i in false_locs] 
        if len(false_locs) > 0:
            raise ValueError("In rows {} your {} are not A, T, C, or G.".format(' '.join(map(str, false_locs)),allele_type))
        else:
            pass

    def validate_binary(self):
        """Validates the binary columns.
        
        Args:
            haps_file(:obj:`numpy array`): Numpy array of the haplotype file
        
        Raises:
            :obj:`ValueError`: A message explaining where the values are not binary
        """
        log_array_haps = np.logical_or(self.haps_file.iloc[:,5:] == '0', self.haps_file.iloc[:,5:] == '1')
        l = (list(zip(*np.where(log_array_haps == False))))
        k  = (1,1)
        true_values = [(i[0]+k[0],i[1]+ k[1]) for i in l]
        if len(true_values) > 0:
            for j in true_values:
                print('In Row: {} and Column: {} your value is not 1 or 0'.format(j[0], j[1]))
            raise ValueError("Please check incorrect values")
        else:
            pass



"""
samp_file = LoadFiles.load_sample(sys.argv[1])
print(samp_file.head())

haps_file = LoadFiles.load_haps(sys.argv[2])
print(haps_file.head())

ilash_obj = LoadFiles(samp_file, haps_file)
print(ilash_obj.__dict__)

FileParser.check_sample_file(ilash_obj)
print(ilash_obj.sample.head())

FileParser.check_same_sample(ilash_obj)

FileParser.check_alleles(ilash_obj.haps_file.iloc[:,[3]].values, ilash_obj.haps_file.iloc[:,[4]].values)

nucfunc = np.vectorize(FileParser.nucleotide_test)

log_array_ref = nucfunc(ilash_obj.haps_file[[3]])
log_array_alt = nucfunc(ilash_obj.haps_file[[4]])
FileParser.validate_nucleotides(log_array_ref, 'Reference Allele')
FileParser.validate_nucleotides(log_array_alt, 'Alternative Allele')
FileParser.validate_binary(ilash_obj)
"""
