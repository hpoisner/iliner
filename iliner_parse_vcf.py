#!/bin/usr/python

import gzip
import numpy as np
import pandas as pd
import os
import sys

class VcfReader(object):
    """Loads vcf file 
       Files must be phased
    """

    def __init__(self, vcf, comments):
        self.vcf = vcf
        self.comments = comments
    
    @staticmethod
    def read_vcf(filename):
        """
        Args: 
            filename(:obj:`filepath`): file path for the VCF

        Raises: 
            :obj:`IOError`: if fiel proovided does not exist

        Returns:
            :obj:`pandas dataframe`: file as a pandas dataframe
            :obj:`int`: number of comments in the VCF file
        """
        if os.path.isfile(filename):
            comp = 'gzip'
            comments = VcfReader.count_comments(filename)
            if filename.endswith('.gz'):
                vcf_file = pd.read_table(filename, compression=comp, skiprows=comments, dtype=str)
                return vcf_file, comments
            else:
                vcf_file = pd.read_table(filename, skiprows=comments, dtype=str)
                return vcf_file, comments
        else:
            raise IOError('VCF file {} does not exist.'.format(filename))

        
    @staticmethod
    def count_comments(filename):
        """
        Args:
            filename(:obj:`filepath`): file path for VCF
        
        Raises:
            :obj:`ValueError`: if file provided does not have comments 

        Returns:
            :obj:`int`: number of comments in VCF
        
        """
        comments = 0
        if filename.endswith('.gz'):
            with gzip.open(filename) as fh:
                for line in fh:
                    if line.startswith(b'##'):
                        comments += 1
                    else:
                        break
            if comments == 0:
                raise ValueError('VCF files must have comments')
        else:
            with open(filename) as fh:
                for line in fh:
                    if line.startswith('##'):
                        comments += 1
                    else:
                        break
            if comments == 0:
                raise ValueError('VCF files must have commnets')
        return comments




class VcfParser(VcfReader):
    """Validates the VCF values
    """
    @staticmethod
    def validate_header(vcf_header):
        """
        Args:
            vcf_header:`list`: header of the VCF file

        Raise:
            :obj:`ValueError`: if the requird columns are missing from the VCF

        """
        required_columns = ['#CHROM','POS','ID','REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        intersected_list = list(set(vcf_header) & set(required_columns))
        if sorted(intersected_list) != sorted(required_columns):
            missing_values = [i for i in required_columns if i not in intersected_list]
            raise ValueError("You are missing {} which are required VCF columns.".format(' '.join(map(str, missing_values))))
        else:
            pass

    def validate_chrom(self):
        """
        Args:
            self

        Raises:
            :obj:`ValueError`: if there is more than one chromosome in the VCF file
        """
        chrom_list = list(set(self.vcf['#CHROM']))
        if len(chrom_list) != 1:
            raise ValueError('Please only pass a VCF with one chromosome')
        else:
            pass

    @staticmethod
    def diff_alt_ref(ref, alt):
        """
        Args:
            ref: Column in the pandas VCF with the reference allele/s
            alt: Column in the pandas VCF with the alternative allele/s
      
        Raises:
            :obj:`ValueError`: if the reference and alternative alleles are the same 
        """
        diff_alleles = np.array(ref == alt)
        true_values = np.where(diff_alleles == True)[0]
        true_list = true_values.tolist()
        true_locs = [int(i)+1 for i in true_list]
        if len(true_locs) > 0:
            raise ValueError("In these rows your reference and alternative alleles are recorded as the same allele: {}.".format(' '.join(map(str, true_locs))))
        else:
            pass

    def remove_indels(self):
        """
        Args:
            self
        
        Returns:
            self.vcf without indels
        """
        mask = (self.vcf['REF'].str.len() == 1) & (self.vcf['ALT'].str.len() == 1)
        self.vcf = self.vcf.loc[mask]

    @staticmethod
    def nucleotide_test(nucleotide):
        """
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
    
    def validate_nucleotides(self, array, allele_type):
        """
        Args:
            self
            array:`boolean array`: an array that shows where the values were not ATCG
            allele_type:`string`: a string of reference or alternative allele

        Raises:
            :obj:`ValueError`: Alleles are not acceptable nucleotides
        """
        comments = int(self.comments) + 2
        false_values = np.where(array == False)[0]
        false_list = false_values.tolist()
        false_locs = [int(i)+comments for i in false_list]
        if len(false_locs) > 0:
            raise ValueError("In rows {} your {} are not A, T, C, or G.".format(' '.join(map(str, false_locs)), allele_type))
        else:
            pass

    def validate_binary(self, binary_array):
        """
        Args:
            self
            binary_array:`boolean array`: an array that shows were values are true or false

        Raises:
            :obj:`ValueError`: Values are not corret.
        """
        comments = int(self.comments) + 2
        l = (list(zip(*np.where(binary_array == False))))
        k = (comments, 1)
        true_values = [(i[0] +k[0], i[1] + k[1]) for i in l]
        if len(true_values) > 0:
            for j in true_values:
                print('In Row: {} and Column: {} your value is not 0|0, 0|1, 1|0, or 1|1'.format(j[0], j[1]))
            raise ValueError('Please check incorrect values')
        else:
            pass
            
class VcfMakeFiles(VcfReader):
    """Makes the pedigree file and the map file
    """
    
    @staticmethod
    def get_people(vcf_header):
        """
        Args:
            vcf_header(:obj:`list`): VCF header as a list

        Returns:
            :obj:`list`: list of the individual ids from the VCF file
        """
        people = vcf_header[9:]
        return people

    @staticmethod
    def start_pedfile(people):
        """
        Args:
            people(:obj:`list`) list of indidivual ids from the VCF file

        Returns:
           :obj:`array`: start of the pedigree file as an array
        """
        people_array = np.asarray(people)
        people_t = people_array.reshape((-1,1))
        ids = np.concatenate((people_t, people_t), axis=1)
        ids_f = np.insert(ids, 2, '0', axis=1)
        ids_m = np.insert(ids_f, 3, '0', axis=1)
        ids_p = np.insert(ids_m, 4, '-9', axis=1)
        ids_s = np.insert(ids_p, 5, '1', axis=1)
        return ids_s

    def convert_haps(self):
        """
        Args:
           self
        
        Returns:
           :obj:`list`: a list of haplotypes
        """
        hapslist = []
        vcf_haps = self.vcf.values.tolist()
        for line in vcf_haps:
            split_nums = []
            ref_all = line[3]
            alt_all = line[4]
            nums = line[9:]
            for i in nums:
                i.split('|')
                split_nums.append(i[0])
                split_nums.append(i[2])
            split_nums = map(lambda x: x.replace('0', ref_all), split_nums)
            split_nums = map(lambda x: x.replace('1', alt_all), split_nums)
            hapslist.append(list(split_nums))
        return hapslist
            
    def make_pedfile(hapslist, pedlen, pedfile_start):
        """
        Args:
            hapslist: list of the haplotype values as nucleotides
            pedlen: length of the pedigree file
            pedfile_start: the first five columns of the pedigree file
         
        Returns:
            :obj:`list`: list that will be written out to be the pedigree file
        """
        lenhaps = int(len(hapslist))
        t_hapslist = np.array(hapslist)
        r_hapslist = t_hapslist.reshape(lenhaps, pedlen, 2)
        c_hapslist = np.concatenate(r_hapslist[0:], axis=1)
        final_list = np.append(pedfile_start, c_hapslist, axis=1)
        return final_list

    def make_gpos_mapfile(self, hapmap, outfile):
        """Produces a map file with genetic position. It does so my interpolating genetic position.
           Code is a modified version of the cited code below
           Citation:
           Author: Joe Pickrell
           Date: 5/29/2018
           Source Code: interpolate_maps.py
           Code Version: Version Commited Jun 19, 2014
           Availability: https://github.com/joepickrell/1000-genomes-genetic-map/blob/master/scripts/interpolate_maps.py

        Args:
            self
            hapmap:file pathway to genetic map file
            outfile: pathway to output file

        Returns:
             Mapfile with genetic position

        """
        mappos = list()
        mapgpos = list()
        
        chromin = self.vcf['#CHROM'].values.tolist()
        rsin = self.vcf['ID'].values.tolist()
        posin = self.vcf['POS'].values.tolist()
        
        line = hapmap.readline()
        line = hapmap.readline()
        while line:
            line = line.strip().split()
            pos = int(line[1])
            gpos = float(line[3])
            mappos.append(pos)
            mapgpos.append(gpos)
            line = hapmap.readline()

        index1 = 0
        index2 = 0
        while index1 < len(posin):
            pos = int(posin[index1])
            rs = rsin[index1]
            chrom = chromin[index1]
            if pos == mappos[index2]:
                ##the 1000 Genomes site was genotopyes as part of the map (comment directly from original code)                              
                ## Tests if physical position from haplotype file is equal to physical postion in hapmap file. Returns corresponding genetic  position                                                                                                                                    
                outfile.write(' '.join([chrom, rs, str(mapgpos[index2]), str(pos)]) +'\n')
                index1 += 1
            elif pos < mappos[index2]:
                ## current position in interpolation before marker                                                                           
                if index2 == 0:
                    ## before the first site in the map (genetic position = 0)                                                               
                    outfile.write(' '.join([chrom, rs, str(mapgpos[index2]), str(pos)]) + '\n')
                    index1 += 1
                else:
                    ## interpolate                                                                                                           
                    ## This part calculates the approximate genetic position if the current physical position from the haplotype file is less than than the                                                                                                                              
                    prevg = mapgpos[index2 - 1]
                    prevpos = mappos[index2]
                    frac = (float(pos) - float(mappos[index2 - 1]))/ (float(prevpos) - float(mappos[index2 - 1]))
                    tmpg = prevg + frac* (mapgpos[index2] - prevg)
                    outfile.write(' '.join([chrom, rs, str(tmpg), str(pos)]) + '\n')
                    index1 +=1
            elif pos > mappos[index2]:
                ## current position in interpolation after marker                                                                            
                if index2 == len(mappos) - 1:
                    ## after the last site in the map (genetic position = maximum in map)                                                    
                    outfile.write(' '.join([chrom, rs, str(mapgpos[index2]), str(pos)]) + '\n')
                    index1 += 1
                else:
                    ##increment the marker                                                                                                   
                    index2 += 1
        outfile.close()


"""
vcf_file, comments = VcfReader.read_vcf(sys.argv[1])
vcf_obj = VcfReader(vcf_file, comments)
vcf_header = vcf_file.columns.values.tolist()
VcfParser.validate_header(vcf_header)
VcfParser.validate_chrom(vcf_obj)
VcfParser.remove_indels(vcf_obj)
VcfParser.diff_alt_ref(vcf_obj.vcf['REF'], vcf_obj.vcf['ALT'])
## update with binary array method
nucfunc = np.vectorize(VcfParser.nucleotide_test)
log_array_ref = nucfunc(vcf_obj.vcf['REF'])
log_array_alt = nucfunc(vcf_obj.vcf['ALT'])
VcfParser.validate_nucleotides(vcf_obj, log_array_ref, 'Reference Allele')
VcfParser.validate_nucleotides(vcf_obj, log_array_alt, 'Alternative Allele')
binary_array = np.isin(vcf_obj.vcf.iloc[:,9:], ['0|0', '0|1', '1|0', '1|1'])
VcfParser.validate_binary(vcf_obj, binary_array)
people = VcfMakeFiles.get_people(vcf_header)
pedfile_start= VcfMakeFiles.start_pedfile(people)
hapslist = VcfMakeFiles.convert_haps(vcf_obj)
pedlen = int(len(hapslist[0])/2)
final_list = VcfMakeFiles.make_pedfile(hapslist, pedlen, pedfile_start)
output = open(sys.argv[3], 'w')
VcfMakeFiles.make_gpos_mapfile(vcf_obj, open(sys.argv[2]), output)
with open(sys.argv[4], "w") as f:
    for line in final_list:
        f.write(' '.join(line) + '\n')
f.close()
"""
