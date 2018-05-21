#!/bin/usr/python

from datetime import datetime
import numpy as np
import os
import sys
from iliner_load_parse import LoadFiles, FileParser


class FileConvert(LoadFiles):
    """FileConvert accpets the validated files.
    """         

    def start_pedfile(self):
        """Returns the first six columns of the pedfile which will be used for the first six columns of the pedigree file

        Args:
            self(:obj:)

        Return:
            :obj:`numpy array`: first six columns of the pedigree file
        """
        sample_file = self.sample.as_matrix()
        ped_f = np.insert(sample_file, 2, '0', axis=1)
        ped_m = np.insert(ped_f, 3, '0', axis=1)
        ped_p = np.insert(ped_m, -1, '-9', axis=1)
        return ped_p

    def convert_haps(self):
        """Converst binary (0,1) to (A,T,C, or G)                                                                                                                                                                                                                            
            Args:                                                                                                                                                            
                haps_file(:obj:`numpy array`): validated haplotype file                                                                                                                                                                             
            Returns:                                                                                                                                                         
                :obj:`list`: list of converted values                                                                                                                        
        """
        haps_list = self.haps_file.values.tolist()
        hapslist = []
        for line in haps_list:
            base1 = line[3]
            base2 = line[4]
            nums = line[5:]
            nums = map(lambda x: x.replace('0', base1), nums)
            nums = map(lambda x: x.replace('1', base2), nums)
            hapslist.append(list(nums))
        return hapslist

    @staticmethod
    def make_pedfile(hapslist, pedlen, pedfile_start):
        """Formats and concatenates the input for the pedigree file                                                                                                      
                                                                                                                                                                             
            Args:                                                                                                                                                            
                hapslist(:obj:`list`): list of converted haplotype values                                                                                                    
                pedlen(:obj:`int`): number of people                                                                                                                         
                pedfile_start(:obj:`numpy array`): first six columns of the pedigree files                                                                                   
                                                                                                                                                                             
            Returns:                                                                                                                                                         
            :obj:`numpy array`: the input for the pedigree file                                                                                                          
        """
        lenhaps = int(len(hapslist))
        t_hapslist = np.array(hapslist)
        r_hapslist = t_hapslist.reshape(lenhaps, pedlen, 2)
        c_hapslist = np.concatenate(r_hapslist[0:], axis=1)
        final_list = np.append(pedfile_start, c_hapslist, axis=1)
        return final_list


    def make_gpos_mapfile(self, hapmap, outfile):
        """Produces a map file with genetic position. It does so by interpolating genetic position. Citation: Joe Pickrell interpolate_maps.py
        
        Args:
            haps_file(:obj:`numpy array`): validated haplotype file
            hapmap(:obj:`file`): file used to compute genetic position
            outfile(:obj:`filename`): name of output file

        Returns:
            :obj:`file`: map file with genetic positon
            output format: [chrom] [rsid] [genetic position] [physical position]
        """
        mappos = list()
        mapgpos = list()

        chrom = 0 
        chromin = self.haps_file[chrom].values.tolist()
        rsid = 1
        rsin = self.haps_file[rsid].values.tolist()
        position = 2
        posin = self.haps_file[position].values.tolist()
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
                ## Tests if physical position from haplotype file is equal to physical postion in hapmap file. Returns corresponding genetic position
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
        


                                  
                            
                
