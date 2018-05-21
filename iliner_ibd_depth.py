
#!/bin/usr/python3

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'figure.max_open_warning': 30})
import numpy as np
import pandas as pd 
import os
import statistics
import sys


class IBDDepth:

    @staticmethod
    def load_files(path):
        """

        Args:
            path(:obj:`filepath`): file path for map file 

        Raises:
            :obj:`IOError`: if file provided does not exist

        Returns:
            :obj:`pandas dataframe`: file as a dataframe
        """
        if os.path.isfile(path):
            input_file = pd.read_table(path, sep = ' ', names=['Chromosome','RSID', 'Genetic_Position', 'Physical_Position'], usecols=['Chromosome','RSID', 'Physical_Position'])
            input_file['index'] = list(range(0,len(input_file)))
            input_file.set_index('RSID',inplace=True)
            return input_file
        else:
            raise IOError("{} does not exist.".format(path))

    @staticmethod
    def load_ilash(i_file):
        """
        
        Args:
            i_file(:obj:`filepath`): file path for iLASH file

        Raises:
            :obj:`IOError`: if file provided does not exist

        Returns:
            :obj:`pandas dataframe`: file as a dataframe
        """
        if os.path.isfile(i_file):
            ilash_file = pd.read_table(i_file, sep='\t', header=None, names=['Person_1', 'Hap_1', 'Person_1', 'HAP_1', 'CHR', 'START', 'STOP', 'SNP1', 'SNP2', 'IBD_Sharing', 'Measure'])
            snps = ilash_file[['SNP1', 'SNP2']]
            return ilash_file, snps
        else:
            raise IOError('{} does not exist'.format(i_file))
                        

    @staticmethod
    def find_ibd(mapfile, snps):
        """
        
        Args:
           mapfile(:obj:`dataframe`): map file as a pandas dataframe
           snps(:obj:`dataframe`): snps from the iLASH file that represent the start and stop values of IBD sharing

        Returns:
            :obj:`pandas dataframe`: mapfile that shows how much ibd sharing exists at the location
        """
        ibdlist = np.zeros(len(mapfile))
        for chunk_num, chunk in snps.groupby(np.arange(len(snps))//(len(snps)//2000)):
            if chunk_num%100==0:
                print('Processing chunk #: ', chunk_num)
            start_idxs = list(mapfile.loc[chunk['SNP1'], 'index'])
            stop_idxs = list(mapfile.loc[chunk['SNP2'], 'index'])
            for i in range(0, len(start_idxs)):
                start_idx = start_idxs[i]
                stop_idx = stop_idxs[i]
                ibdlist[start_idx:(stop_idx+1)]+= 1
        mapfile['IBD_Sharing'] = ibdlist
        return mapfile

    @staticmethod
    def get_stdev(mapfile):
        """

        Args:
            mapfile(:obj:`dataframe`): mapfile with ibd sharing

        Returns:
            four_up:obj:`float`: four standard deviations above provided genome wide mean
            four_down:obj:`float`: four standard deviations below provided genome wide mean
        """
        all_depth_array = np.asarray(mapfile['IBD_Sharing'])
        genome_wide_mean = np.mean(all_depth_array)
        four_up = genome_wide_mean + (4 * np.std(all_depth_array))
        four_down = genome_wide_mean - (4 * np.std(all_depth_array))
        return four_up, four_down
                                
        

    @staticmethod
    def make_depth_plot(mapfile, prefix, four_up, four_down):
        """
        
        Args:
            mapfile(:obj:`dataframe`): map file with ibd
            prefix(:obj;`filepath + prefix string`): path to output with provided prefix
            four_up(:obj:`float`): four standard deviations above provided genome wide mean
            four_down(:obj:`float`): four standard deviations below provided genome wide mean

        Returns:
            figures of ibd distributions
        """
        grouped_map = mapfile.groupby(['Chromosome'])
        for key, item in grouped_map:
            pos_share = grouped_map[['Physical_Position', 'IBD_Sharing']].get_group(key)
            pos = pos_share['Physical_Position'].reset_index(drop=True).tolist()
            shared = pos_share['IBD_Sharing'].reset_index(drop=True).tolist()
            posi = [xi/1000000 for xi in pos]
            fig, ax = plt.subplots()
            ax.scatter(posi, shared, c='xkcd:black')
            ax.set_title('Chromosome ' + str(key))
            ax.set_xlabel('Physical Position (Mb)')
            ax.set_ylabel('Depth of Pairwise IBD Sharing')
            plt.axhline(y=four_up, color='r', linestyle='--')
            plt.axhline(y=four_down, color='r', linestyle='--')
            plt.savefig(prefix + '_chr' + str(key) + '_ibd_depth.png')
            plt.cla()

    @staticmethod
    def remove_outlier(mapfile, ilash_file, prefix, four_up, four_down):
        """

        Args:
            mapfile(:obj:`pandas dataframe`): mapfile with ibd sharing
            ilash_file(:obj:`pandas dataframe`): iLASH file as df
            prefix(:obj:`filepath + string`): filepath with prefix string
            four_up(:obj:`float`): four standard deviations above provided genome wide mean
            four_down(:obj:`float`): four standard deviations below provided genome wide mean

        Returns:
           QCed iLASH file
        """
        big =  mapfile['IBD_Sharing'] >= int(four_up)
        small = mapfile['IBD_Sharing'] <= int(four_down)
        all_outliers_map = mapfile[big | small]
        all_outliers_grouped = all_outliers_map.groupby(['Chromosome'])
        chroms = list(all_outliers_grouped.groups)
        output = open(prefix + '_QC.ilash', 'w')
        ilash_grouped = ilash_file.groupby(['CHR'])
        for chrom, group in ilash_grouped:
            if chrom not in chroms:
                chrom_list = ilash_grouped.get_group(chrom).reset_index().values.tolist()
                for line in range(len(chrom_list)):
                    chrom_list[line].pop(0)
                    chrom_list[line] = [str(i) for i in chrom_list[line]]
                    output.write('\t'.join(chrom_list[line]) + '\n')
            else:
                all_outliers_list = all_outliers_grouped['Physical_Position'].get_group(chrom).reset_index(drop=True).tolist()
                group = group[~group['START'].isin(all_outliers_list)] 
                group = group[~group['START'].isin(all_outliers_list)]
                group_list = group.reset_index(drop=True).values.tolist()
                for group in range(len(group_list)):
                    yes = 0
                    for position in all_outliers_list:
                        if int(group_list[group][5]) < int(position) < int(group_list[group][6]):
                            yes +=1
                    if yes == 0:
                        group_list[group] = [str(i) for i in group_list[group]]
                        output.write('\t'.join(group_list[group]) +'\n')
        output.close()
                
                
        
            
"""
mapfile = IBDDepth.load_files(sys.argv[1])
ilash_file = IBDDepth.load_ilash(sys.argv[2])
mapfile = IBDDepth.find_ibd(mapfile, ilash_file)
four_up, four_down = IBDDepth.get_stdev(mapfile)
print(four_up, four_down)
IBDDepth.make_depth_plot(mapfile, sys.argv[3], four_up, four_down)
    

#ibd_counts = IBDDepth.map_to_dict(mapfile)
#chrom_dict = IBDDepth.make_chrom_dict(ibd_counts)

#chrom_dict = IBDDepth.fill_depth_dict(sys.argv[2], ibd_counts, chrom_dict)
#print(chrom_dict)


ax_title = 'Chromosome ' + str(chromosome_number)
sorted_depth = sorted(ibd_counts.items())

x, y = zip(*sorted_depth)
scale_x = tuple(xi/1000000 for xi in x)
fig, ax = plt.subplots()
ax.scatter(scale_x, y, c = 'xkcd:black')
ax.set_title(ax_title)
ax.set_xlabel('Physical Position (Mb)')
ax.set_ylabel('Depth of Pairwise IBD Sharing')
plt.savefig(sys.argv[3])
"""
        




