#!/bin/usr/python

from collections import defaultdict
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import seaborn as sns
import sys
import os

class IBDStats(object):

    def __init__(self, ilash, population,grouped_ilash=None, stratified_ilash=None):
        self.ilash = ilash
        self.population = population
        self.grouped_ilash = grouped_ilash
        self.stratified_ilash = stratified_ilash

    @staticmethod
    def load_ilash(path):
        """

        Args:
            path(:obj:`filepath`): file path for ilash (should be tab delimited)

        Raises:
            :obj:`IOError`: if file provided does not exist

        Returns:
            :obj:`pandas dataframe`: file as a pd df
        """
        if os.path.isfile(path):
            ilash_file = pd.read_table(path, sep='\t', header=None, usecols=[0,2,9], names = ['ID_1', 'ID_2', 'IBD_Sharing'])
            return ilash_file 
        else:
            raise IOError("{} does not exist.".format(path))

    @staticmethod
    def load_popfile(path):
        """
        Args:
        path(:obj:`filepath`): file path for population file

        Raises:
        :obj:`IOError`: if file provided does not exist
        
        Returns:
        :obj:`pandas dataframe`: file as a pd df
        """
        if os.path.isfile(path):
            pop_file = pd.read_table(path, sep='\t', header=None, names=['Person', 'Population'])
            return pop_file
        else:
            raise IOError("{} does not exist.".format(path))
    

    def remove_self_ibd(self):
        """
        
        Args:
            self object(:obj:`self`): self object with attributes
           
        Modifies:
            :attribute:`ilash`: Removes instances of IBD with self
        """
        self.ilash = self.ilash[self.ilash['ID_1'] != self.ilash['ID_2']]
        

    def stratify_by_pop(self):
        """
        
        Args:
            self object(:obj:`self`): self object with attributes 

        Modifies:
            :attribute:`stratified_ilash`: Creates a dataframe with people stratified by populaton with IBD Sum and Tract Counts
        """
        pop_index = self.population.set_index('Person')
        self.ilash = pd.merge(self.ilash, pop_index, how='left', left_on='ID_1', right_index=True)
        self.ilash = pd.merge(self.ilash, pop_index, how='left', left_on='ID_2', right_index=True)
        self.grouped_ilash = self.ilash.groupby(['Population_x', 'Population_y', 'ID_1', 'ID_2'])['IBD_Sharing'].agg(['sum', 'count'])
        reset_grouped = self.grouped_ilash.reset_index()
        self.stratified_ilash = reset_grouped[reset_grouped['Population_x'] == reset_grouped['Population_y']]

                                              

    def make_dot_plot(self, prefix):
        """

        Args:
            self object(:obj:`self`): self object with attributes
        
        Returns:
           :plot:`Dot Plot`: x: IBD sum per pair, y: # of segments shared 
        """
        sns.set_context('poster')
        pair_plot = sns.pairplot(x_vars = ['sum'], y_vars = ['count'], data = self.stratified_ilash, hue='Population_x', palette = 'viridis', size=6)
        pair_plot.set(xlabel='Pairwise IBD Sharing (cM)', ylabel='Pairwise IBD Tract Count', xscale='log')
        pair_plot.savefig(prefix + '_dot_plot.png')
        plt.clf()
        
    def make_violin_plot(self, prefix):
        """

        Args:
           self object(:obj:`self`): self object with attributes
           prefix(:obj:`filepath + prefix string): filepath with prefix string for output

        Returns:
            A violin plot by population
        """
        sns.set_context('poster')
        g = sns.violinplot('Population_x', 'sum', data = self.stratified_ilash, palette= 'viridis', size=10)
        g.set(xlabel='Populations', ylabel='Pairwise IBD Sharing', yscale='log')
        g.set_xticklabels(g.get_xticklabels(), rotation=30)
        violin_fig = g.get_figure()
        violin_fig.set_size_inches(10.0, 10.0)
        violin_fig.savefig(prefix + '_violin.png')
        plt.clf()
    

    def get_ibd_stats(self, log):
        """
        Args: 
            self object(:obj:`self`): self object with attributes

        Returns:
            :groupby:`Mean Groupby`: Mean of IBD Sharing by population
            :groupby:`Median Groupby`: Median of IBD Sharing by population
        """
        pop_mean = self.stratified_ilash.groupby(['Population_x'])['sum'].mean()
        pop_median = self.stratified_ilash.groupby(['Population_x'])['sum'].median()
        log.logger.info("IBD Sharing Mean:")
        log.logger.info(pop_mean)
        log.logger.info("IBD Sharing Median:") 
        log.logger.info(pop_median)
    
                        
    def get_kw_h(self):
        """
        
        Args:
            self object(:obj:`self`); self object with attributes

        Returns:
            Results from Kruscall-Wallis: H statistic and pvalue
        """
        values_per_group = {col_name: col for col_name, col in self.stratified_ilash.groupby('Population_x')['sum']}
        args = values_per_group.values()
        h_stat = stats.kruskal(*args)
        return h_stat

    def get_ranksum(self, log):
        """
        
        Args:
            self object(:obj:`self`): self object with attributes
        
        Returns:
            For each population a paired Wilcoxen-rank sum tests: a test statistic a two-sided p-value of the test
        """
        values_per_group = {col_name: col for col_name, col in self.stratified_ilash.groupby('Population_x')['sum']}
        populations = list(values_per_group.keys())
        populations_list = list(set([tuple(sorted([i, j])) for i in populations for j in populations]))
        for pop1, pop2 in populations_list:
            if pop1 != pop2:
                result = stats.ranksums(values_per_group[pop1], values_per_group[pop2])
                log.logger.info("Population 1: " + pop1 + " Population 2: " + pop2 + " Statistic: " + str(result.statistic) + " P-value: " + str(result.pvalue))

        
        
    def nCr(n):
        f = math.factorial
        return f(n) // f(2) //f(n-2)
    

    def fraction_ibd_sharing(self):
        """

        Args:
            self object(:obj:`self`): self object with attributes
        
        Returns:
            A heatmap with the fraction of IBD sharing by population.
        """
        bothpops = self.grouped_ilash.reset_index()
        pops = bothpops[['Population_x', 'Population_y']].apply(np.sort, axis=1)
        pops_reduced = pops.groupby(['Population_x', 'Population_y']).size().to_frame('pair_count').reset_index()
        ## Get number of people per group
        pop_count = self.population.groupby('Population')['Person'].count().reset_index()
        ## Calculate choose value
        choose_list = []
        for index, row in pop_count.iterrows():
            for index2, row2, in pop_count.iterrows():
                if index == index2:  ## Within population
                    choosevalue = IBDStats.nCr(row['Person'])
                    choose_list.append([row['Population'], row2['Population'], choosevalue])
                else:  ## Between populations
                    added_count = row['Person'] + row2['Person']
                    bothchoosevalue = IBDStats.nCr(added_count)
                    firstchoosevalue = IBDStats.nCr(row['Person'])
                    secondchoosevalue = IBDStats.nCr(row2['Person'])
                    denom = bothchoosevalue - firstchoosevalue - secondchoosevalue
                    choose_list.append([row['Population'], row2['Population'], denom])
        ## Get chose values 
        choose_df = pd.DataFrame(choose_list, columns = ['Population_x', 'Population_y', 'Choose_Values'])
        ## Merge groups
        merged_pops = pd.merge(pops_reduced, choose_df, how = 'left')
        ## Calculate fraction of sharing
        merged_pops['Fraction_Sharing'] = merged_pops['pair_count']/merged_pops['Choose_Values']
        merged_pivot = merged_pops.pivot(index='Population_x', columns='Population_y', values='Fraction_Sharing')
        mask = merged_pivot.isnull()
        sns.set_context('poster')
        heatmap = sns.heatmap(merged_pivot, mask= mask, annot=True, square=True, cmap='viridis')
        #heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation = 30)
        #heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation = 30)
        return heatmap
        
        

"""
ilash_results = IBDStats.load_files(sys.argv[1])
#print(ilash_results.head())
ilash_pops = IBDStats.load_files(sys.argv[2])
stats_obj = IBDStats(ilash_results, ilash_pops, len(ilash_pops))
print(len(stats_obj.ilash))
IBDStats.remove_self_ibd(stats_obj)
print(len(stats_obj.ilash))
IBDStats.stratify_by_pop(stats_obj)


IBDStats.ibd_by_pairs(stats_obj)
"""
"""

#dot_plot =IBDStats.make_dot_plot(stats_obj)
#dot_plot.savefig(sys.argv[3])
#violin = IBDStats.make_violin_plot(stats_obj)
#violin_fig = violin.get_figure()
#violin_fig.savefig(sys.argv[3])
#print(IBDStats.get_ibd_stats(stats_obj))
#print(IBDStats.get_kw_h(stats_obj))
#pop_pairs = IBDStats.get_pop_pairs(stats_obj)
#IBDStats.get_ranksum(stats_obj, pop_pairs)
IBDStats.fraction_ibd_sharing(stats_obj)

pop_heatmap = IBDStats.fraction_ibd_sharing(stats_obj)
fig = pop_heatmap.get_figure()
fig.set_size_inches(12.0, 13.8)
fig.savefig(sys.argv[3])
"""
