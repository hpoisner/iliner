#!/bin/usr/python3

import argparse
import logging
import numpy as np
import sys
from .iliner_load_parse import LoadFiles, FileParser
from .iliner_make_files import FileConvert
from .iliner_launch_ilash import RunIlash
from .iliner_ibd_depth import IBDDepth
from .iliner_stats import IBDStats
from .iliner_parse_vcf import VcfReader, VcfParser, VcfMakeFiles


class LogFilter(logging.Filter):
    """Code from https://stackoverflow.com/questions/1383254/logging-streamhandler-and-standard-streams
    """
    def __init__(self, level):
        self.level = level

    def filter(self, record):
        return record.levelno < self.level

class Logger(object):
    """Creates a log instance that can be used throughout the program
    """
    def __init__(self):
        """Creates a logging instance for run of code
        """
        ## Creating the looger
        self.logger = logging.getLogger('iLiner_Logger')
        ## Setting the level for the logger
        self.logger.setLevel(logging.DEBUG)
        ## Creating the handler
        stdout_handler = logging.StreamHandler(sys.stdout)
        ## Creating the formatter
        formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
        stdout_handler.setFormatter(formatter)
        stdout_handler.setLevel(logging.DEBUG)
        self.logger.addHandler(stdout_handler)




def main():
    """
    Argument Parser Function that calls the other components of the code
    """
    parser = argparse.ArgumentParser(description = 'Arguments for running the different modules of iLiner')

    parser.add_argument('module', choices=['ilash', 'ibd_depth', 'stats'], help='module choice')
    parser.add_argument('--sample', '-s', type=str, required=False, help='Sample file with file path')
    parser.add_argument('--haps', '-hp', type=str, required=False, help='Phased haplotype file with file path')
    parser.add_argument('--genetic_map', '-gm', type=str, required=False, help='Genetic map file with file path')
    parser.add_argument('--mapfile', '-mf', type=str, required=False, help='Mapfile made by the ilash module')
    parser.add_argument('--outfile_prefix', '-op', type=str, help='Prefix for all of the files produced from this module with file path')
    parser.add_argument('--ilash', '-i', type=str, required=False, help='File path to ilash')
    parser.add_argument('--ilash_output', '-io', type=str, required=False, help='Ilash output, one chromosome for ibd_depth module and all chromosomes for stats module')
    parser.add_argument('--population_file', '-pf', type=str, required=False, help='File with individual ids and population [ID]tab[Population]')
    parser.add_argument('--vcf', '-v', type=str, required=False, help='Phased VCF with file path')
    parser.add_argument('--no_indel', '-ni', action='store_true', required=False, help='Pass if you would like to remove indels from your vcf or haplotype file')

    args = parser.parse_args()
    log = Logger()
    if args.module == 'ilash':
        args_dict = vars(args)
        log.logger.info('You have selected ilash')
        if args_dict['sample'] != None:
            log.logger.info("Parsing your sample and haplotype files")
            samp_file = LoadFiles.load_sample(args.sample)
            log.logger.info('Sample file has been loaded')
            haps_file = LoadFiles.load_haps(args.haps)
            log.logger.info('Haplotype file has been loaded')
            ilash_obj = LoadFiles(samp_file, haps_file)
            FileParser.check_sample_file(ilash_obj)
            log.logger.info('Your sample file has been validated')
            FileParser.check_same_sample(ilash_obj)
            log.logger.info('Your sample file and haplotype file populations are the same')
            if args_dict['no_indel'] == True:
                FileParser.remove_indels(ilash_obj)
            else:
                pass
            FileParser.check_alleles(ilash_obj.haps_file.iloc[:,[3]].values, ilash_obj.haps_file.iloc[:,[4]].values)
            nucfunc = np.vectorize(FileParser.nucleotide_test)
            log_array_ref = nucfunc(ilash_obj.haps_file[[3]])
            log_array_alt = nucfunc(ilash_obj.haps_file[[4]])
            FileParser.validate_nucleotides(log_array_ref, 'Reference Allele')
            FileParser.validate_nucleotides(log_array_alt, 'Alternative Allele')
            log.logger.info('Your reference and alternative alleles have been validated')
            FileParser.validate_binary(ilash_obj)
            log.logger.info('Your haplotypes have been validated')
            log.logger.info("Making your map and pedigree files")
            pedfile_start = FileConvert.start_pedfile(ilash_obj)
            pedlen = int(len(pedfile_start))
            hapmap = open(args.genetic_map)
            mapfile_str = args.outfile_prefix + ".map"
            pedfile_str = args.outfile_prefix + ".ped"
            mapfile = open(mapfile_str, 'w')
            FileConvert.make_gpos_mapfile(ilash_obj, hapmap, mapfile)
            log.logger.info('Your map file has been made')
            hapslist = FileConvert.convert_haps(ilash_obj)
            final_array = FileConvert.make_pedfile(hapslist, pedlen, pedfile_start)
            final_list = final_array.tolist()
            pedfile = open(pedfile_str, 'w')
            for line in range(len(final_list)):
                pedfile.write(' '.join(final_list[line]) + '\n')
            pedfile.close()
            log.logger.info('Your pedigree file has been made')
            log.logger.info('Launching iLASH')
            RunIlash.make_param_file(mapfile_str, pedfile_str, args.ilash)
        if args_dict['vcf'] != None:
            log.logger.info('Your VCF is being parsed')
            vcf_file, comments = VcfReader.read_vcf(args.vcf)
            vcf_obj = VcfReader(vcf_file, comments)
            vcf_header = vcf_file.columns.values.tolist()
            VcfParser.validate_header(vcf_header)
            VcfParser.validate_chrom(vcf_obj)
            if args_dict['no_indels'] == True:
                VcfParser.remove_indels(vcf_ojb)
            else:
                pass
            VcfParser.diff_alt_ref(vcf_obj.vcf['REF'], vcf_obj.vcf['ALT'])
            log.logger.info('Your VCF has been parsed')
            nucfunc = np.vectorize(VcfParser.nucleotide_test)
            log_array_ref = nucfunc(vcf_obj.vcf['REF'])
            log_array_alt = nucfunc(vcf_obj.vcf['ALT'])
            VcfParser.validate_nucleotides(vcf_obj, log_array_ref, 'Reference Allele')
            VcfParser.validate_nucleotides(vcf_obj, log_array_alt, 'Alternative Allele')
            binary_array = np.isin(vcf_obj.vcf.iloc[:,9:], ['0|0', '0|1', '1|0', '1|1'])
            VcfParser.validate_binary(vcf_obj, binary_array)
            people = VcfMakeFiles.get_people(vcf_header)
            pedfile_start = VcfMakeFiles.start_pedfile(people)
            hapslist = VcfMakeFiles.convert_haps(vcf_obj)
            pedlen = int(len(haplist[0]/2))
            final_array = VcfMakeFiles.make_pedfile(hapslist, pedlen, pedfile_start)
            final_list = final_array.tolist()
            mapfile_str = args.outfile_prefix + ".map"
            pedfile_str = args.outfile_prefix + ".ped"
            mapfile = open(mapfile_str, 'w')
            pedfile = open(pedfile_str, 'w')
            hapmap = open(args.genetic_map)
            VcfMakeFiles.make_gpos_mapfile(vcf_obj, hapmap, mapfile)
            log.logger.info('Your map file has been made')
            for line in range(len(final_list)):
                pedfile.write(' '.join(final_list[line]) + '\n')
            pedfile.close()
            log.logger.info('Your pedigree file has been made')
            log.logger.info('Launching iLASH')
            RunIlash.make_param_file(mapfile_str, pedfile_str, args.ilash)
    elif args.module == 'ibd_depth':
        log.logger.info('You have selected ibd_depth')
        mapfile = IBDDepth.load_files(args.mapfile)
        log.logger.info('Your map file has been loaded')
        ilash_output, snps = IBDDepth.load_ilash(args.ilash_output)
        log.logger.info('Your ilash output file has been loaded')
        mapfile = IBDDepth.find_ibd(mapfile, snps)
        four_up, four_down = IBDDepth.get_stdev(mapfile)
        IBDDepth.make_depth_plot(mapfile, args.outfile_prefix, four_up, four_down)
        log.logger.info('Your IBD depth figures have been made')
        log.logger.info('Removing outliers and creating an outlier free file')
        IBDDepth.remove_outlier(mapfile, ilash_output, args.outfile_prefix, four_up, four_down)
    elif args.module == 'stats':
        log.logger.info('You have selected stats')
        ilash_results = IBDStats.load_ilash(args.ilash_output)
        log.logger.info('Your iLASH output has been loaded')
        ilash_pops = IBDStats.load_popfile(args.population_file)
        log.logger.info('Your population file has been loaded')
        stats_obj = IBDStats(ilash_results, ilash_pops)
        IBDStats.remove_self_ibd(stats_obj)
        IBDStats.stratify_by_pop(stats_obj)
        IBDStats.make_dot_plot(stats_obj, args.outfile_prefix)
        log.logger.info('Your pair plot has been made')
        IBDStats.make_violin_plot(stats_obj, args.outfile_prefix)
        log.logger.info('Your violin plot has been made')
        IBDStats.get_ibd_stats(stats_obj, log)
        H, pval = (IBDStats.get_kw_h(stats_obj))
        log.logger.info('Kruskal-Wallis Test')
        log.logger.info('H-statistic: ' + str(H))
        log.logger.info('P-Values: ' + str(pval))
        if pval < 0.05:
            log.logger.info('Reject the null hypothesis - A significant differences exists between groups.')
            log.logger.info('Post-hoc Wilcoxon Rank Sum Test')
            IBDStats.get_ranksum(stats_obj, log)
        else:
            log.logger.info('Fail to reject the hypothesis - A discernable difference does not exist between the groups')
        pop_heatmap = IBDStats.fraction_ibd_sharing(stats_obj)
        fig = pop_heatmap.get_figure()
        fig.set_size_inches(13.0, 13.0)
        fig.savefig(args.outfile_prefix + '_heatmap.png')
        log.logger.info('Your heatmap has been made')


if __name__ == '__main__':
    main()
