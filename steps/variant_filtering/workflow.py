## declare the workflow
from gwf import Workflow # type: ignore
import re, os
from templates import *

gwf = Workflow(defaults={'account': 'ostrich_thermal'})

# set directories for outputs

log_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/logs'
stat_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/stats'
temp_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/temp'
output_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/outputs'

# Check if the directories exists, and create it if not
os.makedirs(log_dir, exist_ok=True)
os.makedirs(stat_dir, exist_ok=True)
os.makedirs(temp_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

vcf_all = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_calling/Aug25/outputs/sorted_Aug25.vcf.gz'

# make random subset of vcf for testing
vcf_sub = gwf.target_from_template(
    name = 'random_subset',
    template = random_subset(
        vcf_file=vcf_all,
        prefix='raw',
        out_dir=temp_dir
        )
    )

stat_subset = gwf.target_from_template(
    name = 'get_stats_subset',
    template=vcf_stats(
        vcf_file=vcf_sub.outputs['subset_file'],
        prefix='pre_filter.subset',
        out_dir=stat_dir
        )
    )
# from this I should get the average DP and QUAL to use in the next steps
# average depth is around 30, quality is abyssmal thus I won't filter on it now (also because by design invariant sites have 0 quality and I want to keep them)
avg_DP = 30

filter_qual = gwf.target_from_template(
    name = 'filter_qual_all',
    template=quality_filter(
        vcf_file=vcf_all,
        prefix='all',
        out_dir=temp_dir,
        #min_qual=20,
        max_missing=0.9,
        min_DP=int(avg_DP/3),
        max_DP=int(avg_DP*2),
        )
    )

# check what quality filtering did
# it won't remove much so do it on small subset again

vcf_sub2 = gwf.target_from_template(
    name = 'subset_filtered',
    template = random_subset(
        vcf_file=filter_qual.outputs['filtered_file'],
        prefix='filtered',
        out_dir=temp_dir
        )
    )

stat_filtered = gwf.target_from_template(
    name = 'get_filtered_stats',
    template=vcf_stats(
        vcf_file=vcf_sub2.outputs['subset_file'],
        prefix='qualfilter_sub',
        out_dir=stat_dir
        )
    )

# for calculating Fst I also need to remove multiallelic sites and indels (but keep monomorphic sites i.e. no MAF filter)
filter_multiallelic = gwf.target_from_template(
    name = 'filter_multiallelic_all',
    template=allele_filter(
        vcf_file=filter_qual.outputs['filtered_file'],
        prefix='all.filtered', # what are you putting in?
        out_dir=output_dir
        )
    )
  
# check what allele filtering did
stat_postallele = gwf.target_from_template(
    name = 'get_stats_post_allele',
    template=vcf_stats_subset(
        vcf_file=filter_multiallelic.outputs['filtered_file'],
        sampling_frq=0.5,
        prefix='post_allelefilter',
        out_dir=stat_dir
        )
    )

 # One can also filter by some biases, I'll try but not sure if I will use it
 # My vcf files don't contain the relevant info fields though, so nevermind for now

""" filter_bias = gwf.target_from_template(
    name = 'filter_bias_all',
    template=bias_filter(
        vcf_file=filter_qual.outputs['filtered_file'],
        prefix='all.filtered', # what are you putting in?
        out_dir=temp_dir,
        #min_qual=20,   
        avgDP=avg_DP,
        qual_div_avgDP=avg_DP/20
        )
    ) """

# check what bias filtering did
""" stat_postbias = gwf.target_from_template(
    name = 'get_stats_post_bias',
    template=vcf_stats_subset(
        vcf_file=filter_bias.outputs['filtered_file'],
        sampling_frq=0.5,
        prefix='post_biasfilter',
        out_dir=stat_dir
        )
    )
 """



# split file into populations for Hardy-Weinberg filter
blacksamples = {'P1878_107', 'P1878_108', 'P1878_109', 'P1878_110', 'P1878_111', 'P1878_112', 'P1878_113', 'P1878_114', 'P1878_115', 'P1878_116'}
bluesamples = {'P1878_117', 'P1878_118', 'P1878_119', 'P1878_120', 'P1878_121', 'P1878_122', 'P1878_123', 'P1878_124', 'P1878_125', 'P1878_126'}
redsamples = {'P1878_127', 'P1878_128', 'P1878_129', 'P1878_130', 'P1878_131', 'P1878_132', 'P1878_133', 'P1878_134', 'P1878_135', 'P1878_136'}

""" HWE_black = gwf.target_from_template(
    name = 'black_HWE,
    template=HWE_filter(
        vcf_file=filter_biallelic.outputs['filtered_file'],
        prefix='black.biallelic.snp',
        samples_list=blacksamples
        temp_dir=temp_dir,
        out_dir=output_dir,
        stat_dir=stat_dir
    )
) """

""" HWE_blue = gwf.target_from_template(
    name = 'blue_HWE,
    template=HWE_filter(
        vcf_file=filter_biallelic.outputs['filtered_file'],
        prefix='blue.biallelic.snp',
        samples_list=bluesamples
        temp_dir=temp_dir,
        out_dir=output_dir,
        stat_dir=stat_dir
    )
) """

""" HWE_red = gwf.target_from_template(
    name = 'red_HWE,
    template=HWE_filter(
        vcf_file=filter_biallelic.outputs['filtered_file'],
        prefix='red.biallelic.snp',
        samples_list=redsamples
        temp_dir=temp_dir,
        out_dir=output_dir,
        stat_dir=stat_dir
    )
) """

# in the end do stats again
""" stat_nrs = gwf.target_from_template(
    name = 'get_vcf_stats',
    template=vcf_stats(
        vcf_file=HWE_black.outputs['filtered_file'],
        prefix='Black_post_HWE',
        out_dir=stat_dir
        )
    ) """