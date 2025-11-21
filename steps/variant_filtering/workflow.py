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

vcf_all = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_calling/Sep25/outputs/sorted_Sep25.vcf.gz'
# this time am using the vcf with only scaffolds > 1000bp, will call this mtt (more than thousand)

# get stats on subset to see what filtering to do
stat_subset = gwf.target_from_template(
    name = 'get_stats_initial',
    template=vcf_stats_subset(
        vcf_file=vcf_all,
        sampling_frq=0.01,
        prefix='pre_filter.mtt',
        out_dir=stat_dir,
        tmp_dir=temp_dir
        )
    )
# from this I should get the average DP and QUAL to use in the next steps
# average depth is around 30, quality is abyssmal thus I won't filter on it now (also because by design invariant sites have 0 quality and I want to keep them)
avg_DP = 30
# Because one sample is less deep but overall fine, I have less strict minDP and maxDP filters.
min_ind_DP = 3
max_ind_DP = 80

filter_qual = gwf.target_from_template(
    name = 'filter_qual_all',
    template=quality_filter(
        vcf_file=vcf_all,
        prefix='mtt',
        out_dir=temp_dir,
        #min_qual=20,
        max_missing=0.9,
        min_DP=int(avg_DP/3),
        max_DP=int(avg_DP*2),
        min_ind_DP=min_ind_DP,
        max_ind_DP=max_ind_DP
        )
    )

# check what quality filtering did
# it won't remove much so do it on small subset again

stat_postqual = gwf.target_from_template(
    name = 'get_stats_post_qual',
    template=vcf_stats_subset(
        vcf_file=filter_qual.outputs['filtered_file'],
        sampling_frq=0.01,
        prefix='post_qualfilter.mtt',
        out_dir=stat_dir,
        tmp_dir=temp_dir
        )
    )

# for calculating Fst I also need to remove multiallelic sites and indels (but keep monomorphic sites i.e. no MAF filter)
filter_multiallelic = gwf.target_from_template(
    name = 'filter_multiallelic_all',
    template=allele_filter(
        vcf_file=filter_qual.outputs['filtered_file'],
        prefix='mtt.filtered', # what are you putting in?
        out_dir=output_dir
        )
    )
  
# check what allele filtering did
# input vcf should be gzipped
stat_postallele = gwf.target_from_template(
    name = 'get_stats_post_allele',
    template=vcf_stats_subset(
        vcf_file=filter_multiallelic.outputs['filtered_file'],
        sampling_frq=0.01,
        prefix='post_allelefilter.mtt',
        out_dir=stat_dir,
        tmp_dir=temp_dir
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
        sampling_frq=0.05,
        prefix='post_biasfilter',
        out_dir=stat_dir
        )
    )
 """



# split file into populations for Hardy-Weinberg filter
# but without the mixed individuals they will mess it up (133 and 122)
blacksamples = "P1878_107,P1878_108,P1878_109,P1878_110,P1878_111,P1878_112,P1878_113,P1878_114,P1878_115,P1878_116"
bluesamples = "P1878_117,P1878_118,P1878_119,P1878_120,P1878_121,P1878_123,P1878_124,P1878_125,P1878_126"
redsamples = "P1878_127,P1878_128,P1878_129,P1878_130,P1878_131,P1878_132,P1878_134,P1878_135,P1878_136"

HWE_black = gwf.target_from_template(
    name = 'black_HWE',
    template=HWE_filter(
        vcf_file=filter_multiallelic.outputs['filtered_file'],
        prefix='black',
        samples_list=blacksamples,
        temp_dir=temp_dir,
        out_dir=output_dir
        )
    )

HWE_blue = gwf.target_from_template(
    name = 'blue_HWE',
    template=HWE_filter(
        vcf_file=filter_multiallelic.outputs['filtered_file'],
        prefix='blue',
        samples_list=bluesamples,
        temp_dir=temp_dir,
        out_dir=output_dir
        )
    )

HWE_red = gwf.target_from_template(
    name = 'red_HWE',
    template=HWE_filter(
        vcf_file=filter_multiallelic.outputs['filtered_file'],
        prefix='red',
        samples_list=redsamples,
        temp_dir=temp_dir,
        out_dir=output_dir
        )
    )

# in the end do stats again for each subspecies
stat_HWE_black = gwf.target_from_template(
    name = 'black_HWE_stats',
    template=vcf_stats_subset(
        vcf_file=HWE_black.outputs['filtered_file'],
        sampling_frq=0.02,
        prefix='Black_post_HWE',
        out_dir=stat_dir,
        tmp_dir=temp_dir
        )
    )

stat_HWE_blue = gwf.target_from_template(
    name = 'blue_HWE_stats',
    template=vcf_stats_subset(
        vcf_file=HWE_blue.outputs['filtered_file'],
        sampling_frq=0.02,
        prefix='Blue_post_HWE',
        out_dir=stat_dir,
        tmp_dir=temp_dir
        )
    )

stat_HWE_red = gwf.target_from_template(
    name = 'red_HWE_stats',
    template=vcf_stats_subset(
        vcf_file=HWE_red.outputs['filtered_file'],
        sampling_frq=0.02,
        prefix='Red_post_HWE',
        out_dir=stat_dir,
        tmp_dir=temp_dir
        )
    )

vcf_list = {HWE_black.outputs['filtered_file'], HWE_blue.outputs['filtered_file'], HWE_red.outputs['filtered_file']}
with open(f'{temp_dir}/vcf_list.txt', 'w') as f:
    for vcf in vcf_list:
        f.write(f"{vcf}\n")
vcf_file = f'{temp_dir}/vcf_list.txt'

# merge vcf files from different populations
merge_vcfs = gwf.target_from_template( 
    name = 'merge_pop_vcfs',
    template=merge_samples(
        vcf_listfile=vcf_file,
        vcf1=HWE_black.outputs['filtered_file'], # because I forgot the indexing earlier
        vcf2=HWE_blue.outputs['filtered_file'],
        vcf3=HWE_red.outputs['filtered_file'],
        prefix='allpops.HWE',
        out_dir=output_dir,
        )
    )


## find out what sites I removed
# what did depth/quality filtering reomove?
out_sites_depth = gwf.target_from_template(
    name = 'find_baddepth_sites',
    template=removed_sites(
        original_vcf=vcf_all,
        filtered_vcf_file=filter_qual.outputs['filtered_file'],
        prefix='qual_depth',
        out_dir=output_dir,
        )
    )
# what did allele filtering remove?
out_sites_allele = gwf.target_from_template(
    name = 'find_indel_multi_sites',
    template=removed_sites(
        original_vcf=vcf_all,
        filtered_vcf_file=filter_multiallelic.outputs['filtered_file'],
        prefix='indel_multiallelic',
        out_dir=output_dir,
        )
    )
# here I kept the original vcf file as comparison because I want to see all indels and multiallelics
# not just the good quality/depth ones

# what did HWE filtering remove?
# for each subspecies
out_sites_HWE_black = gwf.target_from_template(
    name = 'find_black_HWE_sites',
    template=removed_sites(
        original_vcf=filter_multiallelic.outputs['filtered_file'],
        filtered_vcf_file=HWE_black.outputs['filtered_file'],
        prefix='black_HWE',
        out_dir=output_dir,
        )
    )
out_sites_HWE_blue = gwf.target_from_template(
    name = 'find_blue_HWE_sites',
    template=removed_sites(
        original_vcf=filter_multiallelic.outputs['filtered_file'],
        filtered_vcf_file=HWE_blue.outputs['filtered_file'],
        prefix='blue_HWE',
        out_dir=output_dir,
        )
    )
out_sites_HWE_red = gwf.target_from_template(
    name = 'find_red_HWE_sites',
    template=removed_sites(
        original_vcf=filter_multiallelic.outputs['filtered_file'],
        filtered_vcf_file=HWE_red.outputs['filtered_file'],
        prefix='red_HWE',
        out_dir=output_dir,
        )
    )