## declare the workflow
from gwf import Workflow # type: ignore
import numpy as np # type: ignore
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

vcf_all = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_calling/Sc_all/outputs/sorted_Aug25.vcf.gz'

stat_nrs = gwf.target_from_template(
    name = 'get_vcf_stats',
    template=vcf_stats(
        vcf_file=vcf_all,
        prefix='pre_filter',
        out_dir=stat_dir
        )
    )

# from this I should get the average DP and QUAL to use in the next steps
# I will just put in some values for now, but these should be calculated from the vcf stats output
avg_DP_ind = 20
avg_DP = 20
avg_qual = 200

filter_qual = gwf.target_from_template(
    name = 'filter_qual_all',
    template=quality_filter(
        vcf_file='/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_calling/Sc_all/outputs/all.vcf',
        prefix='all',
        out_dir=temp_dir,
        min_site_DP=int(avg_DP/3),
        max_site_DP=int(avg_DP*2),
        min_SNP_DP=int(avg_DP_ind/3),
        max_SNP_DP=int(avg_DP_ind*3)
        )
    )

filter_bias = gwf.target_from_template(
    name = 'filter_bias_all',
    template=bias_filter(
        vcf_file=filter_qual.outputs['filtered_file'],
        prefix='all',
        out_dir=temp_dir,
        min_qual=20,   
        avgDP=avg_DP,
        qual_div_avgDP=avg_DP/20
        )
    )


# in the end do stats again
stat_nrs = gwf.target_from_template(
    name = 'get_vcf_stats',
    template=vcf_stats(
        vcf_file=vcf_all,
        prefix='post_filter',
        out_dir=stat_dir
        )
    )