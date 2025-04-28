## declare the workflow
from gwf import Workflow # type: ignore
import numpy as np # type: ignore
import re
from templates import *

gwf = Workflow(defaults={'account': 'ostrich_thermal'})

# set directories for outputs

log_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/logs'
stat_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/stats'
temp_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/temp'
output_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/outputs'

filter_qual = gwf.target_from_template(
    name = 'filter_qual_all',
    template=quality_filter(
        vcf_file='/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_calling/Sc_all/outputs/all.vcf',
        prefix='all',
        out_dir=temp_dir
        )
    )