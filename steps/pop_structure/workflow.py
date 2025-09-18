## declare the workflow
from gwf import Workflow # type: ignore
import re, os
from templates import *

gwf = Workflow(defaults={'account': 'ostrich_thermal'})

# set directories for outputs

log_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/pop_structure/logs'
stat_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/pop_structure/stats'
temp_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/pop_structure/temp'
output_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/pop_structure/outputs'

# Check if the directories exists, and create it if not
os.makedirs(log_dir, exist_ok=True)
os.makedirs(stat_dir, exist_ok=True)
os.makedirs(temp_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

vcf_filtered = 'path'

# do LD pruning and PCA
LD_PCA = gwf.target_from_template(
    name = 'PCA_all',
    template = PCA(
        vcf_file=vcf_filtered,
        prefix='pca_all',
        temp_dir=temp_dir,
        out_dir=output_dir
        )
    )

