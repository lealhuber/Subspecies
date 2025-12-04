## declare the workflow
from gwf import Workflow # type: ignore
import re, os
from templates import *

gwf = Workflow(defaults={'account': 'ostrich_thermal'})

# set directories for outputs

log_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/targeted_fst/logs'
stat_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/targeted_fst/stats'
temp_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/targeted_fst/temp'
output_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/targeted_fst/outputs'

# Check if the directories exists, and create it if not
os.makedirs(log_dir, exist_ok=True)
os.makedirs(stat_dir, exist_ok=True)
os.makedirs(temp_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

ostrich_gtf = '/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/Struthio_camelus_HiC_augustus.gtf'


### get regions of interest
genesOI = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/targeted_fst/genes_of_interest.txt'
# make file with columns scaffold, start, end from gtf and genesOI
gene_regions = gwf.target_from_template(
    name = 'get_gene_regions',
    template=get_regions(
        gtf=ostrich_gtf,
        genesOI=genesOI,
        add_bp=1000,
        prefix="regionsOI",
        output_dir=temp_dir
        )
    )

### calculate population statistics for regions of interest
vcf_file = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/outputs/allpops.HWE.merged.vcf.gz'
pops_file = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/targeted_fst/populations_file.txt'
script_path = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/fst/scripts/simonmartin_genomics'
# make geno file name on basis of vcf prefix
geno_name = re.sub(r'.*/', '', vcf_file)  # remove path
geno_name = re.sub(r'\.vcf\.gz$', '', geno_name)  # remove .vcf.gz
print(f'geno_name: {geno_name}') # check
print(f'regions file: {gene_regions.outputs["regions_file"]}') # check

pop_stats_targeted = gwf.target_from_template(
    name = 'pop_stats_targeted_fst',
    template=pop_stats(
        vcf_file=vcf_file,
        snp_vcf='/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/pop_structure/temp/indfiltered.qualfiltered.HWE.biallelic.snp.vcf.gz',
        script_path=script_path,
        regions_file=gene_regions.outputs["regions_file"],
        pops_file=pops_file,
        prefix='monobiallelic.HWE',
        geno_prefix=geno_name,
        temp_dir=temp_dir,
        output_dir=output_dir
        )
    )