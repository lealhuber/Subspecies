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
genesOI = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/targeted_fst/genesOI_ageint.txt'
# make file with columns scaffold, start, end from gtf and genesOI
gene_regions = gwf.target_from_template(
    name = 'get_gene_regions',
    template=get_regions_new(
        gtf=ostrich_gtf,
        genesOI=genesOI,
        add_bp=500,
        prefix="regionsOI_ageint_500bp",
        output_dir=temp_dir
        )
    )

## for comparison, get all coding regions
all_coding_regions = gwf.target_from_template(
    name = 'get_coding_regions',
    template=get_coding_regions(
        gtf=ostrich_gtf,
        prefix="Struthio_camelus_HiC",
        output_dir=temp_dir
        )
    )

### calculate population statistics for regions of interest
vcf_file = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/outputs/mtt.filtered.monobiallelic.snp.vcf.gz'
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
        snp_vcf='/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/pop_structure/temp/indfiltered.qualfiltered.biallelic.snp.vcf.gz', # should match vcf file
        script_path=script_path,
        regions_file=all_coding_regions.outputs["coding_regions_file"],
        pops_file=pops_file,
        prefix='monobiallelic.allCDS',
        geno_prefix=geno_name,
        temp_dir=temp_dir,
        output_dir=output_dir
        )
    )

### calculate Tajima's D for regions of interest
# for blue and red separately
tajimasD_targeted_blue = gwf.target_from_template(
    name = 'tajimasD_targeted_blue',
    template=tajimas_d(
        vcf_file=vcf_file,
        regions_file=gene_regions.outputs["regions_file"],
        ind_file='/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/targeted_fst/blue_inds.txt',
        prefix='monobiallelic.500bp_start.ageint.blue3',
        temp_dir=temp_dir,
        out_dir=output_dir
        )
    )

tajimasD_targeted_red = gwf.target_from_template(
    name = 'tajimasD_targeted_red',
    template=tajimas_d(
        vcf_file=vcf_file,
        regions_file=gene_regions.outputs["regions_file"],
        ind_file='/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/targeted_fst/red_inds.txt',
        prefix='monobiallelic.500bp_start.ageint.red3',
        temp_dir=temp_dir,
        out_dir=output_dir
        )
    )

pixy_targeted = gwf.target_from_template(
    name = 'pixy_targeted',
    template=pixy_stats(
        vcf_file=vcf_file,
        regions_file=gene_regions.outputs["regions_file"],
        pops_file=pops_file,
        prefix='monobiallelic.500bp_start.ageint',
        out_dir=output_dir
        )
    )

### calculate population statistics for all coding regions (comparison)
pixy_stats_allCDS = gwf.target_from_template(
    name = 'pixy_stats_allCDS',
    template=pixy_stats(
        vcf_file=vcf_file,
        regions_file=all_coding_regions.outputs["coding_regions_file"],
        pops_file=pops_file,
        prefix='mtt.monobiallelic.allCDS',
        out_dir=output_dir
        )
    )