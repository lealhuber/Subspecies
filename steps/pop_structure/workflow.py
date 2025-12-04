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

vcf_qualfiltered = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/outputs/mtt.filtered.monobiallelic.snp.vcf.gz'

# get stats on the vcf before SNP filtering because I forgot on merged hwe
""" stat_preSNP = gwf.target_from_template(
    name = 'get_stats_HWEmerge',
    template=vcf_stats_subset(
        vcf_file=vcf_qualfiltered,
        sampling_frq=0.01,
        prefix='HWE.merged',
        out_dir=stat_dir,
        temp_dir=temp_dir
        )
    ) """

# make list of individuals to remove
remove_ind = {'P1878_117','P1878_122', 'P1878_129', 'P1878_133', 'P1878_134','P1878_136'} # 117, 129, 134, 136 are related based on previous analysis
# mixes 122 and 133 I already removed during HWE filtering
with open (f'{log_dir}/remove_ind.txt', "w") as f:
    for ind in remove_ind:
        f.write(f'{ind}\n')

# filter for SNPs of sufficient quality and MAC, and remove related individuals
# mac of 1 means remove singletons, I think that's better than doing maf because of varying sample sizes per pop
filter_snp = gwf.target_from_template(
    name = 'filter_snps',
    template=snp_filter(
        vcf_file=vcf_qualfiltered,
        ind_file=f'{log_dir}/remove_ind.txt',
        MAC=1,
        minQ=30,
        prefix='indfiltered.qualfiltered',
        out_dir=temp_dir
        )
    )

# check what snp filtering did
stat_postSNP = gwf.target_from_template(
    name = 'get_stats_SNPfilter',
    template=vcf_stats_subset(
        vcf_file=filter_snp.outputs['filtered_vcf'],
        sampling_frq=0.1,
        prefix='snps',
        out_dir=stat_dir,
        temp_dir=temp_dir
        )
    )


# do LD pruning and PCA
ld_pca = gwf.target_from_template(
    name = 'PCA_relrem',
    template = PCA(
        vcf_file=filter_snp.outputs['filtered_vcf'],
        prefix='pca_noRelnoMix',
        temp_dir=temp_dir,
        out_dir=output_dir
        )
    )

# run ADMIXTURE for K=2-5
Admixture = gwf.target_from_template(
    name = 'Admixture',
    template = admixture(
        bed_file=ld_pca.outputs['Plink_bed'],
        bim_file=ld_pca.outputs['Plink_bim'],
        prefix='admixture_noRelnoMix',
        tmp_dir=temp_dir
    )
)

