#!/bin/env python3
from gwf import AnonymousTarget # type: ignore
import os

## For analysing population structure using PCA and ADMIXTURE
# For this I only need variant sites so as a first step I will filter out invariant sites from the vcf file.
# Secondly PCA and ADMIXTURE assume independence of sites, so I will do LD pruning using plink.
# Because of this I will also remove related individuals during the filtering step.

def snp_filter(vcf_file, ind_file, MAF, minQ, prefix, out_dir):
    """Template for removing invariant sites and certain individuals from already filtered vcf file
    (multiallelic sites and indels should already have been removed and basic quality filtering done)"""
    inputs = {'vcf': vcf_file,
              'individuals': ind_file}
    outputs = {'filtered_vcf': f'{out_dir}/{prefix}.biallelic.snp.vcf.gz'}
    options={
        'cores': 1,
        'memory': '4g',
        'walltime': '03:00:00'
    }
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
    vcftools --gzvcf {vcf_file} \\
    --min-alleles 2 --maf {MAF} --minQ {minQ} \\
    --remove {ind_file} \\
    --recode --recode-INFO-all \\
    --stdout | bgzip -c > {out_dir}/{prefix}.biallelic.snp.vcf.gz
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def vcf_stats_subset(vcf_file, sampling_frq, prefix, out_dir):
    """ Template for vcf stats on a subset of the data.
    Adapted from https://speciationgenomics.github.io/filtering_vcfs/
    """
    inputs = {'file': vcf_file}
    outputs = {'allele_frequency': f'{out_dir}/{prefix}.frq',
               'depth': f'{out_dir}/{prefix}.idepth',
               'mean_depth': f'{out_dir}/{prefix}.ldepth.mean',
               'site_quality': f'{out_dir}/{prefix}.lqual',
               'missing_incv': f'{out_dir}/{prefix}.imiss',
               'missing_site': f'{out_dir}/{prefix}.lmiss',
               'heterozygosity': f'{out_dir}/{prefix}.het',
               'relatedness': f'{out_dir}/{prefix}.relatedness'}
    options={
        'cores': 1,
        'memory': '4g',
        'walltime': '2-00:00:00'
    }
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
    # make random subset
    bcftools view {vcf_file} | vcfrandomsample -r {sampling_frq} > {out_dir}/{prefix}.subset.vcf
    # compress vcf
    bgzip {out_dir}/{prefix}.subset.vcf
    # index vcf
    bcftools index {out_dir}/{prefix}.subset.vcf.gz

    # calculate allele frequency (for sites with max 2 alleles)
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf.gz --freq2 --max-alleles 2 --out {out_dir}/{prefix}
    # mean depth of coverage per individual
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf.gz --depth --out {out_dir}/{prefix}
    # mean depth per site
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf.gz --site-mean-depth --out {out_dir}/{prefix}
    # site quality
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf.gz --site-quality --out {out_dir}/{prefix}
    # genotype quality
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf.gz --extract-FORMAT-info GQ --out {out_dir}/{prefix}
    # proportion of missing data per individual
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf.gz --missing-indv --out {out_dir}/{prefix}
    # proportion of missing data per site
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf.gz --missing-site --out {out_dir}/{prefix}
    # calculate heterozygosity and inbreeding coefficient per individual
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf.gz --het --out {out_dir}/{prefix}
    # caculate relatedness
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf.gz --relatedness --out {out_dir}/{prefix}
    # Count the number of  variants (takes for ever!)
    bcftools view -H {vcf_file} | wc -l > {out_dir}/{prefix}.allCounts
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def PCA(vcf_file, prefix, temp_dir, out_dir):
    """ Template for identifying prune sites usning PLINK.
    I could make family IDs for plink somehow to account for relatedness, maybe later. For now just use only one of the related individuals."""
    inputs = {'file': vcf_file}
    outputs = {'Pruning_output': [f'{temp_dir}/{prefix}.prune.in',f'{temp_dir}/{prefix}.prune.out'],
               'PCA_output': [f'{out_dir}/{prefix}.eigenvec', f'{out_dir}/{prefix}.eigenval'],
               'Plink_bed': [f'{out_dir}/{prefix}.bed', f'{out_dir}/{prefix}.bim', f'{out_dir}/{prefix}.fam']}
    options={
        'cores': 1,
        'memory': '4g',
        'walltime': '2-00:00:00'
    }
    spec = f'''
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    plink --vcf {vcf_file} --double-id --allow-extra-chr \\
        --set-missing-var-ids @:# \\
        --indep-pairwise 50 10 0.1 --out {temp_dir}/{prefix}
    # prune and create pca
    plink --vcf {vcf_file} --double-id --allow-extra-chr --set-missing-var-ids @:# \\
        --extract {temp_dir}/{prefix}.prune.in \\
        --make-bed --pca --out {out_dir}/{prefix}
    echo "END: $(date)"
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def admixture(bed_file, bim_file, prefix, tmp_dir, out_dir):
    """ Template for running admixture on bed files. Fix output names!"""
    inputs = {'bed': bed_file,
              'bim': bim_file}
    outputs = {'CV_error': f'{tmp_dir}/{prefix}.cv.error'}
    options={
        'cores': 1,
        'memory': '16g',
        'walltime': '24:00:00'
    }
    spec = f'''
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
    awk '{{$1="0";print $0}}' {bim_file} > {bim_file}.tmp
    mv {bim_file}.tmp {bim_file}
    for i in {{2..5}}
    do
        admixture --cv {bed_file} $i > {tmp_dir}/log${{i}}.out
    done
    awk '/CV/ {{print $3,$4}}' {tmp_dir}/*out | cut -c 4,7-20 > {tmp_dir}/{prefix}.cv.error
    
    echo "END: $(date)"
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


