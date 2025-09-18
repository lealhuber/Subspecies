#!/bin/env python3
from gwf import AnonymousTarget # type: ignore
import os

# With all these filters, I shouldn't lose invariant sites, so make sure of that when filtering.

def random_subset(vcf_file, prefix, out_dir):
    """ Template for getting a random subset from a vcf file for faster stats testing"""
    inputs = {'file': vcf_file}
    outputs = {'subset_file': f'{out_dir}/{prefix}.subset.vcf.gz'}
    options={
        'cores': 1,
        'memory': '4g',
        'walltime': '24:00:00'
    }
    spec = f'''
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    bcftools view {vcf_file} | vcfrandomsample -r 0.01 > {out_dir}/{prefix}.subset.vcf
    # compress vcf
    bgzip {out_dir}/{prefix}.subset.vcf
    # index vcf
    bcftools index {out_dir}/{prefix}.subset.vcf.gz
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
        'walltime': '3-00:00:00'
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
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf --freq2 --max-alleles 2 --out {out_dir}/{prefix}
    # mean depth of coverage per individual
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf --depth --out {out_dir}/{prefix}
    # mean depth per site
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf --site-mean-depth --out {out_dir}/{prefix}
    # site quality
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf --site-quality --out {out_dir}/{prefix}
    # proportion of missing data per individual
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf --missing-indv --out {out_dir}/{prefix}
    # proportion of missing data per site
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf --missing-site --out {out_dir}/{prefix}
    # calculate heterozygosity and inbreeding coefficient per individual
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf --het --out {out_dir}/{prefix}
    # caculate relatedness
    vcftools --gzvcf {out_dir}/{prefix}.subset.vcf --relatedness --out {out_dir}/{prefix}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def vcf_stats(vcf_file, prefix, out_dir):
    """ Template for vcf stats. Adapted from https://speciationgenomics.github.io/filtering_vcfs/"""
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
        'walltime': '01:00:00'
    }
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
    # Count the number of  variants (takes for ever!)
    # bcftools view -H {vcf_file} | wc -l > {out_dir}/{prefix}.allCounts
    # calculate allele frequency (for sites with max 2 alleles)
    vcftools --gzvcf {vcf_file} --freq2 --max-alleles 2 --out {out_dir}/{prefix}
    # mean depth of coverage per individual
    vcftools --gzvcf {vcf_file} --depth --out {out_dir}/{prefix}
    # mean depth per site
    vcftools --gzvcf {vcf_file} --site-mean-depth --out {out_dir}/{prefix}
    # site quality
    vcftools --gzvcf {vcf_file} --site-quality --out {out_dir}/{prefix}
    # proportion of missing data per individual
    vcftools --gzvcf {vcf_file} --missing-indv --out {out_dir}/{prefix}
    # proportion of missing data per site
    vcftools --gzvcf {vcf_file} --missing-site --out {out_dir}/{prefix}
    # calculate heterozygosity and inbreeding coefficient per individual
    vcftools --gzvcf {vcf_file} --het --out {out_dir}/{prefix}
    # caculate relatedness
    vcftools --gzvcf {vcf_file} --relatedness --out {out_dir}/{prefix}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def allele_filter(vcf_file, prefix, out_dir):
    """ Template for filtering out multiallelic sites and removing indels,
    while keeping monomorphic sites for dxy and pi (no maf filter!).
    For Fst one does not need monomorphic sites but they don't hurt I think/hope."""
    inputs = {'file': vcf_file}
    outputs = {'filtered_file': f'{prefix}.monobiallelic.snp.vcf.gz'}
    options={
        'cores': 1,
        'memory': '32g',
        'walltime': '2-00:00:00'
    }
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
    vcftools --gzvcf {vcf_file} \\
    --remove-indels --max-alleles 2 \\
    --recode --recode-INFO-all \\
    --stdout | bgzip -c > {prefix}.monobiallelic.snp.vcf.gz
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def quality_filter(vcf_file, prefix, out_dir, max_missing, min_DP, max_DP):
    """ Template for filtering out sites passing basic best standart quality filters"""
    inputs = {'file': vcf_file}
    outputs = {'filtered_file': f'{out_dir}/{prefix}.filtered.vcf.gz'}
    options={
        'cores': 1,
        'memory': '32g',
        'walltime': '36:00:00'
    }
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
    vcftools --gzvcf {vcf_file} \\
    --max-missing {max_missing} \\
    --min-meanDP {min_DP} --max-meanDP {max_DP} \\
    # because one sample is less deep but overall fine, don't filter genotype depth
    # --minDP {min_DP} --maxDP {max_DP} \\
    --recode \\
    --stdout | bgzip -c > {out_dir}/{prefix}.filtered.vcf.gz
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def bias_filter(vcf_file, prefix, out_dir, avgDP, qual_div_avgDP):
    """ Template for filtering out sites with mapping quality bias, position bias and strand bias using bcftools
    (suggestions by htslib.org/workflow/filter.html). Only problem my vcfs don't contain any of those info fields..."""
    inputs = {'file': vcf_file}
    outputs = {'filtered_file': f'{out_dir}/{prefix}.nobias.vcf.gz'}
    options={
        'cores': 1,
        'memory': '32g',
        'walltime': '2-00:00:00'
    }
    spec = f'''
    bcftools view -e "MQBZ < -(3.5+4*{qual_div_avgDP}) || \\
        RPBZ > (3+3*{qual_div_avgDP}) || RPBZ < -(3+3*{qual_div_avgDP}) || FORMAT/SP > (40+{avgDP}/2)" \\
        -Oz -o {out_dir}/{prefix}.nobias.vcf.gz {vcf_file}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def HWE_filter(vcf_file, prefix, samples_list, temp_dir, out_dir, stat_dir):
    """ Template for getting samples of one population and removing sites that deviate from HWE within population."""
    inputs = {'file': vcf_file}
    outputs = {'split_file': f'{temp_dir}/{prefix}.vcf.gz',
               'filtered_file': f'{out_dir}/{prefix}.HWE.vcf.gz'}
    options={
        'cores': 1,
        'memory': '32g',
        'walltime': '10:00:00'
    }
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
    bfcftools view -s {samples_list} -Oz -o {temp_dir}/{prefix}.vcf.gz {vcf_file}

    vcftools --gzvcf {temp_dir}/{prefix}.vcf.gz --hwe 0.05 \\
    --recode --recode-INFO-all --stdout | bgzip -c > ${out_dir}/{prefix}.HWE.vcf.gz

    bcftools view -H {out_dir}/{prefix}.HWE.vcf.gz | wc -l \\
    > {stat_dir}/{prefix}.hwe.allCounts
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
