#!/bin/env python3
from gwf import AnonymousTarget # type: ignore
import os

# With all these filters, I shouldn't lose invariant sites, so make sure of that when filtering.

def PCA(vcf_file, prefix, temp_dir, out_dir):
    """ Template for identifying prune sites usning PLINK.
    I could make family IDs for plink somehow to account for relatedness, maybe later. For now just use only one of the related individuals."""
    inputs = {'file': vcf_file}
    outputs = {'Pruning_output': [f'{temp_dir}/{prefix}.prune.in',f'{temp_dir}/{prefix}.prune.out'],
               'PCA_output': [f'{out_dir}/{prefix}.eigenvec', f'{out_dir}/{prefix}.eigenval'],
               'Plink_bed': [f'{out_dir}/{prefix}.bed', f'{out_dir}/{prefix}.bim', f'{out_dir}/{prefix}.fam']}
    protect = [outputs['PCA_output'], outputs['Plink_bed']]
    options={
        'cores': 1,
        'memory': '4g',
        'walltime': '24:00:00'
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
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def admixture(bed_file, bim_file, prefix, tmp_dir, out_dir):
    """ Template for running admixture on bed files. Fix output names!"""
    inputs = {'bed': bed_file,
              'bim': bim_file}
    outputs = {'Q_output': f'{out_dir}/{prefix}.Q',
               'P_output': f'{out_dir}/{prefix}.P',
               'log': f'{out_dir}/{prefix}.log'}
    options={
        'cores': 4,
        'memory': '16g',
        'walltime': '24:00:00'
    }
    spec = f'''
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    # ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
    awk '{{${{1}}="0";print $0}}' {bim_file} > {bim_file}.tmp
    mv {bim_file}.tmp {bim_file}
    for i in {{2..5}}
    do
        admixture --cv {bed_file} $i > {tmp_dir}/log${{i}}.out
    done
    awk '/CV/ {{print $3,$4}}' {tmp_dir}*out | cut -c 4,7-20 > {tmp_dir}/{prefix}.cv.error
    echo "END: $(date)"
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


