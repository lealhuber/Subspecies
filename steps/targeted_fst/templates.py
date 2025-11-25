#!/bin/env python3
from gwf import AnonymousTarget # type: ignore
import os


# definition for making file with regions of interest from gtf and list of gene names
def get_regions(gtf, genesOI, output_dir):
    """Make file with columns scaffold, start, end from gtf and genesOI"""
    inputs = [gtf, genesOI]
    outputs = {'regions_file': f'{output_dir}/regions_of_interest.txt'}
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:10:00',
    }

    spec = f'''
    # plain GTF file (gene name is plain in column 9)
    awk 'NR==FNR {{ g[$1]=1; next }}
     $3=="gene" {{
       name=$9
       # remove anything after a semicolon (in case attributes present) and trim whitespace
       sub(/;.*/,"",name)
       gsub(/^[[:space:]]+|[[:space:]]+$/,"",name)
       if (name in g) print $1 "\t" $4 "\t" $5
     }}' {genesOI} {gtf} > {output_dir}/regions_of_interest.txt
    '''
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec) 

def pop_stats(vcf_file, script_path, regions_file, pops_file, prefix, temp_dir, output_dir):
    """Calculate population statistics for given regions"""
    inputs = [vcf_file, regions_file]
    outputs = {'stats_file': f'{output_dir}/regionsOI_fst.csv.gz'}
    options = {
        'cores': 4,
        'memory': '8g',
        'walltime': '05:00:00',
    }

    spec = f'''
    # make geno file if it doesn't exist
    if [ ! -f {temp_dir}/{prefix}.input.geno.gz ]; then
      python {script_path}/VCF_processing/parseVCFs.py -i {vcf_file} \\
        --threads {options['cores']} -o {temp_dir}/{prefix}.input.geno.gz
    fi
    # calculate population statistics
    python {script_path}/popgenWindows.py --windType predefined --windCoords {regions_file} \\
      --popsFile {pops_file} -p black -p blue -p red \\
      -g {temp_dir}/{prefix}.input.geno.gz -o {output_dir}/{prefix}.regionsOI_fst.csv.gz \\
      -f phased -T {options['cores']}
     '''
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

