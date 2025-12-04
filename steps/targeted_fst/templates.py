#!/bin/env python3
from gwf import AnonymousTarget # type: ignore
import os


# definition for making file with regions of interest from gtf and list of gene names
def get_regions(gtf, genesOI, add_bp, prefix, output_dir):
    """Make file with columns scaffold, start, end from gtf and genesOI"""
    inputs = [gtf, genesOI]
    outputs = {'regions_file': f'{output_dir}/{prefix}.bed',
               'regions_named_file': f'{output_dir}/{prefix}.with_names.bed'
    }
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:10:00',
    }

    spec = f'''
    # plain GTF file (gene name is plain in column 9)
    # produce a 4-column file (scaffold,start,end,gene) then cut out the 4th column for the plain bed
    awk 'NR==FNR {{ g[$1]=1; next }}
     $3=="gene" {{
       name=$9
       # remove anything after a semicolon (in case attributes present) and trim whitespace
       sub(/;.*/,"",name)
       gsub(/^[[:space:]]+|[[:space:]]+$/,"",name)
       if (name in g) {{
         s = $4 - {add_bp}
         if (s < 1) s = 1
         e = $5 + {add_bp}
         print $1 "\\t" s "\\t" e "\\t" name
       }}
     }}' {genesOI} {gtf} > {output_dir}/{prefix}.with_names.bed

    # create the 3-column bed (scaffold,start,end) from the 4-column file
    cut -f1-3 {output_dir}/{prefix}.with_names.bed > {output_dir}/{prefix}.bed
    '''
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec) 

def pop_stats(vcf_file, snp_vcf, script_path, regions_file, pops_file, prefix, geno_prefix, temp_dir, output_dir):
    """Calculate population statistics for given regions"""
    inputs = [vcf_file, regions_file]
    outputs = {'popstats_file': f'{output_dir}/{prefix}.regionsOI_fst.csv',
               'sncounts_file': f'{output_dir}/{prefix}.snp_counts.bed'}
    options = {
        'cores': 4,
        'memory': '24g',
        'walltime': '05:00:00', # needs more time when making geno file
    }

    spec = f'''
    echo "START: $(date)"
	  echo "JobID: $SLURM_JOBID"
    # make geno file if it doesn't exist
    if [ ! -f {temp_dir}/{geno_prefix}.input.geno.gz ]; then
      python {script_path}/VCF_processing/parseVCFs.py -i {vcf_file} \\
        --threads {options['cores']} -o {temp_dir}/{geno_prefix}.input.geno.gz
      echo "Geno file created"
    else
      echo "Geno file already exists, skipping creation"
    fi

    # get the order of scaffolds from the geno file
    if [ ! -f {temp_dir}/scaffold_order.txt ]; then
      zcat {temp_dir}/{geno_prefix}.input.geno.gz | awk '!/^#/ && !seen[$1]++ {{ print $1 }}' > {temp_dir}/scaffold_order.txt
    fi
    # filter out the ones that are not in the geno file (because on too small scaffold)
    awk 'NR==FNR {{keep[$1]=1; next}} ($1 in keep)' {temp_dir}/scaffold_order.txt {regions_file} > {regions_file}.filtered.bed
    # sort the bed file by the order of scaffolds from the geno file
    awk 'NR==FNR {{order[$1]=NR; next}} {{print order[$1] "\\t" $0}}' {temp_dir}/scaffold_order.txt {regions_file}.filtered.bed \\
      | sort -k1,1n -k2,2n \\
      | cut -f2- \\
      > {regions_file}.sorted.bed
    mv {regions_file}.sorted.bed {regions_file}
    echo "finished sorting: $(date)"

    # calculate population statistics
    python {script_path}/popgenWindows.py --windType predefined --windCoords {regions_file} \\
      --popsFile {pops_file} -p black -p blue -p red \\
      -g {temp_dir}/{geno_prefix}.input.geno.gz -o {output_dir}/{prefix}.regionsOI_fst.csv \\
      -f phased -T {options['cores']} --verbose --writeFailedWindows
    echo "Finished population statistics"

    # it will also be nice to know the number of variants per region, do that too, using the SNP vcf file I already have
    bedtools intersect -a {regions_file} -b {snp_vcf} -c > {output_dir}/{prefix}.snp_counts.bed
    echo "END: $(date)"
     '''
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

