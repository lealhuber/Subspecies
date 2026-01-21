#!/bin/env python3
from gwf import AnonymousTarget # type: ignore
import os

# environment: popdiff

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


# improved version that uses start_codon info to decide promoter side
def get_regions_new(gtf, genesOI, add_bp, prefix, output_dir):
    """Make file with columns scaffold, start, end from gtf and genesOI"""
    inputs = [gtf, genesOI]
    outputs = {'regions_file': f'{output_dir}/{prefix}.bed',
               'regions_named_file': f'{output_dir}/{prefix}.with_names.bed',
               'no_start_codon': f'{output_dir}/{prefix}.no_start_codon.txt'
    }
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:10:00',
    }

    spec = f'''
    # plain GTF file (gene name is plain in column 9)
    # produce a 4-column file (scaffold,start,end,gene) where start/end are expanded
    # by add_bp on the promoter side only. We identify the promoter side by comparing
    # the start_codon position (only for transcript_id matching g*.t1) to the gene start.
    awk -F'\\t' 'NR==FNR {{ name=$1; gsub(/^[[:space:]]+|[[:space:]]+$/,"",name); gsub(/\r/,"",name); if (name != "") g[name]=1; next }}
     # record gene coordinates (try gene_id if present, else plain field)
     $3=="gene" {{
       if (match($9,/gene_id "([^"]+)"/,a)) name=a[1];
       else {{ name=$9; sub(/;.*/,"",name); gsub(/^[[:space:]]+|[[:space:]]+$/,"",name); gsub(/\r/,"",name) }}
       if (name in g) {{
         gene_start[name]=$4; gene_end[name]=$5; gene_scaf[name]=$1
         # if we already saw the start_codon for this gene, decide and print
         if (name in start_pos && !(name in printed)) {{
           sp = start_pos[name]; gs = gene_start[name]; ge = gene_end[name]
           if (sp <= gs) {{ s = gs - {add_bp}; if (s < 1) s = 1; e = ge }}
           else {{ s = gs; e = ge + {add_bp} }}
           print gene_scaf[name] "\\t" s "\\t" e "\\t" name "\\t" sp
           printed[name]=1
         }}
       }}
     }}
     # record start_codon position but only for transcript_id g*.t1 and extract gene_id
     $3=="start_codon" {{
       if (match($9,/transcript_id "([^"]+)"/,c) && match($9,/gene_id "([^"]+)"/,b)) {{
         tid=c[1]; gid=b[1]
         gsub(/\r/,"",tid); gsub(/\r/,"",gid)
         if (tid ~ /^g.*\\.t1$/) {{
           if (gid in g) {{
             start_pos[gid]=$4
             # debug: log matched start_codons to stderr so you can confirm matches
             print "START_CODON_MATCHED:" gid "\\t" start_pos[gid] > "/dev/stderr"
             # if gene coords already known, decide and print
             if (gid in gene_start && !(gid in printed)) {{
               gs = gene_start[gid]; ge = gene_end[gid]; sp = start_pos[gid]
               if (sp <= gs) {{ s = gs - {add_bp}; if (s < 1) s = 1; e = ge }}
               else {{ s = gs; e = ge + {add_bp} }}
               print gene_scaf[gid] "\\t" s "\\t" e "\\t" gid "\\t" sp
               printed[gid]=1
             }}
           }}
         }}
       }}
     }}
     END {{
       # fallback: if gene present but no start_codon info led to printing, use promoter on both sides
       for (nm in gene_start) {{
         if ((nm in g) && !(nm in printed)) {{
           s = gene_start[nm] - {add_bp}; if (s < 1) s = 1
           e = gene_end[nm] + {add_bp}
           print gene_scaf[nm] "\\t" s "\\t" e "\\t" nm "\\t" "NA"
         }}
       }}
       # debug: list genes that had no matching start_codon (one per line)
       for (nm in gene_start) {{
         if ((nm in g) && !(nm in start_pos)) {{
           print nm > "{output_dir}/{prefix}.no_start_codon.txt"
         }}
       }}
     }}' {genesOI} {gtf} > {output_dir}/{prefix}.with_names.bed

    # create the 3-column bed (scaffold,start,end) from the 5-column file
    cut -f1-3 {output_dir}/{prefix}.with_names.bed > {output_dir}/{prefix}.bed
    '''
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec) 

# definition for making file with just promoters from the genes of interest
def get_promoters(gtf, genesOI, promo_bp, prefix, output_dir):
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
    # add promo_bp on the promoter side only. We identify the promoter side by comparing
    # the start_codon position (only for transcript_id matching g*.t1) to the gene start.
    awk -F'\\t' 'NR==FNR {{ name=$1; gsub(/^[[:space:]]+|[[:space:]]+$/,"",name); gsub(/\r/,"",name); if (name != "") g[name]=1; next }}
     # record gene coordinates (try gene_id if present, else plain field)
     $3=="gene" {{
       if (match($9,/gene_id "([^"]+)"/,a)) name=a[1];
       else {{ name=$9; sub(/;.*/,"",name); gsub(/^[[:space:]]+|[[:space:]]+$/,"",name); gsub(/\r/,"",name) }}
       if (name in g) {{
         gene_start[name]=$4; gene_end[name]=$5; gene_scaf[name]=$1
         # if we already saw the start_codon for this gene, decide and print
         if (name in start_pos && !(name in printed)) {{
           sp = start_pos[name]; gs = gene_start[name]; ge = gene_end[name]
           if (sp <= gs) {{ s = gs - {promo_bp}; if (s < 1) s = 1; e = gs }}
           else {{ s = ge; e = ge + {promo_bp} }}
           print gene_scaf[name] "\\t" s "\\t" e "\\t" name "\\t" sp
           printed[name]=1
         }}
       }}
     }}
     # record start_codon position but only for transcript_id g*.t1 and extract gene_id
     $3=="start_codon" {{
       if (match($9,/transcript_id "([^"]+)"/,c) && match($9,/gene_id "([^"]+)"/,b)) {{
         tid=c[1]; gid=b[1]
         gsub(/\r/,"",tid); gsub(/\r/,"",gid)
         if (tid ~ /^g.*\\.t1$/) {{
           if (gid in g) {{
             start_pos[gid]=$4
             # debug: log matched start_codons to stderr so you can confirm matches
             print "START_CODON_MATCHED:" gid "\\t" start_pos[gid] > "/dev/stderr"
             # if gene coords already known, decide and print
             if (gid in gene_start && !(gid in printed)) {{
               gs = gene_start[gid]; ge = gene_end[gid]; sp = start_pos[gid]
               if (sp <= gs) {{ s = gs - {promo_bp}; if (s < 1) s = 1; e = gs }}
               else {{ s = ge; e = ge + {promo_bp} }}
               print gene_scaf[gid] "\\t" s "\\t" e "\\t" gid "\\t" sp
               printed[gid]=1
             }}
           }}
         }}
       }}
     }}
     END {{
       # debug: list genes that had no matching start_codon (one per line)
       for (nm in gene_start) {{
         if ((nm in g) && !(nm in start_pos)) {{
           print nm > "{output_dir}/{prefix}.no_start_codon.txt"
         }}
       }}
     }}' {genesOI} {gtf} > {output_dir}/{prefix}.with_names.bed

    # create the 3-column bed (scaffold,start,end) from the 5-column file
    cut -f1-3 {output_dir}/{prefix}.with_names.bed > {output_dir}/{prefix}.bed
    '''
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec) 


def get_coding_regions(gtf, prefix, output_dir):
    """Make file with columns scaffold, start, end for all coding regions from gtf"""
    inputs = [gtf]
    outputs = {'coding_regions_file': f'{output_dir}/{prefix}.coding_regions.bed'}
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:20:00',
    }

    spec = f'''
    # plain GTF file
    # produce a 3-column file (scaffold,start,end) for coding regions
    awk '$3=="CDS" {{ print $1 "\\t" $4 "\\t" $5 }}' {gtf} > {output_dir}/{prefix}.coding_regions.bed
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
      | sort -k1,1n -k3,3n \\
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
    # first make sure the snp_vcf has the same individuals as the pops_file:
    bcftools query -l {snp_vcf} | sort > {temp_dir}/snp_vcf_individuals.txt
    cut -f1 {pops_file} | sort > {temp_dir}/pops_file_ids.txt
    if ! diff {temp_dir}/snp_vcf_individuals.txt {temp_dir}/pops_file_ids.txt > /dev/null; then
      echo "Error: Individuals in SNP VCF and populations file do not match"
      diff {temp_dir}/snp_vcf_individuals.txt {temp_dir}/pops_file_ids.txt || true
      exit 1
    fi
    bedtools intersect -a {regions_file} -b {snp_vcf} -c > {output_dir}/{prefix}.snp_counts.bed
    echo "END: $(date)"
     '''
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def tajimas_d(vcf_file, regions_file, ind_file, prefix, temp_dir, out_dir):
    inputs = [vcf_file, regions_file]
    outputs = {'tajimas_d_file': f'{out_dir}/{prefix}.Tajima.D'}
    options = {
        'cores': 4,
        'memory': '24g',
        'walltime': '05:00:00',
    }

    spec = f'''
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    # vcftools expects a header in bed file so add one to regions_file
    echo -e "CHROM\\tSTART\\tEND" | cat - {regions_file} > {temp_dir}/{prefix}.bed_with_header

    # calculate Tajima's D (made bin very big to capture all variants in each region in one statistic)
    # only for individuals in ind_file
    vcftools --gzvcf {vcf_file} --keep {ind_file} --bed {temp_dir}/{prefix}.bed_with_header --TajimaD 50000 --out {out_dir}/{prefix}

    echo "END: $(date)"
     '''
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def pixy_stats(vcf_file, regions_file, pops_file, prefix, out_dir):
    inputs = [vcf_file, regions_file]
    outputs = {'pixy_dxy': f'{out_dir}/{prefix}.pixy_dxy.txt',
               'pixy_fst': f'{out_dir}/{prefix}.pixy_fst.txt',
               'pixy_pi': f'{out_dir}/{prefix}.pixy_pi.txt',
               'pixy_tajima_d': f'{out_dir}/{prefix}.pixy_tajima_d.txt'
    }
    options = {
        'cores': 4,
        'memory': '24g',
        'walltime': '10:00:00',
    }

    spec = f'''
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    # calculate pixy stats
    pixy --stats pi dxy fst tajima_d --vcf {vcf_file} --populations {pops_file} --bed_file {regions_file} --n_cores {options['cores']} \\
      --output_folder {out_dir} --output_prefix {prefix}.pixy --fst_type hudson

    echo "END: $(date)"
     '''
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)