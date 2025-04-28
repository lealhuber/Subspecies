#!/bin/env python3
from gwf import AnonymousTarget # type: ignore
import os, glob
from gwf import Workflow # type: ignore
import numpy as np # type: ignore

gwf = Workflow(defaults={'account': 'ostrich_thermal'})

def FastQC(input_fasta_pair1, input_fasta_pair2, sample_name, out_dir):
    """Template for running FastQC on reads.
    Whether the val is needed in the QC files is just a quick, ugly fix, sorry. Make nice later!!"""
    inputs = {'read1': input_fasta_pair1,
              'read2': input_fasta_pair2}
    outputs = {'QC_files': [f'{out_dir}/{sample_name}_val_1_fastqc.html',
                            f'{out_dir}/{sample_name}_val_1_fastqc.zip',
                            f'{out_dir}/{sample_name}_val_2_fastqc.html',
                            f'{out_dir}/{sample_name}_val_2_fastqc.zip']}
    options = {
        'cores': 4,
        'memory': '16g',
        'walltime': '09:00:00'
    }

    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
    fastqc -t 4 -o {out_dir} {input_fasta_pair1} {input_fasta_pair2}
	'''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def adapterremoval(input_fasta_pair1, input_fasta_pair2, sample_name, out_dir, min_length):
    """Template for paired end adapter removal.
    Trim galore guesses adapters automatically and trims ends with Phred scores below 20."""

    inputs = {'read1': input_fasta_pair1,
              'read2': input_fasta_pair2}
    outputs = {'pair1': f'{out_dir}/{sample_name}_val_1.fq.gz',
               'pair2': f'{out_dir}/{sample_name}_val_2.fq.gz',
               'reports': [f'{out_dir}/{sample_name}_1.fastq.gz_trimming_report.txt',
                      f'{out_dir}/{sample_name}_2.fastq.gz_trimming_report.txt']}
    options = {
        'cores': 8,
		'memory': '32g',
		'walltime': '12:00:00'
    }

    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

    trim_galore --cores {options['cores']} --paired --gzip \\
        --length {min_length} --trim-n \\
        --o {out_dir} --basename {sample_name} \\
       {input_fasta_pair1} {input_fasta_pair2}
	'''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def bwa_index(ref_genome, out_dir):
    """Template for indexing a genome with `bwa index`."""
    inputs = {'reference': ref_genome}
    outputs = {'bwa': [f'{out_dir}/{os.path.basename(ref_genome)}.amb',
					   f'{out_dir}/{os.path.basename(ref_genome)}.ann',
					   f'{out_dir}/{os.path.basename(ref_genome)}.pac',
					   f'{out_dir}/{os.path.basename(ref_genome)}.bwt',
					   f'{out_dir}/{os.path.basename(ref_genome)}.sa']}
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '09:00:00'
    }

    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
    bwa index -p {out_dir}/{os.path.basename(ref_genome)} -a bwtsw {ref_genome}
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mapping_pairdend(input_fasta_pair1,input_fasta_pair2,sample_name,ref_genome,out_dir):
    """Template of mapping a pair of fasta files to a reference genome using bwa mem"""
    inputs = {'read1': input_fasta_pair1,
              'read2': input_fasta_pair2,
              'reference': ref_genome}
    outputs = {'alignment': f'{out_dir}/{sample_name}.mapped.bam'}
    options = {
        'cores': 16,
        'memory': '32g',
        'walltime': '09:00:00'
    }
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

    bwa mem \\
        -t {options['cores']} \\
        -R "@RG\\tID:{sample_name}\\tSM:{sample_name}" \\
        {ref_genome} {input_fasta_pair1} {input_fasta_pair2} \\
    | samtools sort \\
        --threads {options['cores']} \\
        -n -o {out_dir}/{sample_name}.mapped.bam -
        
    echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"    
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def mark_dups_samtools(alignment_file, sample_name, out_dir):
    """Template: Mark duplicate alignments using samtools markdup.
    There, the -m optiton adds a mate score tag, """
    inputs = {'alignment': alignment_file}
    outputs = {'markdup': f'{out_dir}/{sample_name}.markdup.bam',
               'bai': f'{out_dir}/{sample_name}.markdup.bam.bai',
               'stats': f'{out_dir}/{sample_name}.markdup.bam.stats'}
    options = {
		'cores': 18,
		'memory': '60g',
		'walltime': '24:00:00'
	}
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

    samtools fixmate \\
        --threads {options['cores'] - 1} \\
        		-m \\
		--output-fmt BAM \\
		{alignment_file} \\
		- \\
	| samtools sort \\
		--threads {options['cores'] - 1} \\
		--output-fmt BAM \\
		- \\
	| samtools markdup \\
		--threads {options['cores'] - 1} \\
		--output-fmt BAM \\
		-s \\
		-f {out_dir}/{sample_name}.markdup.bam.stats \\
		- \\
		{out_dir}/{sample_name}.markdup.bam
	
	samtools index \\
		--threads {options['cores'] - 1} \\
		-b \\
		{out_dir}/{sample_name}.markdup.bam \\
		{out_dir}/{sample_name}.markdup.bam.bai

    echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def samtools_stats(alignment_file, sample_name, out_dir):
    """ Template to create various mapping statistics using samtools"""
    inputs = {'alignment': alignment_file}
    outputs = {'stats': [f'{out_dir}/{sample_name}.idxstats',
                         f'{out_dir}/{sample_name}.flagstat',
                         f'{out_dir}/{sample_name}.coverage',
                         f'{out_dir}/{sample_name}.stats']}
    options = {
		'cores': 8,
		'memory': '32g',
		'walltime': '12:00:00'
	}
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

    samtools idxstats \\
		--threads {options['cores']} \\
		{alignment_file} \\
		> {out_dir}/{sample_name}.idxstats
	
	samtools flagstat \\
		--threads {options['cores']} \\
		{alignment_file} \\
		> {out_dir}/{sample_name}.flagstat

	samtools coverage \\
		-o {out_dir}/{sample_name}.coverage \\
		{alignment_file}
	
	samtools stats \\
		--threads {options['cores']} \\
		--coverage 1,1000,1 \\
		{alignment_file} \\
		> {out_dir}/{sample_name}.stats

    echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def samtools_filter(alignment_file, sample_name, out_dir, min_mq):
    """ Template for filtering out low quality reads and duplicates using samtools"""
    inputs = {'alignment': alignment_file}
    outputs = {'filtered': f'{out_dir}/{sample_name}.dedup.filtered.bam',
               'bai': f'{out_dir}/{sample_name}.dedup.filtered.bam.bai'}
    protect = [outputs['filtered'], outputs['bai']]
    options = {
        'cores': 16,
		'memory': '64g',
		'walltime': '24:00:00'
    }
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

    samtools view \\
		--threads {options['cores'] - 1} \\
		--bam \\
        --exclude-flags 0x400 \\
		--min-MQ {min_mq} \\
		--output {out_dir}/{sample_name}.dedup.filtered.bam \\
		{alignment_file}

	samtools index \\
		--threads {options['cores'] - 1} \\
		--bai \\
		--output {out_dir}/{sample_name}.dedup.filtered.bam.bai \\
		{out_dir}/{sample_name}.dedup.filtered.bam
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
        
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

def qc_qualimap(alignment_file, sample_name, out_dir):
    """ Template for running qualimap qc of alignments because samtools stats is not good enough apparently.
    Add mkdir {out_dir}/{sample_name} in there if you do it again from scrach """

    inputs = {'alignment': alignment_file}
    outputs = {'html': f'{out_dir}/{sample_name}/qualimapReport.html',
			   'raw': f'{out_dir}/{sample_name}/genome_results.txt'}
    protect = [outputs['html'], outputs['raw']]
    options = {
		'cores': 16,
		'memory': '64g',
		'walltime': '24:00:00'
	}
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

    export _JAVA_OPTIONS="-Djava.awt.headless=true -Xmx{options['memory']}"
	
	qualimap bamqc \\
		-nt {options['cores']} \\
		-bam {alignment_file} \\
		-outdir {out_dir}/{sample_name} \\
		-outformat PDF:HTML \\
		-outfile {out_dir}/{sample_name}/qualimapReport.pdf \\
		--java-mem-size={options['memory']}

    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)


# very easy test workflow: save head of each file pair concatenated together

def save_head(input_fasta_pair1,input_fasta_pair2,out_dir):
    """Template for stupidly saving the output of head of both reads of a pair"""
    inputs = {'read1': input_fasta_pair1, 'read2': input_fasta_pair2}
    outputs = [f'{out_dir}/head_both']
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '00:30:00',
    }
    spec = f'''
    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
    echo "Input1: ${input_fasta_pair1}"
    echo "Output: ${out_dir}"
    (head -n 1 {input_fasta_pair1} ; head -n 1 {input_fasta_pair2}) > {out_dir}/head_both
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
