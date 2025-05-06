#!/bin/env python3
from gwf import AnonymousTarget # type: ignore
import os

def add_RG(input_bam, sample_name):
    """add read groups because I forgot them like an ass. Also requires reindexing"""
    inputs={'bam_files': input_bam}
    outputs={'RG_bams': input_bam,
             'RG_bais': f'{input_bam}.bai'}
    options={
    'cores': 4,
    'memory': '32g',
    'walltime': '04:00:00'
    }
    spec=f'''
    samtools addreplacerg \\
        --threads {options['cores']} \\
        -r ID:{sample_name} \\
        -r SM:{sample_name} \\
        {input_bam} \\
        | samtools sort \\
            --threads {options['cores']} \\
            -o {input_bam} -

        samtools index {input_bam}
    '''
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
def parse_fasta(fastaFile: str):
	"""
	Parses :format:`FASTA` file of reference genome returning all sequence names and lengths paired in a list of dictionaries.
	
	::
	
		return [{'sequence_name': str, 'sequence_length': int}, ...]
	
	:param str fastaFile:
		Sequence file in :format:`FASTA` format.
	"""
	fastaList = []
	seqName = None
	length = 0
	with open(fastaFile, 'r') as fasta:
		for entry in fasta:
			entry = entry.strip()
			if entry.startswith(">"):
				if seqName:
					fastaList.append({'sequence_name': seqName, 'sequence_length': length})
					length = 0
				entry = entry.split(" ", 1)
				seqName = entry[0][1:]
			else:
				length += len(entry)
		fastaList.append({'sequence_name': seqName, 'sequence_length': length})
	return fastaList

def padding_calculator(parseFasta: list, size: int  | None = 500000):
	"""
	Calculates proper 0 padding for numbers in **partition_chrom**.

	:param list parseFasta:
		List of dictionaries produced by **parse_fasta**.
	:param int size:
		Size of partitions. Default 500kb.
	"""
	num = 1
	for chrom in parseFasta:
		wholeChunks = chrom['sequence_length'] // size
		num += (wholeChunks + 1)
	return len(str(num))

def partition_chrom(parseFasta: list, size: int = 500000, nPad: int = 5):
	"""
	Partitions :format:`FASTA` file parsed with **parse_fasta**.
	
	Uses the list of dictionaries from **parse_fasta** to creates a list of dictionaries
	containing partition number, sequence name, start and end postion (0 based).
	By default the partition size is 500kbs.
	
	::
	
		return [{'num': int, 'region': str, 'start': int, 'end': int}, ...]
	
	:param list parseFasta:
		List of dictionaries produced by **parse_fasta**.
	:param int size:
		Size of partitions. Default 500kb.
	"""
	chromPartition = []
	num = 1
	for chrom in parseFasta:
		wholeChunks = chrom['sequence_length'] // size
		partialChunk = chrom['sequence_length'] - wholeChunks * size
		start = 0
		for chunk in range(wholeChunks):
			end = start + size
			chromPartition.append({'num': f'{num:0{nPad}}', 'region': chrom['sequence_name'], 'start': start, 'end': end})
			start = end
			num += 1
		if partialChunk:
			chromPartition.append({'num': f'{num:0{nPad}}', 'region': chrom['sequence_name'], 'start': start, 'end': start + partialChunk})
			num += 1
	return chromPartition


#CALL chr vcfs 
def freebayes_CHR_vcf(files: list, reference_genome: str, population_match, temp_path: str, output_name: str, region, start, end):
    """Calling all specified indiviudals to CHR vcfs"""
    inputs={'vcfs': files}
    outputs={'genome_chr_vcf': temp_path + output_name}
    options={
    'cores': 1,
    'memory': '144g',
    'walltime': '2-00:00:00'
    }
    spec="""
    freebayes \
        -f {reference_genome} \
        --report-monomorphic \
        -b {bam_files} \
		--populations {populations} \
        -r {region}:{start}-{end} \
        > {output}
    """.format(bam_files=' '.join(files), reference_genome=reference_genome, populations=population_match, region=region, start=start, end=end, output = temp_path + output_name)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

#concat CHR vcfs -> genome vcf
def concat_CHR_vcf(files: list, species_name: str, path: str):
    """Concat CHR vcfs -> genome vcf"""
    output_file = '{path}{name}.vcf'.format(path=path, name=species_name)
    inputs={'vcfs': files}
    outputs={'genome_vcf': output_file}
    options={
    'cores': 1,
    'memory': '8g'
    }
    spec="""
    bcftools concat \
        {VCFs} \
        > {output}
    """.format(VCFs=' '.join(files), output=output_file)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

#Sort genome vcf -> sorted.vcf
def Sort_vcf(vcf: str, species_name: str, path: str, samples: list):
    """Sort vcf"""
    output_file = '{path}sorted_{name}.vcf'.format(path=path, name=species_name)
    inputs={'genome_vcf': vcf}
    outputs={'sorted_genome_vcf': output_file}
    options={
    'cores': 1,
    'memory': '8g',
    'walltime': '2-00:00:00'
    }
    spec="""
    bcftools view -s {samples} {vcf} | bcftools sort > {output}
    """.format(vcf=vcf, output=output_file, samples=','.join(samples))
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# Zip VCF 
def zip_file(file: str, species_name: str):
    """zip genome file -> .vcf.gz """
    output_file = '{path}/sorted_{name}.vcf.gz'.format(path=os.path.dirname(file), name=species_name)
    inputs={'genome_vcf': file}
    outputs={'vcf_gz': output_file}
    protect=outputs['vcf_gz']
    options={
    'cores': 8,
    'memory': '64g',
    'walltime': '2-00:00:00'
    }
    spec="""
    bcftools view -Oz -o {output} {VCF}
    bcftools index {output}
    bcftools stats {output} > {name}.stats
    """.format(VCF=file, output=output_file, name=species_name)
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

# Stats on genome_vcf, individual stats 
def vcf_stats(file: str, sample: str, stats_path: str, this_stat:str):
    """Statistics on genome_vcf.gz, individual stats"""
    output_file = '{path}/{sample}.stats'.format( \
        path=stats_path, sample=sample)
    inputs={'vcf_gz': file}
    outputs={'vcf_pop_stats': this_stat}
    protect=outputs['vcf_pop_stats']
    options={
    'cores': 1,
    'memory': '16g',
    'walltime': '10:00:00'
    }
    spec="""
    #in order to summarrize the individuals stats with multiqc - the vcf needs to be "different" 

    #stats pr. sample
    bcftools view -s {sample} {VCF} > {path}{sample}.vcf
    bcftools view -Oz -o {path}{sample}.vcf.gz {VCF}
    bcftools index -c {path}{sample}.vcf.gz
    bcftools stats -s {sample} {path}{sample}.vcf.gz > {out}

    rm {path}{sample}.vcf
    rm {path}{sample}.vcf.gz
    """.format(VCF=file, path=stats_path, sample=sample, out=output_file)
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)

# Multiqc genome_stats, individual stats 
def multiqc_vcf_stats(files: str, species_name: str, result_path: str):
    """Multiqc vcf_stats, individual stats -> .html """
    output_file = '{path}/{name}.html'.format(path=result_path, name=species_name)
    inputs={'vcf_pop_stats': files}
    outputs={'multiqc': output_file}
    protect=outputs['multiqc']
    options={
    'cores': 1,
    'memory': '8g'
    }
    spec="""
    cd {path} || exit
    multiqc . 
    mv multiqc_report.html {output}
    """.format(path=os.path.dirname(files[-1]), output=output_file)
    return AnonymousTarget(inputs=inputs, outputs=outputs, protect=protect, options=options, spec=spec)