#!/bin/env python3
from gwf import AnonymousTarget # type: ignore
import os

#CALL chr vcfs 
def freebayes_CHR_vcf(files: list, reference_genome: str, temp_path: str, output_name: str, chromosome: str):
    """Calling all specified indiviudals to CHR vcfs"""
    inputs={'vcfs': files}
    outputs={'genome_chr_vcf': temp_path + output_name}
    options={
    'cores': 1,
    'memory': '100g',
    'walltime': '3-00:00:00'
    }
    spec="""
    freebayes \
        -f {reference_genome} \
        --report-monomorphic \
        -b {bam_files} \
        -r {chr} \
        > {output}
    """.format(bam_files=' '.join(files), reference_genome=reference_genome, chr=chromosome, output=temp_path + output_name)
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
    options={
    'cores': 8,
    'memory': '64g',
    'walltime': '2-00:00:00'
    }
    spec="""
    bcftools view -Oz -o {output} {VCF}
    bcftools index {output}
    bcftools stats {output} > {name}.states
    """.format(VCF=file, output=output_file, name=species_name)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# States on genome_vcf, individual stats 
def vcf_stats(file: str, sample: str, stats_path: str, this_stat:str):
    """Statistics on genome_vcf.gz, individual stats"""
    output_file = '{path}/{sample}.stats'.format( \
        path=stats_path, sample=sample)
    inputs={'vcf_gz': file}
    outputs={'vcf_pop_stats': this_stat}
    options={
    'cores': 1,
    'memory': '16g',
    'walltime': '10:00:00'
    }
    spec="""
    #in order to summarrize the individuals states with multiqc - the vcf needs to be "different" 

    #stats pr. sample
    bcftools view -s {sample} {VCF} > {path}{sample}.vcf
    bcftools view -Oz -o {path}{sample}.vcf.gz {VCF}
    bcftools index -c {path}{sample}.vcf.gz
    bcftools stats -s {sample} {path}{sample}.vcf.gz > {out}

    rm {path}{sample}.vcf
    rm {path}{sample}.vcf.gz
    """.format(VCF=file, path=stats_path, sample=sample, out=output_file)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# Multiqc genome_stats, individual stats 
def multiqc_vcf_stats(files: str, species_name: str, result_path: str):
    """Multiqc vcf_stats, individual stats -> .html """
    output_file = '{path}/{name}.html'.format(path=result_path, name=species_name)
    inputs={'vcf_pop_stats': files}
    outputs={'multiqc': output_file}
    options={
    'cores': 1,
    'memory': '8g'
    }
    spec="""
    cd {path} || exit
    multiqc . 
    mv multiqc_report.html {output}
    """.format(path=os.path.dirname(files[-1]), output=output_file)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)