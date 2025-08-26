## declare the workflow
from gwf import Workflow # type: ignore
import numpy as np # type: ignore
import re
from templates import *

gwf = Workflow(defaults={'account': 'ostrich_thermal'})

import os # module provides functions for interacting with operating system, eg with files
import glob # module for regex patterns for finding file patterns or names

#############################
### index ref genome
#############################

out_dir = '/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/bwa_indexed'

index = gwf.target_from_template(
    name='IndexGenome',
    template=bwa_index(
        ref_genome='/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/bwa_indexed/Struthio_camelus_HiC.fasta',
        out_dir='/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/bwa_indexed'
    )
)

base_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/data_re'
subspecies = ["black", "blue", "red"]

for SUBSPEC in subspecies:
    # set directories for output and for QC
    outDir = f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/mapping/outputs/{SUBSPEC}'
    outQC = f'/faststorage/project/ostrich_thermal/BACKUP/subpop_dna/data/individuals_fq/{SUBSPEC}/FastQC'
    outQC2 = f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/mapping/QC/{SUBSPEC}'
    # make list of all fasta files
    file_list = glob.glob(f"{base_dir}/{SUBSPEC}/*.fastq.gz", recursive=True)
    sample_names = set() # Use a set to store unique sample names
    for file in file_list:
        filename = os.path.basename(file)
        sample_name = re.sub(r"_[12]\.fastq\.gz$", "", filename)  # Remove _1.fastq.gz or _2.fastq.gz
        sample_names.add(sample_name) # Add to set (duplicates are ignored)
    sample_names = sorted(sample_names) # sort it for niceness
    # print(sample_names)
    sampleIDs = set() # need later for merging samples run on different platform units

    # start of actual workflow going through each sample
    for sample in sample_names:

        read1 = f'{base_dir}/{SUBSPEC}/{sample}_1.fastq.gz'
        read2 = f'{base_dir}/{SUBSPEC}/{sample}_2.fastq.gz'

        trimming = gwf.target_from_template(
            name = f'Trimming_{SUBSPEC}_{sample}',
            template=adapterremoval(
                input_fasta_pair1=read1,
                input_fasta_pair2=read2,
                sample_name = sample,
                out_dir=f'{outDir}/trimmed_reads',
                min_length=60
            )
        )

        ReadQC = gwf.target_from_template(
            name=f'ReadQC_{SUBSPEC}_{sample}',
            template=FastQC(
                input_fasta_pair1=trimming.outputs['pair1'],
                input_fasta_pair2=trimming.outputs['pair2'],
                sample_name = sample,
                out_dir=f'{outQC2}/FastQC'
            )
        )

        sampleID = "_".join(sample.split("_")[3:])  # Remove the first 3 chunks separated by "_"
        sampleIDs.add(sampleID) # add to set (ignores duplicates)
        DT = sample.split("_")[1] # Only keep the second chunk
        PU = "_".join([sample.split("_")[i] for i in (0, 2)]) # keep and join first and third chunk

        # add proper RG tags before use!
        map_paired = gwf.target_from_template(
            name=f'mapping_{SUBSPEC}_{sample}',
            template=mapping_pairdend(
                input_fasta_pair1=trimming.outputs['pair1'],
                input_fasta_pair2=trimming.outputs['pair2'],
                sample_name=sample,
                ref_genome='/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/bwa_indexed/Struthio_camelus_HiC.fasta',
                out_dir=outDir
            )
        )

        
        addRG = gwf.target_from_template(
            name=f'addRG_{SUBSPEC}_{sample}',
            template=add_RG(
                alignment_file=map_paired.outputs['alignment'],
                sample_name=sample,
                SM=sampleID, PU=PU, DT=DT,
                out_dir=outDir
            )
        )

    print(f"SampleIDs for {SUBSPEC}: {sampleIDs}")  # <-- Add this line to check

    for sampleID in sampleIDs:

        # Find .merged.bam files in outDir matching sampleID
        bam_files = [f for f in glob.glob(f"{outDir}/*{sampleID}*.RG.bam") if os.path.isfile(f)]
        if len(bam_files) >= 2:
            file1, file2 = bam_files[:2]
            merge_samples = gwf.target_from_template(
                name=f'merge_{SUBSPEC}_{sampleID}',
                template=merge_bam(
                    alignment_file1=file1,
                    alignment_file2=file2,
                    sampleID=sampleID,
                    out_dir=outDir
                )
            )
            input_bam = merge_samples.outputs['bam'] # to ensure the next step has an input
        elif len(bam_files) == 1:
            input_bam = bam_files[0] # if only one bam file, use it directly
        else:
            print(f"No BAM files found for {SUBSPEC} sampleID {sampleID}, skipping.")
            continue  # No BAM files found, skip

        mark_dups = gwf.target_from_template(
            name=f'dups_{SUBSPEC}_{sampleID}',
            template=mark_dups_samtools(
                alignment_file=merge_samples.outputs['bam'],
                sample_name=sampleID,
                out_dir=outDir
            )
        )

        pre_filter_stats = gwf.target_from_template(
            name=f'prefilter_stats_{SUBSPEC}_{sampleID}',
            template=samtools_stats(
                alignment_file=mark_dups.outputs['markdup'],
                sample_name=sampleID,
                out_dir=outDir
            )
        )

        filter_bam = gwf.target_from_template(
            name=f'filter_{SUBSPEC}_{sampleID}',
            template=samtools_filter(
                alignment_file=mark_dups.outputs['markdup'],
                sample_name=sampleID,
                out_dir=outDir,
                min_mq=20
            )
        )

        qc_bam = gwf.target_from_template(
            name=f'mapingQC_{SUBSPEC}_{sampleID}',
            template=qc_qualimap(
                alignment_file=filter_bam.outputs['filtered'],
                sample_name=sampleID,
                out_dir=outQC2
            )
        )
