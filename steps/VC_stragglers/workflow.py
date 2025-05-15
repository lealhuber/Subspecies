#!/bin/env python3
from gwf import Workflow, AnonymousTarget # type: ignore
import glob
from templates import * # type: ignore
import os

#Accunt details
gwf = Workflow(defaults={'account': 'ostrich_thermal'})


#______GENERAL VARIABLES___CHANGEABLE! ___#___ OBS STORTING JOB MOVED TO FILER_WORKFLOW

#Change pr spp.
#Species (spp) information
subsp = 'Sc_subsp'

###########################################
#Reference genome REPLACE PR SPP !!!!! 
###########################################

reference_genome='/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/bwa_indexed/Struthio_camelus_HiC.fasta'

#The path to where the folder containing .bam files (files can be in subfolders)
bam_path =f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/addRG/outputs/'

#Folders for pipeline outputs. If they do not exsist, they will be created.  

#path to Where you want the temp files
temps =f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/VC_stragglers/temp/'

#path to where logfiles shoulde be outputtet 
out_folder =f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/VC_stragglers/logs/'

#RESULTS folder, where should the population VCF be outputtet 
results =f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/VC_stragglers/outputs/'

#path to where the stats should be
stats_path=f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/VC_stragglers/stats'


#____CREATED VAIABLES__DO NOT CHANGE!!  --- UNLESS structure changes!!! ____# 

#_____________________________________________ DONT CHANGE !!!! 
#OUTPUT directories: *change for new USER
#Check if present - if not create: Temporary folder (and spp year folder) + Results folders(),    

###### dir structure !  -- Should not be changed if one wishes to use the filterVCF workflow :D

#Temporary dir
chr_vcf_path = temps
# Check if the directory exists, and create it if not
os.makedirs(chr_vcf_path, exist_ok=True)
# print(f"Directory '{chr_vcf_path}' created or already exists.")

#Results dir
# Check if the directory exists, and create it if not
os.makedirs(results, exist_ok=True)
# print(f"Directory '{results}' created or already exists.")


#out dir
# Check if the directory exists, and create it if not
os.makedirs(out_folder, exist_ok=True)
# print(f"Directory '{out_folder}' created or already exists.")

# Stats dir
# Check if the directory exists, and create it if not
os.makedirs(stats_path, exist_ok=True)
# print(f"Directory '{stats_path}' created or already exists.")


# Bam: Use glob to find matching files - Create a bam list path and file. 
# Search all subdirectories in the given directory for matching bams 
# search_pattern = os.path.join(bam_path, '**', bam)
# bam_list = glob.glob(search_pattern, recursive=True)
 
#Above is the prior way of creating the list. 
#This does not work in all cases, therefore I have made the below loop - if the structure changes this should change aswell. 

# ######Create a list with bam_paths          ------- OBS the nameing! some are called _filtered some are called .filtered
# # Initialize an empty list to store the paths to all bam files  
bam_list = []
# Iterate over the directory tree
for root, dirs, files in os.walk(bam_path):
    # Check each file in the current directory
    for file in files:
        # Check if the file ends with *.bam
        if file.endswith(('filtered.bam', 'filtered.bam.bai')):
            # I want the samples to have the subspecies name in front so they will be easier later to distinguish, so I rename them here
            # Get the parent directory name (black, blue, red)
            # but because I ran it again like an ass I will break and fix it now
            parent_dir = os.path.basename(root)
            # Check if the file already starts with the parent directory name to avoid double-prefixing
            if file.startswith("_"):
                new_name = file[1:]  # Remove only the first character
                old_path = os.path.join(root, file)
                new_path = os.path.join(root, new_name)
                os.rename(old_path, new_path)
                file = new_name  # Update the file variable so the renamed path is appended
            if file.endswith('filtered.bam'):
                # Append the absolute path of the file to the list, but only of the bam file
                bam_list.append(os.path.join(root, file))


# # Sort the bam list by the ID
sorted_bam_list =  sorted(bam_list, key=lambda x: x[-22:])

# Full path for the bam.list file 
path_bam_list = os.path.join(out_folder, '{}_bam.list'.format(subsp))

bam_files=sorted_bam_list

# print(f'bam list: ,{bam_files}') # to check if all there

# getting samle names. 
samples = [path.split('/')[-1].split('.')[0] for path in sorted_bam_list] # this will break easily but for now it should work
# print(f'Sample names: , {samples}')


# Example BAM file path
bam_file_path = sorted_bam_list[0]

## make populations file
with open(f'{out_folder}populations_file.txt', "w") as f:
    for sample in samples:
        prefix = sample.split("_")[0]
        f.write(f"{sample}\t{prefix}\n")

populations = f'{out_folder}populations_file.txt'

#Get the CHR names from the header of the reference genome

with open(reference_genome, 'r') as file:
    chromosomes = [line[1:].strip() for line in file if line.startswith('>')]


# Writing the called CHR to SPP_CHR.list
with open('{}{spp}_CHR.list'.format(out_folder, spp=subsp), 'w') as file:
        file.write('\n'.join(chromosomes))
# print(f'Chromosomes: {chromosomes}')

## calculate start and endpoints for partitioning chromosomes from reference genome
PARTITION_SIZE = 1000000
sequences = parse_fasta(reference_genome)
with open(f'{out_folder}reference_sequences.txt', 'w') as outfile:
	outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
# Partitions reference genome
nPadding = padding_calculator(parseFasta=sequences, size=PARTITION_SIZE)
partitions = partition_chrom(parseFasta=sequences, size=PARTITION_SIZE, nPad=nPadding)
selected_numbers = ['0063', '0069'] # choose the chunks to re-run
partitions = [entry for entry in partitions if entry["num"] in selected_numbers]
with open(f'{out_folder}reference_partitions.{PARTITION_SIZE}bp.txt', 'w') as outfile:
    outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
print(f'Partitions: {partitions}')




#####################################################################################################
#___________________________________________________________________________________________________#

#-----------------------------4 jobs from bams to the VCF -------------------------------------------

#___________________________________________________________________________________________________#
#####################################################################################################


#gwf 1 job: Call CHR VCFs from bam files.
all_vcfs = []
for part in partitions:
    target_name = 'call_vcf_' + subsp + part['num'] # Construct the unique target name
    this_vcf = subsp + '_' + part['num'] + '.vcf'
    full_path_vcf = os.path.join(temps, this_vcf)
    all_vcfs.append(full_path_vcf)
    job1 = gwf.target_from_template(
        name=target_name,
        template=freebayes_CHR_vcf(
            files=bam_files,
            temp_path=temps,
            output_name=this_vcf,
            reference_genome=reference_genome,
            population_match=populations,
            region=part['region'],
            start=part['start'],
            end=part['end']
        )
    )

print(f"All VCFs: {all_vcfs}")

