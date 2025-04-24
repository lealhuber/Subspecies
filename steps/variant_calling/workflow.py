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
subsp = 'red'

###########################################
#Reference genome REPLACE PR SPP !!!!! 
###########################################

reference_genome='/faststorage/project/ostrich_thermal/BACKUP/ostrich_reference/Struthio_camelus_HiC/bwa_indexed/Struthio_camelus_HiC.fasta'

#The path to where the folder containing .bam files (files can be in subfolders)
bam_path =f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/mapping/outputs/{subsp}'

#Folders for pipeline outputs. If they do not exsist, they will be created.  

#path to Where you want the temp files
temps =f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_calling/{subsp}/temp/'

#path to where logfiles shoulde be outputtet 
out_folder =f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_calling/{subsp}/logs/'

#RESULTS folder, where should the population VCF be outputtet 
results =f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_calling/{subsp}/vcf_outputs/'

#path to where the states should be
stats_path=f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_calling/{subsp}/stats'


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
        if file.endswith('filtered.bam'):
            # Append the absolute path of the file to the list
            bam_list.append(os.path.join(root, file))


# # Sort the bam list by the ID
sorted_bam_list =  sorted(bam_list, key=lambda x: x[-22:])

# Full path for the bam.list file 
path_bam_list = os.path.join(out_folder, '{}_bam.list'.format(subsp))

bam_files=sorted_bam_list


####   #writing the above bam list to a file - deletes any previously deleted file
####   OBS If parts needs to be run again (e.g last two jobs) ! this needs to be masked out 
####   - for a full run it dosent matter!!  

# # Check if the BAM file exists
# if os.path.exists(bam_files):
#     print("bamlist found. Deleting...")
#     os.remove(bam_files)  # Delete the file if it exists

# # Write each element on a separate line
# with open(bam_files, 'w') as file:
#     for i, element in enumerate(sorted_bam_list):
#         file.write(element)
#         if i < len(sorted_bam_list) - 1:
#             file.write('\n')

# print("New bamlist file created.")


# getting samle names. 
samples = [path.split('/')[-1].split('.')[0] for path in sorted_bam_list] # this will break easily but for now it should work


# Example BAM file path
bam_file_path = sorted_bam_list[0]

#Get the CHR names from the header of the reference genome

with open(reference_genome, 'r') as file:
    chromosomes = [line[1:].strip() for line in file if line.startswith('>')]


# Writing the called CHR to SPP_CHR.list
with open('{}{spp}_CHR.list'.format(out_folder, spp=subsp), 'w') as file:
        file.write('\n'.join(chromosomes))



#####################################################################################################
#___________________________________________________________________________________________________#

#-----------------------------4 jobs from bams to the VCF -------------------------------------------

#___________________________________________________________________________________________________#
#####################################################################################################



#gwf 1 job: Call CHR VCFs from bam files.
all_vcfs = []
for i, chr in enumerate(chromosomes):
    target_name = 'call_vcf_' + subsp + chr # Construct the unique target name
    this_vcf = subsp + '_' + chr + '.vcf'
    full_path_vcf = os.path.join(temps, this_vcf)
    all_vcfs.append(full_path_vcf)
    job1 = gwf.target_from_template(
        name=target_name,
        template=freebayes_CHR_vcf(
            files=bam_files,
            temp_path=temps,
            output_name=this_vcf,
            chromosome=chr,
            reference_genome=reference_genome
        )
    )


#DEBUG HERE !! 
# Add this print statement to check the files in the list before the second job 
# removed because annoying :D 
# print(f"All VCFs: {all_vcfs}")

#gwf 2 job: concat CHR vcf -> genome vcf 
job2 = gwf.target_from_template(
    name='concat_CHR_vcf' + subsp,
    template=concat_CHR_vcf(
        files=all_vcfs,
        species_name=subsp,
        path=results
        )
    )

# gwf 3 job: sort file by samle list, and chr_pos sorting
job3 = gwf.target_from_template(
    name='sort_vcf' + subsp,
    template=Sort_vcf(
        vcf=job2.outputs['genome_vcf'],
        species_name=subsp,
        samples=samples, 
        path=results
        )
    )

# gwf 4 job: zip vcf
job4 = gwf.target_from_template(
    name='zip_vcf' + subsp,
    template=zip_file(
        file=job3.outputs['sorted_genome_vcf'],
        species_name=subsp,  
        )
    )


# gwf 5 job: Raw vcf stats pr sample
all_stats = []
for i, sample in enumerate(samples):
    this_stat='{path}{sample}.stats'.format( \
            path=stats_path, sample=sample)
    all_stats.append(this_stat)
    target_name = 'vcf_stats_' + sample  # Construct the unique target name
    job5 = gwf.target_from_template(
        name=target_name,
        template=vcf_stats(
            file=job4.outputs['vcf_gz'],
            stats_path=stats_path,
            sample=sample, 
            this_stat=this_stat  # Use the sample from list   
        )
    )



#gwf 6 job: Multiqc genome_stats, individual stats 

job6 = gwf.target_from_template(
    name='multiqc' + subsp,
    template=multiqc_vcf_stats(
        files=all_stats,
        species_name=subsp, 
        result_path=results
    )
)
