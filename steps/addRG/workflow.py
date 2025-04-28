## declare workflow
from gwf import Workflow # type: ignore
import os

gwf = Workflow(defaults={'account': 'ostrich_thermal'})

#The path to where the folder containing .bam files (files can be in subfolders)
bam_path =f'/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/mapping/outputs/'
out_dir = '/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/addRG/outputs'

bam_list = []
# Iterate over the directory tree
for root, dirs, files in os.walk(bam_path):
    # Check each file in the current directory
    for file in files:
        if file.endswith('filtered.bam'):
            # Append the absolute path of the file to the list, but only of the bam file
            bam_list.append(os.path.join(root, file))


# # Sort the bam list by the ID
sorted_bam_list =  sorted(bam_list, key=lambda x: x[-22:])
# print(f'bam files with path: {sorted_bam_list}')

for sample in sorted_bam_list:
    sample_name = os.path.basename(sample).split('.')[0]
    # print(f'sample names: {sample_name}')
    gwf.target(f'addRG_{sample_name}', #name of the target
           cores=4,
           memory='16gb',
           walltime='06:00:00',	
           inputs= {'bamfile': sample}, 
           outputs=[f'{out_dir}/{sample_name}.RG.filtered.bam',
                    f'{out_dir}/{sample_name}.RG.filtered.bam.bai']) << '''
    samtools addreplacerg \\
        --threads 4 \\
        -r ID:{sample_name} \\
        -r SM:{sample_name} \\
        {sample} \\
        | samtools sort \\
            --threads 4 \\
            -o {out_dir}/{sample_name}.RG.filtered.bam -
        samtools index {out_dir}/{sample_name}.RG.filtered.bam
    '''.format(sample=sample, sample_name=sample_name, out_dir=out_dir)
    