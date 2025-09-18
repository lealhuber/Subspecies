#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=12:00:00

# run in steps/fst/scripts/

# make geno file
python simonmartin_genomics/VCF_processing/parseVCFs.py -i /faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/outputs/input.vcf.gz \
    --threads 8 -o ../temp/input.geno.gz
# loop through different window sizes
seq = (10000 25000 50000 100000)
for w in ${seq[@]}; do
    step=$(($w/2))
    echo "Window size: $w, step size: $step"
    python simonmartin_genomics/popgenWindows.py -w $w -s $step -m 5000 -g ../temp/input.geno.gz -o ../outputs/w${w}.csv.gz -f phased -T 8 --popsFile populations_file.txt
done