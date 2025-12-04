#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=8
#SBATCH --mem=5g
#SBATCH --time=12:00:00

# run in steps/fst/scripts/
vcf_ini=/faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/outputs/allpops.HWE.merged.vcf.gz
# index if not already indexed
# bcftools index --threads 8 $vcf_ini
# echo "VCF indexed"
# make geno file
python simonmartin_genomics/VCF_processing/parseVCFs.py -i $vcf_ini \
    --threads 8 -o ../temp/atp2.HWE.input.geno.gz
echo "Geno file created"
# loop through different window sizes
seq=(10000 25000 50000 100000)
for w in ${seq[@]}; do
    # step=$(($w/2))
    # echo "Window size: $w, step size: $step"
    python simonmartin_genomics/popgenWindows.py -w $w -m 500 -g ../temp/atp2.HWE.input.geno.gz -o ../outputs/atp2.w_$w.csv.gz -f phased -T 8 -p black -p blue -p red --popsFile populations_file.txt
done
echo "All window sizes completed"
echo "Ended at $(date)"