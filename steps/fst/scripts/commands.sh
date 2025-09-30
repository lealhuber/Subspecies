#!/bin/bash
#SBATCH --account=ostrich_thermal
#SBATCH --cpus-per-task=8
#SBATCH --mem=5g
#SBATCH --time=12:00:00

# run in steps/fst/scripts/
# index if not already indexed
bcftools index --threads 8 /faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/outputs/all.filtered.monobiallelic.snp.vcf.gz
# make geno file
python simonmartin_genomics/VCF_processing/parseVCFs.py -i /faststorage/project/ostrich_thermal/people/leah/Subspecies/steps/variant_filtering/outputs/all.filtered.monobiallelic.snp.vcf.gz \
    --threads 8 -o ../temp/atp2.HWE.input.geno.gz
echo "Geno file created"
# loop through different window sizes
seq=(10000 25000 50000 100000)
for w in ${seq[@]}; do
    # step=$(($w/2))
    # echo "Window size: $w, step size: $step"
    python simonmartin_genomics/popgenWindows.py -w $w -m 500 -g ../temp/atp1.noHWE.input.geno.gz -o ../outputs/atp2.w_$w.csv.gz -f phased -T 8 -p black -p blue -p red --popsFile populations_file.txt
done