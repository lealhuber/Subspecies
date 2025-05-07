# Subspecies project
This project is part of the ChamberTempRNA project, but because it involves very different analyses I am treating it as its own. We have individual and pooled WGS data from three subspecies of ostriches that used to be raised on the Oudtshoorn research farm. The aims of this project are
1) finding genes related to thermal adaptation by identifying signatures of selection between the subspecies
2) identifying promoter regions of the candidate genes from the ChamberTempRNA project and looking if there are variants
For this, we will do variant calling and subsequently do an Fst analysis
First everything for the individual sequencing.
## Mapping
First I removed adapters and trimmed reads using TrimGalore, discarding reads with a Phred score below 20 and removing Ns. Then reads were mapped using bwa mem with default settings, adding sample names as read groups. Using samtools, duplicates were marked and, after some statistics, discarded. I also filtered out reads with a MQ below 20. QC results can be found in the respective subspecies folders red, black and blue in mapping/outputs/.
## Variant calling
Joint calling with Nathalies script. It splits into scaffolds and the scaffolds into 1 MB chunks and then calls variants with freebayes, also keeping monomorphic sites. The --populations option was used to inform freebayes of subspecies membership. Then the vcf files were concatenated into one file for statistics and filtering.
## Variant filtering
Variants were filtered for basic best practice quality of 30, a mean depth of at least a third of average depth (--> x) and at most twice average depth (--> x). Then...
## Fst analysis