#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def freebayes_population_set_workflow(configFile: str = glob.glob('*config.y*ml')[0]):
	"""
	Workflow: Create :format:`VCF` file for each sample in configuration.
	
	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""
	# --------------------------------------------------
	#                  Configuration
	# --------------------------------------------------
	
	CONFIG = yaml.safe_load(open(configFile))
	ACCOUNT: str = CONFIG['account']
	TAXONOMY: str | None = CONFIG['taxonomicGroup'].lower() if CONFIG['taxonomicGroup'] else None
	SPECIES_NAME: str = CONFIG['speciesName']
	REFERENCE_GENOME: str = CONFIG['referenceGenomePath']
	INTERGENIC_BED: str | None = CONFIG['intergenicBedFile']
	REPEATS_BED: str | None = CONFIG['repeatsBedFile']
	WORK_DIR: str = CONFIG['workingDirectoryPath'][:len(CONFIG['workingDirectoryPath']) - 1] if CONFIG['workingDirectoryPath'].endswith('/') else CONFIG['workingDirectoryPath']
	OUTPUT_DIR: str | None = (CONFIG['outputDirectoryPath'][:len(CONFIG['outputDirectoryPath']) - 1] if CONFIG['outputDirectoryPath'].endswith('/') else CONFIG['outputDirectoryPath']) if CONFIG['outputDirectoryPath'] else None
	PARTITION_SIZE: int = CONFIG['partitionSize'] if CONFIG['partitionSize'] else 100000
	FREEBAYES_SETTINGS: dict = CONFIG['freebayesSettings']
	FREEBAYES_PLOIDY: int | None = FREEBAYES_SETTINGS['samplePloidy'] if FREEBAYES_SETTINGS['samplePloidy'] else 100
	FREEBAYES_BESTN: int | None = FREEBAYES_SETTINGS['bestNAlleles'] if FREEBAYES_SETTINGS['bestNAlleles'] else 3
	FREEBAYES_MINALTFRC: float | int | None = FREEBAYES_SETTINGS['minAlternateFraction'] if FREEBAYES_SETTINGS['minAlternateFraction'] else 0
	FREEBAYES_MINALTCNT: int | None = FREEBAYES_SETTINGS['minAlternateCount'] if FREEBAYES_SETTINGS['minAlternateCount'] else 2
	VCF_MEM: int | None = FREEBAYES_SETTINGS['memory'] if FREEBAYES_SETTINGS['memory'] else 80
	VCF_TIME: str | None = str(FREEBAYES_SETTINGS['time']) if FREEBAYES_SETTINGS['time'] else '48:00:00'
	SAMPLE_LIST: list = CONFIG['sampleList']
	BATCHSETTINGS: dict = CONFIG['batchSettings']
	NBATCHES: int = BATCHSETTINGS['numberOfBatches'] if BATCHSETTINGS['numberOfBatches'] else 0
	BATCHNR: int | None = BATCHSETTINGS['currentBatchNumber'] if BATCHSETTINGS['currentBatchNumber'] else None

	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	)
	
	if os.path.exists(f'reference_partitions.{PARTITION_SIZE}bp.txt') and os.path.exists('reference_sequences.txt'):
		# If files exists reads data directly from files
		# Loads reference genome partitioning
		with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'r') as infile:
			partitions = [{'num': entry.split(sep='\t')[0].strip(), 'region': entry.split(sep='\t')[1].strip(), 'start': entry.split(sep='\t')[2].strip(), 'end': entry.split(sep='\t')[3].strip()} for entry in infile]
			nPadding = len(str(sum(1 for line in partitions)))
		# Loads list of contigs in reference genome
		with open('reference_sequences.txt', 'r') as infile:
			contigs = [{'contig': entry.split(sep='\t')[0].strip()} for entry in infile]
	else:
		# If files don't exist, generate data and write files
		sequences = parse_fasta(REFERENCE_GENOME)
		with open(f'reference_sequences.txt', 'w') as outfile:
			outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in sequences))
		# Partitions reference genome
		nPadding = padding_calculator(parseFasta=sequences, size=PARTITION_SIZE)
		partitions = partition_chrom(parseFasta=sequences, size=PARTITION_SIZE, nPad=nPadding)
		with open(f'reference_partitions.{PARTITION_SIZE}bp.txt', 'w') as outfile:
			outfile.write('\n'.join('\t'.join(str(i) for i in entry.values()) for entry in partitions))
		# Creates list of contigs in reference genome
		contigs = [{'contig': contig['sequence_name']} for contig in sequences]

	# When the argument line for any command becomes too long it cannot be executed. This can become an issue when jobs have too many dependencies.
	# For this workflow in can occur when cancatenating massively parallellised task, eg. create the VCF parts.
	# To alleviate the issue, multiple concatenation jobs will be created so that no concatenation has more than 5000 dependencies.
	nPartitions = len(partitions)
	segmentSize = 5000
	nSegments = int(round(nPartitions / segmentSize, 0) + 1) if (nPartitions / segmentSize > round(nPartitions / segmentSize, 0)) else int(round(nPartitions / segmentSize, 0))

	topDir = f'{WORK_DIR}/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}/vcf' if TAXONOMY else f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}/vcf'
	topOut = f'{OUTPUT_DIR}/vcf/{TAXONOMY.replace(" ", "_")}/{SPECIES_NAME.replace(" ", "_")}' if TAXONOMY else f'{OUTPUT_DIR}/vcf/{SPECIES_NAME.replace(" ", "_")}'

	setupDict = {group['groupName'].lower().replace(' ', '_'): {'name': group['groupName'].lower().replace(' ', '_'),
				  												'status': group['groupStatus'].lower() if group['groupStatus'] else 'i',
																'minDP': group['groupMinDP'],
																'highQuality': [sample for sample in group['highQualityBamFileList'] if sample],
																'allQuality': [sample for sample in group['highQualityBamFileList'] if sample] + [sample for sample in group['lowQualityBamFileList'] if sample]}
				for group in SAMPLE_LIST if group['highQualityBamFileList'] or group['lowQualityBamFileList']}

	indexReferenceGenome = gwf.target_from_template(
		name=f'index_reference_genome_{SPECIES_NAME.replace(" ", "_")}',
		template=index_reference_genome(
			referenceGenomeFile=REFERENCE_GENOME,
			outputDirectory=topDir
		)
	)

	# Small, half-baked job batch system in case the workflow would produce to many jobs to queue at once.
	# ------------------------------------------------------ #
	# Batched branch. Includes only generation of VCF parts. #
	# ------------------------------------------------------ #
	if NBATCHES != 0 and BATCHNR:
		batchSize = int(len(partitions)/NBATCHES)
		firstJob = (BATCHNR - 1) * batchSize
		lastJob = BATCHNR * batchSize
		if NBATCHES == BATCHNR:
			if len(partitions) - batchSize * NBATCHES != 0:
				lastJob += len(partitions) - batchSize * NBATCHES
		
		for group in setupDict:
			for sample in setupDict[group]['allQuality']:
				freebayesPartitionSingle = gwf.map(
					name=name_freebayes_partition_single,
					template_func=freebayes_partition_single,
					inputs=partitions[firstJob:lastJob],
					extra={'referenceGenomeFile': indexReferenceGenome.outputs['symlink'],
						   'bamFile': sample,
						   'outputDirectory': topDir,
						   'groupName': setupDict[group]['name'],
						   'sampleName': os.path.basename(os.path.dirname(sample)),
						   'ploidy': FREEBAYES_PLOIDY,
						   'bestNAlleles': FREEBAYES_BESTN,
						   'minAlternateFraction': FREEBAYES_MINALTFRC,
						   'minAlternateCount': FREEBAYES_MINALTCNT,
						   'memory': VCF_MEM,
						   'time': VCF_TIME}
				)

	# -------------------------------------- #
	# Non-batched branch. Includes all jobs. #
	# -------------------------------------- #
	else:
		# Create bed file of intergenic regions without repetitive regions
		if INTERGENIC_BED and REPEATS_BED:
			bedExcludeOverlapRepeats = gwf.target_from_template(
				name=f'intergenic_exluding_repeats_bed',
				template=bed_exclude_overlap(
					mainBedFile=INTERGENIC_BED,
					subtractionBedFile=REPEATS_BED,
					outputDirectory=topDir,
					speciesName=SPECIES_NAME
				)
			)

		elif INTERGENIC_BED and not REPEATS_BED:
			extractSoftmaskedIntervals = gwf.target_from_template(
				name=f'extract_repetitive_intervals',
				template=extract_softmasked_intervals(
					referenceGenomeFile=REFERENCE_GENOME,
					outputDirectory=topOut if OUTPUT_DIR else topDir
				)
			)

			bedExcludeOverlapRepeats = gwf.target_from_template(
				name=f'intergenic_exluding_repeats_bed',
				template=bed_exclude_overlap(
					mainBedFile=INTERGENIC_BED,
					subtractionBedFile=extractSoftmaskedIntervals.outputs['bed'],
					outputDirectory=topDir,
					speciesName=SPECIES_NAME
				)
			)

		for group in setupDict:
			# Create depth distribution and site counts
			if not setupDict[group]['highQuality'] == setupDict[group]['allQuality'] and setupDict[group]['highQuality']:
				depthDistributionHighQuality = gwf.target_from_template(
					name=f'depth_distribution_high_quality_{setupDict[group]['name']}',
					template=depth_distribution(
						bamFiles=setupDict[group]['highQuality'],
						minCoverageThreshold=setupDict[group]['minDP'],
						mode=0 if setupDict[group]['status'] == 'i' else 1,
						entryName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.highQuality',
						outputDirectory=topOut if OUTPUT_DIR else topDir,
						outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.highQuality.depth'
					)
				)

				sharedSitesWithinThresholdBedHighQuality = gwf.target_from_template(
					name=f'depth_threshold_bed_high_quality_{setupDict[group]['name']}',
					template=shared_sites_within_threshold_bed(
						bamFiles=setupDict[group]['highQuality'],
						depthDistributionTsv=depthDistributionHighQuality.outputs['tsv'],
						outputDirectory=topOut if OUTPUT_DIR else topDir,
						outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.highQuality'
					)
				)

				if INTERGENIC_BED:
					siteCountAllHighQuality = gwf.target_from_template(
						name=f'site_count_all_high_quality_{setupDict[group]['name']}',
						template=site_count_region(
							bamFiles=setupDict[group]['highQuality'],
							depthDistributionTsv=depthDistributionHighQuality.outputs['tsv'],
							bedFile=None,
							siteType='all',
							outputDirectory=topDir,
							outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.highQuality'
						)
					)

					siteCountIntergenicHighQuality = gwf.target_from_template(
						name=f'site_count_intergenic_high_quality_{setupDict[group]['name']}',
						template=site_count_region(
							bamFiles=setupDict[group]['highQuality'],
							depthDistributionTsv=depthDistributionHighQuality.outputs['tsv'],
							siteType='intergenic',
							outputDirectory=topDir,
							outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.highQuality',
							bedFile=INTERGENIC_BED
						)
					)

					siteCountRegionIntergenicExclRepeatsHighQuality = gwf.target_from_template(
						name=f'site_count_intergenic_excl_repeats_high_quality_{setupDict[group]['name']}',
						template=site_count_region(
							bamFiles=setupDict[group]['highQuality'],
							depthDistributionTsv=depthDistributionHighQuality.outputs['tsv'],
							siteType='intergenic_excl_repeats',
							outputDirectory=topDir,
							outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.highQuality',
							bedFile=bedExcludeOverlapRepeats.outputs['bed']
						)
					)

					mergeSiteTablesHighQuality = gwf.target_from_template(
						name=f'merge_site_tables_high_quality_{setupDict[group]['name']}',
						template=merge_site_tables(
							siteTables=[siteCountAllHighQuality.outputs['sitetable'],
										siteCountIntergenicHighQuality.outputs['sitetable'],
										siteCountRegionIntergenicExclRepeatsHighQuality.outputs['sitetable']],
							outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.highQuality',
							outputDirectory=topOut if OUTPUT_DIR else topDir
						)
					)

			depthDistributionAllQuality = gwf.target_from_template(
				name=f'depth_distribution_all_quality_{setupDict[group]['name']}',
				template=depth_distribution(
					bamFiles=setupDict[group]['allQuality'],
					minCoverageThreshold=setupDict[group]['minDP'],
					mode=0 if setupDict[group]['status'] == 'i' else 1,
					entryName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.allQuality',
					outputDirectory=topOut if OUTPUT_DIR else topDir,
					outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.allQuality.depth'
				)
			)

			sharedSitesWithinThresholdBedAllQuality = gwf.target_from_template(
				name=f'depth_threshold_bed_all_quality_{setupDict[group]['name']}',
				template=shared_sites_within_threshold_bed(
					bamFiles=setupDict[group]['allQuality'],
					depthDistributionTsv=depthDistributionAllQuality.outputs['tsv'],
					outputDirectory=topOut if OUTPUT_DIR else topDir,
					outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.allQuality'
				)
			)

			if INTERGENIC_BED:
				siteCountAllAllQuality = gwf.target_from_template(
					name=f'site_count_all_all_quality_{setupDict[group]['name']}',
					template=site_count_region(
						bamFiles=setupDict[group]['allQuality'],
						depthDistributionTsv=depthDistributionAllQuality.outputs['tsv'],
						bedFile=None,
						siteType='all',
						outputDirectory=topDir,
						outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.allQuality'
					)
				)

				siteCountIntergenicAllQuality = gwf.target_from_template(
					name=f'site_count_intergenic_all_quality_{setupDict[group]['name']}',
					template=site_count_region(
						bamFiles=setupDict[group]['allQuality'],
						depthDistributionTsv=depthDistributionAllQuality.outputs['tsv'],
						siteType='intergenic',
						outputDirectory=topDir,
						outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.allQuality',
						bedFile=INTERGENIC_BED
					)
				)

				siteCountRegionIntergenicExclRepeatsAllQuality = gwf.target_from_template(
					name=f'site_count_intergenic_excl_repeats_all_quality_{setupDict[group]['name']}',
					template=site_count_region(
						bamFiles=setupDict[group]['allQuality'],
						depthDistributionTsv=depthDistributionAllQuality.outputs['tsv'],
						siteType='intergenic_excl_repeats',
						outputDirectory=topDir,
						outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.allQuality',
						bedFile=bedExcludeOverlapRepeats.outputs['bed']
					)
				)

				mergeSiteTablesAllQuality = gwf.target_from_template(
					name=f'merge_site_tables_all_quality_{setupDict[group]['name']}',
					template=merge_site_tables(
						siteTables=[siteCountAllAllQuality.outputs['sitetable'],
									siteCountIntergenicAllQuality.outputs['sitetable'],
									siteCountRegionIntergenicExclRepeatsAllQuality.outputs['sitetable']],
						outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.allQuality',
						outputDirectory=topOut if OUTPUT_DIR else topDir
					)
				)

			# Create vcf parts, concatenate parts into whole vcfs, normalize vcfs, and merge vcfs
			vcfSingleListHighQuality = []
			vcfSingleListAllQuality = []
			for sample in setupDict[group]['allQuality']:
				sampleName = os.path.basename(os.path.dirname(sample))
				freebayesPartitionSingle = gwf.map(
					name=name_freebayes_partition_single,
					template_func=freebayes_partition_single,
					inputs=partitions,
					extra={'referenceGenomeFile': indexReferenceGenome.outputs['symlink'],
						   'bamFile': sample,
						   'outputDirectory': topDir,
						   'groupName': setupDict[group]['name'],
						   'sampleName': sampleName,
						   'ploidy': FREEBAYES_PLOIDY,
						   'bestNAlleles': FREEBAYES_BESTN,
						   'minAlternateFraction': FREEBAYES_MINALTFRC,
						   'minAlternateCount': FREEBAYES_MINALTCNT,
						   'memory': VCF_MEM,
						   'time': VCF_TIME}
				)

				if nSegments <= 1:
					concatenateFreebayesSingle = gwf.target_from_template(
						name=f'concatenate_freebayes_vcf_{setupDict[group]['name']}_{sampleName.replace("-", "_")}',
						template=concat_vcf(
							files=collect(freebayesPartitionSingle.outputs, ['vcf'])['vcfs'],
							outputName=f'{sampleName}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall',
							outputDirectory=f'{topOut}/{setupDict[group]['name']}/{sampleName}' if OUTPUT_DIR else f'{topDir}/raw_vcf/{setupDict[group]['name']}/{sampleName}',
							compress=True
						)
					)

				else:
					segmentList = []
					start = 0
					end = segmentSize
					collection = collect(freebayesPartitionSingle.outputs, ['vcf'])['vcfs']

					for i in range(nSegments):
						concatenateFreebayesSingleSegment = gwf.target_from_template(
							name=f'concatenate_freebayes_vcf_{setupDict[group]['name']}_{sampleName.replace("-", "_")}_segment_{i+1}',
							template=concat_vcf(
								files=collection[start : end],
								outputName=f'{sampleName}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}.segment{i+1}',
								outputDirectory=f'{topDir}/raw_vcf/{setupDict[group]['name']}/{sampleName}/tmp',
								compress=True
							)
						)

						segmentList.append(concatenateFreebayesSingleSegment.outputs['concat_file'])
						if i < nSegments - 1:
							start = end
							end += segmentSize
						elif i == nSegments - 1:
							start = end
							end = nPartitions

					concatenateFreebayesSingle = gwf.target_from_template(
						name=f'concatenate_freebayes_vcf_{setupDict[group]['name']}_{sampleName.replace("-", "_")}_complete',
						template=concat_vcf(
							files=segmentList,
							outputName=f'{sampleName}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall',
							outputDirectory=f'{topOut}/{setupDict[group]['name']}/{sampleName}' if OUTPUT_DIR else f'{topDir}/raw_vcf/{setupDict[group]['name']}/{sampleName}',
							compress=True
						)
					)

				# normVcfSingle = gwf.target_from_template(
				# 	name=f'normalize_vcf_{setupDict[group]['name']}_{sampleName.replace("-", "_")}',
				# 	template=norm_vcf(
				# 		vcfFile=concatenateFreebayesSingle.outputs['concat_file'],
				# 		referenceGenomeFile=indexReferenceGenome.outputs['symlink'],
				# 		outputName=f'{sampleName}.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall',
				# 		outputDirectory=f'{topDir}/raw_vcf/{setupDict[group]['name']}/{sampleName}'
				# 	)
				# )
					
				# # Collect concatenated single population VCF files
				# if sample in setupDict[group]['highQuality']:
				# 	vcfSingleListHighQuality.append(normVcfSingle.outputs['vcf'])
				# vcfSingleListAllQuality.append(normVcfSingle.outputs['vcf'])

				# Collect concatenated single population VCF files
				if sample in setupDict[group]['highQuality']:
					vcfSingleListHighQuality.append(concatenateFreebayesSingle.outputs['concat_file'])
				vcfSingleListAllQuality.append(concatenateFreebayesSingle.outputs['concat_file'])

			# if not setupDict[group]['highQuality'] == setupDict[group]['allQuality'] and setupDict[group]['highQuality']:				
			# 	if len(vcfSingleListHighQuality) > 1:
			# 		mergeVcfSingleHighQuality = gwf.target_from_template(
			# 			name=f'merge_vcf_single_high_quality_{setupDict[group]['name']}',
			# 			template=merge_vcf_no_duplicates(
			# 				vcfFiles=vcfSingleListHighQuality,
			# 				outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.highQuality.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall.norm',
			# 				outputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/raw_vcf/{setupDict[group]['name']}'
			# 			)
			# 		)

			# if len(vcfSingleListAllQuality) > 1:
			# 	mergeVcfSingleAllQuality = gwf.target_from_template(
			# 		name=f'merge_vcf_single_all_quality_{setupDict[group]['name']}',
			# 		template=merge_vcf_no_duplicates(
			# 			vcfFiles=vcfSingleListAllQuality,
			# 			outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.allQuality.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall.norm',
			# 			outputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/raw_vcf/{setupDict[group]['name']}'
			# 		)
			# 	)
			
			if not setupDict[group]['highQuality'] == setupDict[group]['allQuality'] and setupDict[group]['highQuality']:				
				if len(vcfSingleListHighQuality) > 1:
					mergeNormRmdupVcfSingleHighQuality = gwf.target_from_template(
						name=f'merge_norm_rmDup_vcf_single_high_quality_{setupDict[group]['name']}',
						template=merge_norm_no_duplicates_vcf(
							vcfFiles=vcfSingleListHighQuality,
							referenceGenomeFile=indexReferenceGenome.outputs['symlink'],
							outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.highQuality.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall',
							outputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/raw_vcf/{setupDict[group]['name']}'
						)
					)

			if len(vcfSingleListAllQuality) > 1:
				mergeNormRmdupVcfSingleAllQuality = gwf.target_from_template(
					name=f'merge_norm_rmDup_vcf_single_all_quality_{setupDict[group]['name']}',
					template=merge_norm_no_duplicates_vcf(
						vcfFiles=vcfSingleListAllQuality,
						referenceGenomeFile=indexReferenceGenome.outputs['symlink'],
						outputName=f'{species_abbreviation(SPECIES_NAME)}.{setupDict[group]['name']}.{'ingroup' if setupDict[group]['status'] == 'i' else 'outgroup'}.allQuality.freebayes_n{FREEBAYES_BESTN}_p{FREEBAYES_PLOIDY}_minaltfrc{FREEBAYES_MINALTFRC}_minaltcnt{FREEBAYES_MINALTCNT}_singlecall',
						outputDirectory=f'{topOut}/{setupDict[group]['name']}' if OUTPUT_DIR else f'{topDir}/raw_vcf/{setupDict[group]['name']}'
					)
				)

	print(f'Intermediary files will be place at: {topDir}/')
	print(f'Output files will be placed at: {topOut if OUTPUT_DIR else topDir}/')
	
	return gwf
