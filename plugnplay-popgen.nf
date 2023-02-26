#!/usr/bin/env nextflow

nextflow.enable.dsl = 2




// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	// Input channels
	ch_vcfs = Channel
		.fromPath ( params.samplesheet )
		.splitCsv ( header: true )
		.map { row -> tuple( row.sample_id, row.species, row.population, row.prep_type, row.platform, file(row.raw_vcf) ) }

	// Workflow steps
	SNP_QUAL_FILTERING (
		ch_vcfs
	)
	
	// SINGLE_VCF_STATS (
	// 	SNP_QUAL_FILTERING.out.vcf
	// )
	
	MERGE_VCFS (
		SNP_QUAL_FILTERING.out.vcf
			.map { sample_id, species, population, prep_type, platform, raw_vcf -> tuple( file(raw_vcf), species, prep_type ) }
			.filter { file(it[0]).countLines() > 0 }
			.groupTuple( by: [1,2] ),
		SNP_QUAL_FILTERING.out.index.collect()
	)
	
	FILTER_SNPS_BY_SAMPLE_COUNT {
		MERGE_VCFS.out.vcf
	}
	
	RECORD_SNP_FILTERS {
		FILTER_SNPS_BY_SAMPLE_COUNT.out.vcf
	}
	
	// MULTI_VCF_STATS (
	// 	FILTER_SNPS_BY_SAMPLE_COUNT.out.vcf
	// )
	
	MAKE_POP_MAP (
		ch_vcfs
			.map { sample_id, species, population, prep_type, platform, raw_vcf -> 
				sample_id.replaceAll('0', '') + "_" + sample_id.replaceAll('0', '') + ", ${species}, ${population}" }
			.collect()
			.collectFile( name: "pop_map.txt", newLine: true ),
		FILTER_SNPS_BY_SAMPLE_COUNT.out.vcf
			.map { vcf, species, prep, sample_size -> vcf }
			.collect()
	)
	
	SELECT_SFS_PROJECTION (
		FILTER_SNPS_BY_SAMPLE_COUNT.out.vcf,
		MAKE_POP_MAP.out
	)
	
	SFS_ESTIMATION (
		FILTER_SNPS_BY_SAMPLE_COUNT.out.vcf,
		MAKE_POP_MAP.out,
		SELECT_SFS_PROJECTION.out
	)
	
	VISUALIZE_SFS (
		SFS_ESTIMATION.out.sfs
	)
	
	BUILD_STAIRWAY_PLOT_BLUEPRINT (
		SFS_ESTIMATION.out.sfs
	)
	
	BUILD_STAIRWAY_PLOT_SCRIPT (
		BUILD_STAIRWAY_PLOT_BLUEPRINT.out
	)
	
	STAIRWAY_PLOT ( 
		BUILD_STAIRWAY_PLOT_SCRIPT.out
	)
	
	// POP_STRUCTURE_PCA ( )
	// 
	// ADMIXTURE_PLOT ( )

}
// --------------------------------------------------------------- //




// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //

// working with single nucleotide polymorphisms
params.snp_results = params.results + "/" + "1_per_sample_SNPs"
params.single_sample_snp_stats = params.snp_results + "/" + "snp_stats"
params.merged_snps = params.results + "/" + "2_per_species_SNPs"
params.multisample_snp_stats = params.merged_snps + "/" + "snp_stats"

// Additional, optional analyses
params.analyses = params.results + "/" + "3_analyses"
params.sfs_results = params.analyses + "/" + "site_frequency_spectra"
params.angsd_results = params.analyses + "/" + "ANGSD_outputs"
params.angsd_sfs_plots = params.angsd_results + "/" + "SFS_plots"
params.stairway_plots = params.analyses + "/" + "Stairway_plots"
params.stairway_plot_scripts = params.stairway_plots + "/" + "stairwayplot_scripts"
params.admixture_results = params.analyses + "/" + "admixture_plots"
params.pca_results = params.analyses + "/" + "PCA_plots"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //


process SNP_QUAL_FILTERING {
	
	tag "${sample}"
	publishDir params.snp_results, pattern: "*.vcf.gz", mode: 'copy'
	publishDir params.snp_results, pattern: "*.tbi", mode: 'copy'
	
	cpus 2
	
	input:
	tuple val(sample), val(species), val(pop), val(prep), val(platform), path(vcf)
	
	output:
	tuple val(sample), val(species), val(pop), val(prep), val(platform), path("${sample}_${prep}_filtered.vcf.gz"), emit: vcf
	path "*.tbi", emit: index
	
	script:
	"""
	
	vcftools --gzvcf ${vcf} \
	--max-alleles 2 \
	--minQ ${params.min_quality} \
	--minDP ${params.min_depth} \
	--remove-indels \
	--remove-filtered-all \
	--recode --stdout \
	| bgzip -c > "${sample}_${prep}_filtered.vcf.gz" && \
	tabix -p vcf "${sample}_${prep}_filtered.vcf.gz"
	
	"""
	
}

process SINGLE_VCF_STATS {

	publishDir params.single_sample_snp_stats, mode: 'copy'
	
	input:
	tuple val(sample), val(species), val(pop), val(prep), val(platform), path(filtered_vcf)
	
	output:
	path "*.stats"
	
	script:
	"""
	
	bcftools stats \
	-F ${params.reference} \
	-s - ${filtered_vcf} > ${filtered_vcf}.stats
	
	"""

}


process MERGE_VCFS {
	
	tag "${prep} ${species}"
	
	cpus 8
	
	input:
	tuple path(vcfs), val(species), val(prep)
	path vcf_indices
	
	output:
	tuple path("*.vcf.gz"), val(species), val(prep), emit: vcf
	
	shell:
	'''
	
	bcftools merge \
	--merge snps \
	--output-type z \
	--threads !{task.cpus} \
	--output !{species}_multisample_!{prep}.vcf.gz \
	*!{species}*.vcf.gz
	
	'''
	
}


process FILTER_SNPS_BY_SAMPLE_COUNT {
	
	tag "${prep} ${species}"
	publishDir params.merged_snps, pattern: "*.vcf.gz", mode: 'copy'
	publishDir params.merged_snps, pattern: "*.tbi", mode: 'copy'
	
	cpus 2
	
	input:
	tuple path(vcf), val(species), val(prep)
	
	output:
	tuple path("${species}_${prep}_filtered.vcf.gz"), val(species), val(prep), env(sample_size), emit: vcf
	path "*.tbi", emit: index
	
	script:
	"""
	
	vcftools --gzvcf ${vcf} \
	--max-alleles 2 \
	--max-missing-count ${params.max_snp_missingness} \
	--remove-indels \
	--remove-filtered-all \
	--recode --stdout \
	| bgzip -c > "${species}_${prep}_filtered.vcf.gz" && \
	tabix -p vcf "${species}_${prep}_filtered.vcf.gz" && \
	sample_size=`bcftools query -l ${species}_${prep}_filtered.vcf.gz | wc -l`
	
	"""
	
}


process RECORD_SNP_FILTERS {
	
	tag "${prep} ${species}"
	
	publishDir params.merged_snps, mode: 'copy', pattern: '*.txt'
	if ( params.stairwayplot == true ){ publishDir params.sfs_results, mode: 'copy', pattern: '*.txt' }
	
	cpus 1
	
	input:
	tuple path(vcf), val(species), val(prep), val(sample_size)
	
	output:
	path "*.txt"
	
	shell:
	'''
	sample_ids=`bcftools query -l !{vcf}`
	minInd="$((!{sample_size} - !{params.max_snp_missingness}))"
	
	touch !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "VCF FILTER SETTINGS APPLIED TO !{prep} !{species} SNPS" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "----------------------------------------------------------" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minor allele frequency: !{params.minor_allele_frequency}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum samples: ${minInd} / !{sample_size}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum variant quality score: !{params.min_quality}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum Depth: !{params.min_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Maximum Depth: !{params.max_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Sample IDs included: " >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo ${sample_ids} >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	'''
	
}


process MULTI_VCF_STATS {
	
	publishDir params.multisample_snp_stats, mode: 'copy'
	
	input:
	tuple path(vcf), val(species), val(prep), val(sample_size)
	
	output:
	path "*.stats"
	
	script:
	"""
	
	bcftools stats \
	-F ${params.reference} \
	-s - ${vcf} > ${vcf}.stats
	
	"""
	
}


process MAKE_POP_MAP {
	
	input:
	path pop_map
	path merged_vcfs
	
	output:
	path "*_pop_map.txt"
	
	script:
	if ( params.whole_species_mode == true ){
		"""
		# Create a list of different species in this dataset. Each species will get its own pop map.
		# Then, take the first two columns of pop map, which are sample ID and species ID.
		# Finally, use for loops to make a population map for each species, and check that each sample
		# is actually present in the VCF.
		cut -d',' -f2 ${pop_map} | sort | uniq > species.txt && \
		cut -f 1,2 -d',' ${pop_map} | tr ', ' '\t' > tmp_pop_map1.txt && \
		for i in `cat species.txt`;
		do
			bcftools query -l \${i}_*_filtered.vcf.gz > \${i}_samples.txt && \
			grep -w \$i tmp_pop_map1.txt > tmp_pop_map2.txt && \
			touch \${i}_pop_map.txt && \
			for j in `cat \${i}_samples.txt`;
			do
				if grep -qw \$j tmp_pop_map2.txt; then
					grep -w \$j tmp_pop_map2.txt >> \${i}_pop_map.txt
				fi
			done
		done
		"""
	} else {
		"""
		# Create a list of different species in this dataset. Each species will get its own pop map.
		# Then, take the first and third columns of pop map, which are sample ID and population ID.
		# Finally, use for loops to make a population map for each species, and check that each sample
		# is actually present in the VCF.
		cut -d',' -f2 ${pop_map} | sort | uniq > species.txt && \
		cut -f 1,3 -d',' ${pop_map} | tr ', ' '\t' > tmp_pop_map1.txt && \
		for i in `cat species.txt`;
		do
			bcftools query -l \${i}_*_filtered.vcf.gz > \${i}_samples.txt && \
			grep -w \$i tmp_pop_map1.txt > tmp_pop_map2.txt && \
			touch \${i}_pop_map.txt && \
			for j in `cat \${i}_samples.txt`;
			do
				if grep -qw \$j tmp_pop_map2.txt; then
					grep -w \$j tmp_pop_map2.txt >> \${i}_pop_map.txt
				fi
			done
		done
		"""
	}
}


process SELECT_SFS_PROJECTION {
	
	tag "${prep} ${species}"
	
	input:
	tuple path(snps), val(species), val(prep), val(sample_size)
	each path(pop_maps)
	
	output:
	env projection
	
	script:
	"""
	
	/usr/local/bin/easysfs/easySFS.py -avf \
	-i ${snps} \
	-p ${species}_pop_map.txt \
	--preview > projection_options.txt && \
	select_optimal_easySFS_projection.R && \
	projection=`cat projection.txt`
	
	"""
}


process SFS_ESTIMATION {
	
	tag "${prep} ${species}"
	
	publishDir params.sfs_results, mode: 'copy', overwrite: true
	
	input:
	tuple path(snps), val(species), val(prep), val(sample_size)
	each path(pop_maps)
	val projection
	
	output:
	tuple path("*.sfs"), val(species), val(prep), val(sample_size), emit: sfs
	
	script:
	"""
	
	/usr/local/bin/easysfs/easySFS.py -avf \
	-i ${snps} \
	-p ${species}_pop_map.txt \
	-o . \
	--proj ${projection}
	
	mv dadi/${species}.sfs ./${species}_${params.date}.sfs
	
	"""
}


process VISUALIZE_SFS {
	
	tag "${prep} ${species}"
	
	publishDir params.sfs_results, mode: 'copy', overwrite: true
	
	input:
	tuple path(sfs), val(species), val(prep), val(sample_size)
	
	output:
	path "*.pdf"
	
	script:
	"""
	SFS_plotting.R ${sfs} ${species} ${params.date}
	"""
}

process BUILD_STAIRWAY_PLOT_BLUEPRINT {
	
	tag "${prep} ${species}"
	
	publishDir params.stairway_plot_scripts, mode: 'copy', overwrite: true
	
	input:
	tuple path(sfs), val(species), val(prep), val(sample_size)
	
	output:
	path "*.blueprint"
	
	script:
	"""
	
	create_stairwayplot_blueprint.R ${sfs} \
	${species} \
	${species} \
	${prep} \
	${sample_size} \
	${params.genome_length} \
	${params.year_per_generation} \
	${params.mutation_rate} \
	${params.random_seed} \
	${params.whether_folded} \
	/usr/local/bin/stairway_plot_v2.1.1/stairway_plot_v2.1.1/files/stairway_plot_es
	
	"""
	
}

process BUILD_STAIRWAY_PLOT_SCRIPT {
	
	publishDir params.stairway_plot_scripts, mode: 'copy', overwrite: true
	
	input:
	path blueprint
	
	output:
	path "*.sh"
	
	script:
	"""
	java -cp /usr/local/bin/stairway_plot_v2.1.1/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder ${blueprint}
	"""
	
}

process STAIRWAY_PLOT {
	
	publishDir params.stairway_plots, pattern: '*.final.summary.pdf', mode: 'copy', overwrite: true
	publishDir params.stairway_plots, pattern: '*.final.summary.png', mode: 'copy', overwrite: true
	publishDir params.stairway_plots, pattern: '*.final.summary', overwrite: true
	
	input:
	path blueprint_script
	
	output:
	path "*"
	
	when:
	params.stairwayplot == true
	
	script:
	"""
	bash ${blueprint_script}
	"""
	
}

process POP_STRUCTURE_PCA {
	
	publishDir params.params.pca_results, mode: 'copy'
	
	input:
	tuple path(snps), val(species), val(prep)
	
	output:
	path "*.pdf"
	
	when:
	params.principal_component_analysis == true
	
	script:
	"""
	
	plink --vcf ${snps} --double-id --allow-extra-chr \
	--set-missing-var-ids @:#\$1,\$2 \
	--indep-pairwise 50 10 0.1 --out ${species}_${prep}_${params.date}
	
	plink --vcf ${snps} --double-id --allow-extra-chr --set-missing-var-ids @:#\$1,\$2 \
	--extract ${species}_${prep}_${params.date}.prune.in \
	--make-bed --pca --out ${species}_${prep}_${params.date}
	
	plot_PCA.R
	
	"""
	
}

process ADMIXTURE_PLOT {
	
	input:
	tuple path(snps), val(species), val(prep)
	
	output:
	path "*."
	
	when:
	params.admixture_plot == true
	
	script:
	output_prefix = species "_" + prep + "_admix"
	"""
	# Create a temporary plink file for ADMIXTURE
	plink --vcf ${snps} --recode --allow-extra-chr --out ${output_prefix}
	
	# Run ADMIXTURE on the plink file
	admixture --cv ${output_prefix}.ped \
	${params.admixture_plot_K} | tee ${output_prefix}.log
	
	# Plot ADMIXTURE results
	plot_admixture.R ${output_prefix} ${params.admixture_plot_K}
	"""
}


// --------------------------------------------------------------- //
