#!/usr/bin/env nextflow

nextflow.enable.dsl = 2




// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	// Input channels
		
	ch_vcfs = Channel
		.fromPath ( params.samplesheet )
		.splitCsv ( header: true )
		.map { row -> tuple( row.sample_id, row.species, row.population, row.prep_type, row.platform, path(raw_vcf) ) }
	
	ch_reference = Channel
		.fromPath ( params.reference )

	// Workflow steps
	SNP_FILTERING (
		ch_vcfs
	)
	
	SINGLE_VCF_STATS (
		SNP_FILTERING.out
	)
	
	MERGE_VCFS (
		SNP_FILTERING.out.vcf
			.filter { it[3] == "WGS" }
			.collect { unique(it[1]) }
			.mix (
				SNP_FILTERING.out.vcf
				.filter { it[3] == "RAD" }
				.collect { unique(it[1]) }
			)
	)
	
	MULTI_VCF_STATS (
		MERGE_VCFS.out
	)
	
	SFS_ESTIMATION (
		MERGE_VCFS.out,
		ch_vcfs
			.map { sample_id, species, population, prep_type, platform, raw_vcf -> sample_id, population }
			.collect { unique(it[1]) }
	)
	
	VISUALIZE_SFS (
		SFS_ESTIMATION.out.sfs
	)
	
	BUILD_STAIRWAY_PLOT_SCRIPT (
		SFS_ESTIMATION.out.sfs,
		ch_vcfs
			.map { sample_id, species, population, prep_type, platform, raw_vcf -> sample_id, population }
			.collect { unique(it[1]) }
			.count()
	)
	
	STAIRWAY_PLOT ( 
		BUILD_STAIRWAY_PLOT_SCRIPT.out
	)
	
	POP_STRUCTURE_PCA ( )
	
	ADMIXTURE_PLOT ( )

}
// --------------------------------------------------------------- //




// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //

// working with single nucleotide polymorphisms
params.snp_results = params.results + "/" + "4_per_sample_SNPs"
params.single_sample_snp_stats = params.snp_results + "/" + "snp_stats"
params.merged_snps = params.results + "/" + "5_per_species_SNPs"
params.multisample_snp_stats = params.merged_snps + "/" + "snp_stats"

// Additional, optional analyses
params.analyses = params.results + "/" + "analyses"
params.sfs_results = params.analyses + "/" + "site_frequency_spectra"
params.angsd_results = params.analyses + "/" + "ANGSD_outputs"
params.angsd_sfs_plots = params.angsd_results + "/" + "SFS_plots"
params.stairway_plots = params.analyses + "/" + "Stairway_plots"
params.stairway_plot_run_date = params.analyses + "/" + "Stairway_plots" + "/" + ${params.date}
params.admixture_results = params.analyses + "/" + "admixture_plots"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //


process SNP_FILTERING {
	
	tag "${sample}"
	publishDir params.snp_results, pattern: "*.vcf.gz", mode: 'copy'
	publishDir params.snp_results, pattern: "*.tbi", mode: 'copy'
	
	container 'docker://quay.io/biocontainers/vcftools:0.1.16--he513fc3_4'
	
	cpus 2
	
	input:
	tuple val(sample), val(pop), val(prep), val(platform), path(vcf)
	
	output:
	tuple val(sample), val(species), val(pop), val(prep), val(platform), path("*_filtered.vcf.gz"), emit: vcf
	path("*.tbi"), emit: index
	
	script:
	"""
	
	vcftools --gzvcf ${vcf} \
	--max-alleles 2 \
	--mac ${params.minor_allele_count} \
	--max-missing-count ${params.max_snp_missingness} \
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
	
	container 'bcftools container'
	
	input:
	path vcf
	
	output:
	path "*.stats"
	
	script:
	"""
	
	bcftools stats \
	-F ${params.reference} \
	-s - ${vcf} > ${vcf}.stats
	
	"""
	
}


process MERGE_VCFS {
	
	publishDir params.merged_snps, mode: 'copy'
	
	container 'docker://quay.io/biocontainers/bcftools:1.12--h3f113a9_0'
	
	cpus 8
	
	input:
	tuple val(sample), val(species), val(pop), val(prep), val(platform), path(vcf_files)
	
	output:
	tuple path("*${species}*.vcf.gz"), val(species), val(pop), val(prep)
	
	script:
	"""
	
	bcftools merge \
	--merge snps \
	--output-type z \
	--threads 8 \
	--output ${species}_multisample_${prep}.vcf.gz
	*${species}*.vcf.gz
	
	"""
	
}


process MULTI_VCF_STATS {
	
	publishDir params.multisample_snp_stats, mode: 'copy'
	
	container 'bcftools container'
	
	input:
	tuple path(vcf), val(species), val(pop), val(prep)
	
	output:
	path "*.stats"
	
	script:
	"""
	
	bcftools stats \
	-F ${params.reference} \
	-s - ${vcf} > ${vcf}.stats
	
	"""
	
}


process SFS_ESTIMATION {
	
	publishDir params.sfs_results, mode: 'copy', overwrite: true
	
	input:
	tuple path(snps), val(species), val(pop), val(prep)
	path pop_map
	
	output:
	tuple path("*.sfs"), val(species), val(pop), val(prep), emit: sfs
	path "vcf_filter_settings_*.txt", emit: vcf_filters
	
	shell:
	'''
	
	sample_size=`bcftools query command to count samples in a merged vcf`
	minInd="$((${sample_size} - !{params.max_snp_missingness}))"
	
	/absolute/path/to/easySFS/in/docker/container/easySFS.py -avf \
	-i !{snps} \
	-p !{pop_map} \
	-o . \
	-f \
	--proj ${sample_size}
	
	mv dadi/!{species}.sfs ./!{species}_!{params.date}.sfs
	
	touch vcf_filter_settings_!{params.date}.txt
	echo "Minor allele frequency: !{params.minor_allele_frequency}" >> vcf_filter_settings_!{params.date}.txt
	echo "Minimum samples: ${minInd}" >> vcf_filter_settings_!{params.date}.txt
	echo "Minimum variant quality score: !{params.min_quality}" >> vcf_filter_settings_!{params.date}.txt
	echo "Minimum Depth: !{params.min_depth}" >> vcf_filter_settings_!{params.date}.txt
	echo "Maximum Depth: !{params.max_depth}" >> vcf_filter_settings_!{params.date}.txt
	
	'''
}


process VISUALIZE_SFS {
	
	publishDir params.sfs_results, mode: 'copy', overwrite: true
	
	input:
	tuple path(sfs), val(species), val(pop), val(prep)
	
	output:
	path "*.pdf"
	
	script:
	"""
	SFS_plotting.R ${sfs} ${params.date}
	"""
}

process BUILD_STAIRWAY_PLOT_SCRIPT {
	
	input:
	tuple path(sfs), val(species), val(pop), val(prep)
	val sample_size
	
	output:
	path "*.sh"
	
	when:
	params.stairwayplot == true
	
	script:
	"""
	
	java -cp path/to/stairway_plot/files/stairway_plot_es Stairbuilder *.blueprint
	
	create_stairwayplot_script.R ${sfs} \
	${species} \
	${pop} \
	${prep} \
	${sample_size} \
	${params.genome_length} \
	${params.year_per_generation} \
	${params.mutation_rate} \
	${params.random_seed} \
	${params.whether_folded} \
	/absolute/docker/container/path/to/stairway_plot/files/stairway_plot_es
	
	"""
	
}

process STAIRWAY_PLOT {
	
	publishDir params.stairway_plot_run_date, pattern: '*.pdf', mode: 'move', overwrite: true
	publishDir params.stairway_plot_run_date, pattern: '*.png', mode: 'move', overwrite: true
	
	input:
	path blueprint_script
	
	output:
	path "*"
	
	script:
	"""
	bash ${blueprint_script}
	"""
	
}

process POP_STRUCTURE_PCA {}

process ADMIXTURE_PLOT {}






// --------------------------------------------------------------- //
