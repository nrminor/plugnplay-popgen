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
	
	ch_reference = Channel
		.fromPath ( params.reference )

	// Workflow steps
	SNP_FILTERING (
		ch_vcfs
	)
	
	// SINGLE_VCF_STATS (
	// 	SNP_FILTERING.out.vcf
	// )
	
	MERGE_VCFS (
		SNP_FILTERING.out.vcf
			.map { sample_id, species, population, prep_type, platform, raw_vcf -> tuple( file(raw_vcf), species, prep_type ) }
			.filter { file(it[0]).countLines() > 0 }
			.groupTuple( by: [1,2] ),
		SNP_FILTERING.out.index.collect()
	)
	
	// MULTI_VCF_STATS (
	// 	MERGE_VCFS.out.vcf
	// )
	
	MAKE_POP_MAP (
		ch_vcfs
			.map { sample_id, species, population, prep_type, platform, raw_vcf -> "${sample_id}, ${species}, ${population}" }
			.collect()
			.collectFile( name: "pop_map.txt", newLine: true )
	)
	
	SFS_ESTIMATION (
		MERGE_VCFS.out.vcf,
		MAKE_POP_MAP.out
	)
	
	VISUALIZE_SFS (
		SFS_ESTIMATION.out.sfs
	)
	
	BUILD_STAIRWAY_PLOT_SCRIPT (
		SFS_ESTIMATION.out.sfs
	)
	
	// STAIRWAY_PLOT ( 
	// 	BUILD_STAIRWAY_PLOT_SCRIPT.out
	// )
	
	// POP_STRUCTURE_PCA ( )
	// 
	// ADMIXTURE_PLOT ( )

}
// --------------------------------------------------------------- //




// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //

// working with single nucleotide polymorphisms
params.snp_results = params.results + "/" + "per_sample_SNPs"
params.single_sample_snp_stats = params.snp_results + "/" + "snp_stats"
params.merged_snps = params.results + "/" + "per_species_SNPs"
params.multisample_snp_stats = params.merged_snps + "/" + "snp_stats"

// Additional, optional analyses
params.analyses = params.results + "/" + "analyses"
params.sfs_results = params.analyses + "/" + "site_frequency_spectra"
params.angsd_results = params.analyses + "/" + "ANGSD_outputs"
params.angsd_sfs_plots = params.angsd_results + "/" + "SFS_plots"
params.stairway_plots = params.analyses + "/" + "Stairway_plots"
params.stairway_plot_run_date = params.analyses + "/" + "Stairway_plots" + "/" + params.date
params.admixture_results = params.analyses + "/" + "admixture_plots"
params.pca_results = params.analyses + "/" + "PCA_plots"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //


process SNP_FILTERING {
	
	tag "${sample}"
	publishDir params.snp_results, pattern: "*.vcf.gz", mode: 'copy'
	publishDir params.snp_results, pattern: "*.tbi", mode: 'copy'
	
	cpus 2
	
	input:
	tuple val(sample), val(species), val(pop), val(prep), val(platform), path(vcf)
	
	output:
	tuple val(sample), val(species), val(pop), val(prep), val(platform), path("*_filtered.vcf.gz"), emit: vcf
	path "*.tbi", emit: index
	
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
	
	input:
	tuple val(sample), val(species), val(pop), val(prep), val(platform), path(filtered_vcf)
	
	output:
	path "*.stats"
	
	script:
	"""
	
	bgzip -d -c ${filtered_vcf} > ${sample}_${prep}_filtered.vcf && \
	bcftools stats \
	-F ${params.reference} \
	-s - ${filtered_vcf} > ${filtered_vcf}.stats
	
	"""
	
}


process MERGE_VCFS {
	
	tag "${prep} ${species}"
	
	publishDir params.merged_snps, mode: 'copy'
	
	cpus 8
	
	input:
	tuple path(vcf_files), val(species), val(prep)
	path index_files
	
	output:
	tuple path("*.vcf.gz"), val(species), val(prep), env(sample_size), emit: vcf
	path "*.txt"
	
	shell:
	'''
	
	bcftools merge \
	--merge snps \
	--output-type z \
	--threads !{task.cpus} \
	--output !{species}_multisample_!{prep}.vcf.gz \
	*!{species}*.vcf.gz
	
	sample_ids=`bcftools query -l !{species}_multisample_!{prep}.vcf.gz`
	sample_size=`bcftools query -l !{species}_multisample_!{prep}.vcf.gz | wc -l`
	minInd="$((${sample_size} - !{params.max_snp_missingness}))"
	
	touch !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "VCF FILTER SETTINGS APPLIED TO !{prep} !{species} SNPS" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "----------------------------------------------------------" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minor allele frequency: !{params.minor_allele_frequency}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum samples: ${minInd} / ${sample_size}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum variant quality score: !{params.min_quality}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum Depth: !{params.min_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Maximum Depth: !{params.max_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Sample IDs included: ${sample_ids}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	
	'''
	
}


process MULTI_VCF_STATS {
	
	publishDir params.multisample_snp_stats, mode: 'copy'
	
	input:
	tuple path(vcf), val(species), val(prep), val(sample_size)
	
	output:
	path "*.stats"
	
	shell:
	'''
	
	sample_ids=`bcftools query -l !{vcf}`
	minInd="$((${sample_size} - !{params.max_snp_missingness}))"
	
	bcftools stats \
	-F ${params.reference} \
	-s - ${vcf} > ${vcf}.stats
	
	touch !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "VCF FILTER SETTINGS APPLIED TO !{prep} !{species} SNPS" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "----------------------------------------------------------" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minor allele frequency: !{params.minor_allele_frequency}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum samples: ${minInd} / !{sample_size}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum variant quality score: !{params.min_quality}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum Depth: !{params.min_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Maximum Depth: !{params.max_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Sample IDs included: ${sample_ids}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	
	'''
	
}


process MAKE_POP_MAP {
	
	input:
	path pop_map
	
	output:
	path "*.txt"
	
	script:
	if ( params.whole_species_mode == true ){
		"""
		# take first two columns of pop map
		cut -d',' -f2 ${pop_map} | sort | uniq > species.txt
		cut -f 1,2 -d',' ${pop_map} > tmp_pop_map.txt
		for i in `cat species.txt`;
		do
			grep -w \$i tmp_pop_map.txt > \${i}_pop_map.txt
		done
		rm tmp_pop_map.txt
		"""
	} else {
		"""
		# take first and third columns of pop map
		cut -d',' -f2 ${pop_map} | sort | uniq > species.txt
		cut -f 1,2 -d',' ${pop_map} > tmp_pop_map.txt
		for i in `cat species.txt`;
		do
			grep -w \$i tmp_pop_map.txt > \${i}_pop_map.txt
		done
		rm tmp_pop_map.txt
		"""
	}
}


process SFS_ESTIMATION {
	
	tag "${prep} ${species}"
	
	publishDir params.sfs_results, mode: 'copy', overwrite: true
	
	input:
	tuple path(snps), val(species), val(prep), val(sample_size)
	each path(pop_maps)
	
	output:
	tuple path("*.sfs"), val(species), val(prep), val(sample_size), emit: sfs
	path "*_vcf_filter_settings_*.txt", emit: vcf_filters
	
	shell:
	'''
	
	sample_ids=`bcftools query -l !{snps}`
	minInd="$((!{sample_size} - !{params.max_snp_missingness}))"
	
	/usr/local/bin/easysfs/easySFS.py -avf \
	-i !{snps} \
	-p !{species}_pop_map.txt \
	-o . \
	-f \
	--proj !{sample_size}
	
	mv dadi/!{species}.sfs ./!{species}_!{params.date}.sfs
	
	touch !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "VCF FILTER SETTINGS APPLIED TO !{prep} !{species} SNPS" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "----------------------------------------------------------" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minor allele frequency: !{params.minor_allele_frequency}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum samples: ${minInd} / !{sample_size}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum variant quality score: !{params.min_quality}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Minimum Depth: !{params.min_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Maximum Depth: !{params.max_depth}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	echo "Sample IDs included: ${sample_ids}" >> !{prep}_!{species}_vcf_filter_settings_!{params.date}.txt
	
	'''
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

process BUILD_STAIRWAY_PLOT_SCRIPT {
	
	tag "${prep} ${species}"
	
	input:
	tuple path(sfs), val(species), val(prep), val(sample_size)
	
	output:
	path "*.sh"
	
	when:
	params.stairwayplot == true
	
	script:
	"""
	
	create_stairwayplot_script.R ${sfs} \
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
	
	java -cp /usr/local/bin/stairway_plot_v2.1.1/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder ${species}_${species}_${prep}.blueprint
	
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
