#!/usr/bin/env nextflow

/*
Nextflow script for main NARST screening workflow
*/

/*
FUNCTIONS
*/
import nextflow.splitter.CsvSplitter

def fetchRunAccessions( tsv ) {

    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )

    splitter.parseHeader( reader )

    List<String> run_accessions = []
	List<String> run_accessions_filtered = []
    Map<String,String> row

    while( row = splitter.fetchRecord( reader ) ) {

       run_accessions.add( row['SRR'] )
    }
	
	run_accessions_filtered = run_accessions.minus('')
	
    return run_accessions_filtered
}

/*
PROCESSES
*/
process READMETRICS {
	tag "CG-Pipeline on $wgs_id"
	publishDir (
		path: "${params.outdir}/${wgs_id}/readmetrics/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.endsWith('.txt')) "$filename"
					else null
		}
	)
	containerOptions {
		'-B $PWD:$PWD'
		'--workdir $PWD'
	}	
	errorStrategy = 'ignore'
	
	input:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genus), val(species)
	
	output:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), path('readmetrics.txt')
	
	script:
	genome_size = ''
	min_coverage = ''
	
	if (genus == 'Salmonella') {
		genome_size = 5000000
		min_coverage = 30
	} else if (genus == 'Escherichia') {
		genome_size = 5000000
		min_coverage = 40
	} else if (genus == 'Shigella') {
		genome_size = 5000000
		min_coverage = 40
	} else if (genus == 'Campylobacter') {
		genome_size = 1800000
		min_coverage = 20
	} else if (genus == 'Vibrio') {
		genome_size = 5000000
		min_coverage = 40
	} else if (genus == 'Yersinia') {
		genome_size = 4700000
		min_coverage = 40
	} else {
		genome_size = 5000000
		min_coverage = 40
	}
	"""
	run_assembly_readMetrics.pl *.fastq.gz --fast --numcpus 16 -e $genome_size  | sort -k3,3n > readmetrics.txt
	"""
}

process COVERAGECUTOFF {
	tag "Calculate coverage cut-off for $wgs_id"
	errorStrategy = 'ignore'
	
	input:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), path(readmetrics)
	
	output:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), stdout
	
	script:
	"""
	cat ${readmetrics} | awk 'BEGIN { FS = "\\t"} {if (\$1 != "File"){print \$1, \$9/10}}' | sed -e 's/_[1,2].fastq.gz//'| awk '{a[\$1] += \$2} END{for (i in a) print i, a[i]}' | awk '{print \$2}'
	"""
	
}

process ASSEMBLE {
	tag "Shovill on $wgs_id"
	publishDir (
		path: "${params.outdir}/${wgs_id}/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.matches("shovill")) "$filename"
					else if (filename.matches("lowCoverageIsolates.tsv")) "$filename"
					else null
		}
	)
	beforeScript "TMPDIR=''"
	errorStrategy = 'ignore'
	
	input:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), val(covcutoff)
	
	output:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), val(covcutoff), path('shovill/contigs.fa'), path("shovill")
	
	script:
	def value = covcutoff.toFloat()
	def mincov = min_coverage / 10
	
	if (value >= mincov)
		"""
	shovill --outdir "shovill" --R1 ${read1} --R2 ${read2} --mincov ${min_coverage} --trim --namefmt ${wgs_id}_contig%05d --depth 100
		"""
	else if (value < mincov)
		"""
	printf "${wgs_id}\tNARMSWF5.0" > lowCoverageIsolates.tsv
		"""
}

process ARANDPLASMIDSCREEN {
	tag "Staramr on ${wgs_id}"
	publishDir (
		path: "${params.outdir}/${wgs_id}/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.matches("staramr")) "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), val(covcutoff), path(shovill_contigs)
	
	output:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), val(covcutoff), path(shovill_contigs), path("staramr")
	
	script:
	if (genus == 'Salmonella')
		"""
		staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
		staramr search --database databases --pointfinder-organism "salmonella" --exclude-genes-file ${params.staramr_salmonella_exclude} --pid-threshold 90 --percent-length-overlap-resfinder 50  -o "staramr" ${shovill_contigs}
		"""
	else if (genus == 'Escherichia')
		"""
		staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
		staramr search --database databases --exclude-genes-file ${params.staramr_escherichia_exclude} --pid-threshold 90 --percent-length-overlap-resfinder 50 -o "staramr" ${shovill_contigs}
		"""
	else if (genus == 'Shigella')
		"""
		staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
		staramr search --database databases --exclude-genes-file ${params.staramr_escherichia_exclude} --pid-threshold 90 --percent-length-overlap-resfinder 50 -o "staramr" ${shovill_contigs}
		"""
	else if (genus == 'Campylobacter' && species == 'jejuni')
		"""
		staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
		staramr search --database databases --pointfinder-organism "campylobacter" --pid-threshold 90 --percent-length-overlap-resfinder 50 --exclude-genes-file ${params.staramr_campylobacter_exclude} -o "staramr" ${shovill_contigs}
		"""
	else if (genus == 'Campylobacter' && species == 'coli')
		"""
		staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
		staramr search --database databases --pid-threshold 90 --percent-length-overlap-resfinder 50 --exclude-genes-file ${params.staramr_campylobacter_exclude} -o "staramr" ${shovill_contigs}
		"""
	else if (genus == 'Campylobacter' && species == 'upsaliensis')
		"""
		staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
		staramr search --database databases --pid-threshold 90 --percent-length-overlap-resfinder 50 --exclude-genes-file ${params.staramr_campylobacter_exclude} -o "staramr" ${shovill_contigs}
		"""
	else if (genus == 'Vibrio' && species == 'cholerae')
		"""
		staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
		staramr search --database databases --pid-threshold 90 --percent-length-overlap-resfinder 50 --no-exclude-genes -o "staramr" ${shovill_contigs}
		"""
	else if (genus == 'Vibrio' && species == 'parahaemolyticus')
		"""
		staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
		staramr search --database databases --pid-threshold 90 --percent-length-overlap-resfinder 50 --no-exclude-genes -o "staramr" ${shovill_contigs}
		"""	
	else if (genus == 'Yersiniae')
		"""
		staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
		staramr search --database databases --pid-threshold 90 --percent-length-overlap-resfinder 50 --no-exclude-genes -o "staramr" ${shovill_contigs}
		"""
	else
		"""
		staramr db update --resfinder-commit ${params.staramr_resfinder_commit} --pointfinder-commit ${params.staramr_pointfinder_commit} --plasmidfinder-commit ${params.staramr_plasmidfinder_commit} databases
		staramr search --database databases --pid-threshold 90 --percent-length-overlap-resfinder 50 --no-exclude-genes -o "staramr" ${shovill_contigs}
		"""
}

process AMRFINDER {
	tag "AMRFinderPlus on ${wgs_id}"
	publishDir (
		path: "${params.outdir}/${wgs_id}/amrfinderplus/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.matches("mutational.tsv")) "$filename"
					else if(filename.matches("${wgs_id}.amrfinder.fasta")) "$filename"
					else if(filename.matches("${wgs_id}.flank5.fasta")) "$filename"
					else if(filename.matches("amrfinder.tsv")) "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	
	input:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), val(covcutoff), path(shovill_contigs)
	
	output:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), val(covcutoff), path(shovill_contigs), path("mutational.tsv"), path("${wgs_id}.amrfinder.fasta"), path("${wgs_id}.flank5.fasta"), path("amrfinder.tsv")
	
	script:
	if (genus == 'Salmonella')
		"""
		amrfinder --nucleotide ${shovill_contigs} --organism Salmonella --mutation_all "mutational.tsv" --name ${wgs_id} --nucleotide_output "${wgs_id}.amrfinder.fasta" --nucleotide_flank5_output "${wgs_id}.flank5.fasta" --nucleotide_flank5_size 100 --ident_min -1 -o "amrfinder.tsv"
		"""
	else if (genus == 'Escherichia')
		"""
		amrfinder --nucleotide ${shovill_contigs} --organism Escherichia --mutation_all "mutational.tsv" --name ${wgs_id} --nucleotide_output "${wgs_id}.amrfinder.fasta" --nucleotide_flank5_output "${wgs_id}.flank5.fasta" --nucleotide_flank5_size 100 --ident_min -1 -o "amrfinder.tsv"
		"""
	else if (genus == 'Shigella')
		"""
		amrfinder --nucleotide ${shovill_contigs} --organism Escherichia --mutation_all "mutational.tsv" --name ${wgs_id} --nucleotide_output "${wgs_id}.amrfinder.fasta" --nucleotide_flank5_output "${wgs_id}.flank5.fasta" --nucleotide_flank5_size 100 --ident_min -1 -o "amrfinder.tsv"
		"""
	else if (genus == 'Campylobacter')
		"""
		amrfinder --nucleotide ${shovill_contigs} --organism Campylobacter --mutation_all "mutational.tsv" --name ${wgs_id} --nucleotide_output "${wgs_id}.amrfinder.fasta" --nucleotide_flank5_output "${wgs_id}.flank5.fasta" --nucleotide_flank5_size 100 --ident_min -1 -o "amrfinder.tsv"
		"""
	else if (genus == 'Vibrio' && species == 'cholerae')
		"""
		amrfinder --nucleotide ${shovill_contigs} --organism Vibrio_cholerae --mutation_all "mutational.tsv" --name ${wgs_id} --nucleotide_output "${wgs_id}.amrfinder.fasta" --nucleotide_flank5_output "${wgs_id}.flank5.fasta" --nucleotide_flank5_size 100 --ident_min -1 -o "amrfinder.tsv"
		"""
	else
		"""
		amrfinder --nucleotide ${shovill_contigs} --mutation_all "mutational.tsv" --name ${wgs_id} --nucleotide_output "${wgs_id}.amrfinder.fasta" --nucleotide_flank5_output "${wgs_id}.flank5.fasta" --nucleotide_flank5_size 100 --ident_min -1 -o "amrfinder.tsv"
		"""
}

process MUTATIONAL {
	tag "Ariba on ${wgs_id}"
	publishDir (
		path: "${params.outdir}/${wgs_id}/",
		mode: 'copy',
		saveAs: { filename ->
					if (filename.matches("${params.ariba_outdir}")) "$filename"
					else null
		}
	)
	errorStrategy = 'ignore'
	beforeScript "TMPDIR=''"
	
	input:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), val(covcutoff), path(shovill_contigs)
	
	output:
	tuple val(srr_id), path(read1), path(read2), val(wgs_id), val(genome_size), val(min_coverage), val(genus), val(species), val(covcutoff), path(shovill_contigs), path("${params.ariba_outdir}")
	
	script:
	if (genus == 'Salmonella')
		"""
		echo "FROM THE MOMENT I UNDERSTOOD THE WEAKNESS OF MY FLESH, IT DISGUSTED ME"
		mkdir ${params.ariba_outdir}
		"""
	else if (genus == 'Escherichia')
		"""
		ariba run ${params.ariba_ecoli_database} ${read1} ${read2} ${params.ariba_outdir}
		"""
	else if (genus == 'Shigella')
		"""
		ariba run ${params.ariba_ecoli_database} ${read1} ${read2} ${params.ariba_outdir}
		"""
	else if (genus == 'Campylobacter' && species == 'jejuni')
		"""
		mkdir ${params.ariba_outdir}
		ariba run ${params.ariba_cjejuni_database} ${read1} ${read2} ${params.ariba_outdir}/mutational_23S
		"""
	else if (genus == 'Campylobacter' && species == 'coli')
		"""
		mkdir ${params.ariba_outdir}
		ariba run ${params.ariba_ccoli_database_gyrA} ${read1} ${read2} ${params.ariba_outdir}/mutational_gyrA
		ariba run ${params.ariba_ccoli_database_23S} ${read1} ${read2} ${params.ariba_outdir}/mutational_23S
		"""
	else if (genus == 'Campylobacter' && species == 'upsaliensis')
		"""
		mkdir ${params.ariba_outdir}
		ariba run ${params.ariba_cupsaliensis_database_gyrA} ${read1} ${read2} ${params.ariba_outdir}/mutational_gyrA
		ariba run ${params.ariba_cupsaliensis_database_23S} ${read1} ${read2} ${params.ariba_outdir}/mutational_23S
		"""
	else if (genus == 'Vibrio' && species == 'cholerae')
		"""
		ariba run ${params.ariba_vcoli_database} ${read1} ${read2} ${params.ariba_outdir}
		"""
	else if (genus == 'Vibrio' && species == 'parahaemolyticus')
		"""
		ariba run ${params.ariba_vparahaemolyticus_database} ${read1} ${read2} ${params.ariba_outdir}
		"""	
	else if (genus == 'Yersiniae')
		"""
		ariba run ${params.ariba_yersiniae_database} ${read1} ${read2} ${params.ariba_outdir}
		"""
	else
		"""
		echo "EVEN IN DEATH, I SERVE THE OMNISSIAH"
		mkdir ${params.ariba_outdir}
		"""
}

/*
WORKFLOW
*/
workflow {
	if (params.local_reads) {
		Channel
			.fromPath(params.import_sheet, checkIfExists:true)
			.splitCsv(sep: '\t', header:true, strip:true)
			.map{ row -> tuple(row.Read, row.WGS, row.Genus, row.Species) }
			.view()
			.set{import_ch}
		Channel
			.fromFilePairs(params.reads, flat:true, checkIfExists:true)
			.join(import_ch, by: 0)
			.view()
			.set{reads_ch}
	}
	else {
	accessions = fetchRunAccessions(params.import_sheet)
	//println accessions
	
	Channel
		.fromPath(params.import_sheet, checkIfExists:true)
		.splitCsv(sep: '\t', header:true, strip:true)
		.map{ row -> tuple(row.SRR, row.WGS, row.Genus, row.Species) }
		//.view()
		.set{import_ch}

	Channel
		.fromSRA(accessions, apiKey:params.ncbi_api_key, cache: false)
		//.view()
		.set{sra_ch}
	
	sra_ch
		.map { tuple ->
		def srr = tuple[0]
		def read1 = tuple[1][0]
		def read2 = tuple[1][1]
		
		[srr, read1, read2]
		}
		.join(import_ch, by: 0)
		//.view()
		.set{reads_ch}
	}

	readmetrics_ch = READMETRICS(reads_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4,5,6,7,8] }
		//.view()
	coveragecutoff_ch = COVERAGECUTOFF(readmetrics_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4,5,6,7,8] }
		//.view()
	shovill_ch = ASSEMBLE(coveragecutoff_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4,5,6,7,8,9] }
		//.view()
	if (params.amrfinder) {
		staramr_ch = AMRFINDER(shovill_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4,5,6,7,8,9] }
		//.view()
	}
	else {
	staramr_ch = ARANDPLASMIDSCREEN(shovill_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4,5,6,7,8,9] }
		//.view()		
	}
	pointmutation_ch = MUTATIONAL(staramr_ch)
		.map{ outTuple -> outTuple[0,1,2,3,4,5,6,7,8,9] }
		//.view()
}