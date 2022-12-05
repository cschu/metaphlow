process run_metaphlan3 {
	
	input:
	tuple val(sample), path(fastqs)
	path(mp3_db)

	output:
	tuple val(sample), path("metaphlan3_tables/${sample.id}.mp3.txt"), emit: mp3_table
	// tuple val(sample), path("${sample.id}.bowtie2.bz2"), emit: mp3_bt2
	
	script:
	def mp3_params = "--index \$(cat \$(readlink ${mp3_db})/mpa_latest) --bowtie2db \$(readlink ${mp3_db}) --input_type fastq --nproc ${task.cpus} --tmp_dir tmp/"
	def mp3_input = ""
	def bt2_out = "--bowtie2out ${sample.id}.bowtie2.bz2"

	
	if (fastqs instanceof Collection && fastqs.size() == 2) {
		mp3_input = "${sample.id}_R1.fastq.gz,${sample.id}_R2.fastq.gz"
	} else if (fastqs instanceof Collection && fastqs.size() == 3) {
		mp3_input = "${sample.id}_R1.fastq.gz,${sample.id}_R2.fastq.gz,${sample.id}.singles_R1.fastq.gz"
	} else {
		mp3_input = "${fastqs}"
	}

	"""
	mkdir -p tmp/ metaphlan3_tables/

	metaphlan ${mp3_input} ${mp3_params} ${bt2_out} -o metaphlan3_tables/${sample.id}.mp3.txt
	"""
}

process combine_metaphlan3 {

	input:
	tuple val(sample), path(bt2)

	output:
	tuple val(sample), path("metaphlan3/${sample.id}/${sample.id}.mp3.txt"), emit: mp3_table

	script:
	def mp3_params = "--input_type bowtie2out --nproc ${task.cpus} --tmp_dir tmp/"
	def bt2_out = "--bowtie2out ${sample.id}.bowtie2.bz2"
	def mp3_input = "${sample.id}.bowtie2.bz2,${sample.id}.singles.bowtie2.bz2"
	"""
	mkdir -p metaphlan3/${sample.id}/

	metaphlan ${mp3_input} ${mp3_params} -o metaphlan3/${sample.id}/${sample.id}.mp3.txt
	"""
}


process collate_metaphlan3_tables {

	input:
	path(tables)

	output:
	path("metaphlan3_abundance_table.txt")

	script:
	"""
	mkdir -p metaphlan3/

	merge_metaphlan_tables.py ${tables} > metaphlan3_abundance_table.txt
	"""

}