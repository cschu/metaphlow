process run_metaphlan4 {
	
	input:
	tuple val(sample), path(fastq)
	path(mp4_db)

	output:
	// tuple val(sample), path("mp4/${sample.id}/${sample.id}.mp4.txt"), emit: mp4_table
	tuple val(sample), path("${sample.id}.bowtie2.bz2"), emit: mp4_bt2
	
	script:
	def mp4_params = "--bowtie2db ${mp4_db} --input_type fastq --nproc ${task.cpus} --tmp_dir tmp/"
	def mp4_input = ""
	def bt2_out = "--bowtie2out ${sample.id}.bowtie2.bz2"
	if (!sample.is_paired) {
		mp4_input = "${fastq}"
	} else {
		mp4_input = "${sample.id}_R1.fastq.gz,${sample.id}_R2.fastq.gz"
	}

	"""
	mkdir -p tmp/

	metaphlan ${mp4_input} ${mp4_params} ${bt2_out} -o ${sample.id}.mp4.txt
	"""
}

process combine_metaphlan4 {

	input:
	tuple val(sample), path(bt2)

	output:
	tuple val(sample), path("metaphlan4/${sample.id}/${sample.id}.mp4.txt"), emit: mp4_table

	script:
	def mp4_params = "--input_type bowtie2out --nproc ${task.cpus} --tmp_dir tmp/"
	def bt2_out = "--bowtie2out ${sample.id}.bowtie2.bz2"
	"""
	mkdir -p mp4/${sample.id}/

	cat ${bt2} | bzip2 -dc - | bzip2 -c - > tmp/bowtie2out.bz2

	cat ${bt2} | metaphlan tmp/bowtie2out.bz2 ${mp4_params} -o metaphlan4/${sample.id}/${sample.id}.mp4.txt
	rm -rvf tmp/
	"""
}


process collate_metaphlan4_tables {

	input:
	path(tables)

	output:
	path("metaphlan4_abundance_table.txt")

	script:
	"""
	mkdir -p metaphlan4/

	merge_metaphlan_tables.py ${tables} > metaphlan4_abundance_table.txt
	"""

}