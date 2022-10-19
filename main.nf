#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"
include { db_filter; db2bed3 } from "./nevermore/modules/align/helpers"

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir

if (!params.fastq_input_pattern) {
	params.fastq_input_pattern = "**[._]{fastq.gz,fq.gz}"
}
def fastq_input_pattern = input_dir + "/" + params.fastq_input_pattern


workflow {

	fastq_input(
		Channel.fromPath(fastq_input_pattern)
	)
	
	fastq_ch = fastq_input.out.fastqs

	nevermore_main(fastq_ch)

	if (!params.skip_dbfilter) {
		db2bed3(params.gq_db)
		db_filter(nevermore_main.out.alignments, db2bed3.out.db)
	}

	if (!params.skip_profiling) {

		

		//gffquant_flow(nevermore_main.out.alignments)
		gffquant_flow(db_filter.out.bam)

	}

}
