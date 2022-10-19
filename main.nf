#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { nevermore_prep_align } from "./nevermore/workflows/align"
include { fastq_input } from "./nevermore/workflows/input"
include { run_metaphlan4 } from "./nevermore/modules/profilers/metaphlan4"

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir

if (!params.fastq_input_pattern) {
	params.fastq_input_pattern = "**[._]{fastq.gz,fq.gz}"
}
def fastq_input_pattern = input_dir + "/" + params.fastq_input_pattern

params.skip_alignment = true

workflow {

	fastq_input(
		Channel.fromPath(fastq_input_pattern)
	)
	
	fastq_ch = fastq_input.out.fastqs

	nevermore_main(fastq_ch)

	fastq_ch = nevermore_prep_align.out.fastqs
	fastq_ch.view()


	if (!params.skip_profiling) {

		run_metaphlan4(fastq_ch, params.mp4_db)

	}

}
