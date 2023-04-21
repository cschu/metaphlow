#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { nevermore_prep_align } from "./nevermore/workflows/align"
include { fastq_input } from "./nevermore/workflows/input"
include { run_metaphlan4; combine_metaphlan4; collate_metaphlan4_tables } from "./nevermore/modules/profilers/metaphlan4"
include { run_metaphlan3; combine_metaphlan3; collate_metaphlan3_tables } from "./nevermore/modules/profilers/metaphlan3"
include { reduce_metaphlan_profiles; generate_humann_joint_index; run_humann3 } from "./nevermore/modules/profilers/humann3"


def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir

if (!params.fastq_input_pattern) {
	params.fastq_input_pattern = "**[._]{fastq.gz,fq.gz,fastq.bz2,fq.bz2}"
}
def fastq_input_pattern = input_dir + "/" + params.fastq_input_pattern

params.skip_alignment = true

workflow {

	fastq_input(
		// Channel.fromPath(fastq_input_pattern)
		Channel.fromPath(input_dir + "/*", type: "dir")
	)
	
	fastq_ch = fastq_input.out.fastqs

	nevermore_main(fastq_ch)

	fastq_ch = nevermore_main.out.fastqs
	// fastq_ch.view()


	if (!params.skip_profiling) {

		fastq_ch = fastq_ch
			.map { sample, fastqs ->
				sample_id = sample.id.replaceAll(/\.singles$/, "")
				return tuple(sample_id, fastqs)
			}
			.groupTuple()
			// .map { sample_id, fastqs -> return tuple(sample_id, fastqs.flatten()) }
			// .groupTuple(sort: true)
			.map { sample_id, fastqs ->
				def meta = [:]
				meta.id = sample_id				
				return tuple(meta, fastqs.flatten())
			}

		// fastq_ch.view()

		run_metaphlan4(fastq_ch, params.mp4_db)
		
		// mp4_ch = run_metaphlan4.out.mp4_bt2
		// 	.map { sample, bt2 ->
		// 		sample_id = sample.id.replaceAll(/\.singles$/, "")
		// 		return tuple(sample_id, bt2)
		// 	}
		// 	.groupTuple(sort: true)
		// 	.map { sample_id, bt2 ->
		// 		def meta = [:]
		// 		meta.id = sample_id				
		// 		return tuple(meta, bt2)
		// 	}

		// combine_metaphlan4(mp4_ch)
		mp4_tables_ch = run_metaphlan4.out.mp4_table
			.map { sample, table -> return table }

		collate_metaphlan4_tables(mp4_tables_ch.collect())


		if (params.run_humann3) {
			reduce_metaphlan_profiles(
				collate_metaphlan4_tables.out.mp4_abundance_table,
				"max"
			)

			generate_humann_joint_index(
				reduce_metaphlan_profiles.out.mp_reduced_profiles,
				params.humann_nuc_db
			)

			humann_input_ch = fastq_ch
				.join(run_metaphlan4.out.mp4_table, remainder: false)

			run_humann3(
				humann_input_ch,
				generate_humann_joint_index.out.joint_bt2_index,
				generate_humann_joint_index.out.chocophlan_db
			)

		}

		
		if (params.run_metaphlan3) {
			run_metaphlan3(fastq_ch, params.mp3_db)
			mp3_tables_ch = run_metaphlan3.out.mp3_table
			.map { sample, table -> return table }

			collate_metaphlan3_tables(mp3_tables_ch.collect())
		}

	}

}
