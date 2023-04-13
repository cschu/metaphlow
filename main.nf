#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { nevermore_prep_align } from "./nevermore/workflows/align"
include { fastq_input } from "./nevermore/workflows/input"
include { run_metaphlan4; combine_metaphlan4; collate_metaphlan4_tables } from "./nevermore/modules/profilers/metaphlan4"
include { run_metaphlan3; combine_metaphlan3; collate_metaphlan3_tables } from "./nevermore/modules/profilers/metaphlan3"

include { samestr } from "./nevermore/workflows/samestr"
include { run_samestr_convert; run_samestr_merge; run_samestr_filter; run_samestr_stats; run_samestr_compare; run_samestr_summarize } from "./nevermore/modules/profilers/samestr"


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

		if (params.run_samestr) {
			samestr_input_ch = run_metaphlan4.out.mp4_sam
				.join(run_metaphlan4.out.mp4_table)
				.map { sample, sam, profile -> return tuple(sam, profile) }
				.collect()

			run_samestr_convert(
				samestr_input_ch,
				params.samestr_marker_db
			)

			grouped_npy_ch = run_samestr_convert.out.sstr_npy
				.flatten()
				.map { file ->
						def species = file.name.replaceAll(/[.].*/, "")
						return tuple(species, file)
				}
				.groupTuple(sort: true)
            
			run_samestr_merge(grouped_npy_ch)
			run_samestr_filter(
			 	run_samestr_merge.out.sstr_npy,
			 	params.samestr_marker_db
			)
			run_samestr_stats(run_samestr_filter.out.sstr_npy)
			run_samestr_compare(run_samestr_filter.out.sstr_npy)

			// symlink all sstr_compare/mp_profiles
			run_samestr_summarize(
				run_samestr_compare.out.sstr_compare.collect(),
				samestr_input_ch
					.map { sam, profile -> return profile }
					.collect()
			)

			// samestr(
			// 	samestr_input_ch
			// 	// .join(run_metaphlan4.out.mp4_table)
			// 	// .map { sample, sam, profile ->
			// 	// 	return tuple(sam, profile)
			// 	// }
			// 	// .collect()
			// )
		}

		
		if (params.run_metaphlan3) {
			run_metaphlan3(fastq_ch, params.mp3_db)
			mp3_tables_ch = run_metaphlan3.out.mp3_table
			.map { sample, table -> return table }

			collate_metaphlan3_tables(mp3_tables_ch.collect())
		}

	}

}
