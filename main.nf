#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { fastq_input } from "./nevermore/workflows/input"
include { run_metaphlan4; combine_metaphlan4; collate_metaphlan4_tables } from "./nevermore/modules/profilers/metaphlan4"
include { run_metaphlan3; combine_metaphlan3; collate_metaphlan3_tables } from "./nevermore/modules/profilers/metaphlan3"
include { run_motus } from "./nevermore/modules/profilers/motus"
include { humann3 } from "./metaphlow/workflows/humann3"
include { samestr_full; samestr_post_convert; samestr_post_merge } from "./metaphlow/workflows/samestr"

// include { ss_load_convert_tarball } from "./metaphlow/modules/samestr/samestr_load"


def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir

if (!params.fastq_input_pattern) {
	params.fastq_input_pattern = "**[._]{fastq.gz,fq.gz,fastq.bz2,fq.bz2}"
}
def fastq_input_pattern = input_dir + "/" + params.fastq_input_pattern


workflow metaphlow_upstream {
	take:
		input_fastq_ch

	main:
		fastq_input(
			// Channel.fromPath(fastq_input_pattern),
			input_fastq_ch,
			Channel.of(null)
		)

		fastq_input_ch = fastq_input.out.fastqs

		fastq_input_ch.dump(pretty: true, tag: "fastq_input_ch")
		nevermore_main(fastq_input_ch)

		fastq_ch = nevermore_main.out.fastqs

		fastq_ch = fastq_ch
			.map { sample, fastqs ->
				sample_id = sample.id.replaceAll(/\.singles$/, "")
				return tuple(sample_id, fastqs)
			}
			.groupTuple(size: ((params.single_end_libraries) ? 1 : 2), remainder: true)
			.map { sample_id, fastqs ->
				def meta = [:]
				meta.id = sample_id				
				return tuple(meta, [fastqs].flatten())
			}

		alignments_ch = Channel.empty()
		tax_profiles_ch = Channel.empty()

		if (params.sstr_profiler == "motus") {

			run_motus(fastq_ch, params.motus_db)

			alignments_ch = run_motus.out.motus_bam
			tax_profiles_ch = run_motus.out.motus_profile

		} else {

			run_metaphlan4(fastq_ch, params.mp4_db)

			collate_metaphlan4_tables(
				run_metaphlan4.out.mp4_table
					.map { sample, table -> return table }
					.collect()
			)

			if (params.run_humann3) {
				humann3(
					collate_metaphlan4_tables.out.mp4_abundance_table,
					run_metaphlan4.out.mp4_table,
					fastq_ch
				)
			}

			if (params.run_metaphlan3) {
				run_metaphlan3(fastq_ch, params.mp3_db)
				mp3_tables_ch = run_metaphlan3.out.mp3_table
					.map { sample, table -> return table }

				collate_metaphlan3_tables(mp3_tables_ch.collect())
			}

			alignments_ch = run_metaphlan4.out.mp4_sam
			tax_profiles_ch = run_metaphlan4.out.mp4_table
			
		}
	
	emit:
		alignments = alignments_ch
		tax_profiles = tax_profiles_ch

}


workflow {

	if (params.run_mode == "full") {

		metaphlow_upstream(Channel.fromPath(fastq_input_pattern))		

		if (params.run_samestr) {
			samestr_full(metaphlow_upstream.out.alignments, metaphlow_upstream.out.tax_profiles)
		}
      
	} else {

		// mp4_tables = Channel.fromPath(input_dir + "/mp4_profiles/**.mp4.txt")
		mp4_tables = Channel.fromPath(input_dir + "/**.mp4.txt")
			.map { file ->
				def meta = [:]
				meta.id = file.name.replaceAll(/\.txt$/, "")
				return [ meta, file ]
			}
		
		if (params.run_mode == "samestr_convert") {

			convert_input = Channel.empty()

			if (params.load_convert_tarball) {

				convert_input = Channel.fromPath(params.load_convert_tarball)

			} else {

				mp4_alignments = Channel.fromPath(input_dir + "/**.sam.bz2")
					.map { file ->
						def meta = [:]
						meta.id = file.name.replaceAll(/\.sam\.bz2$/, "")
						return [ meta, file ]
					}

				convert_input = mp4_alignments

			}

			samestr_full(convert_input, mp4_tables)

		} else {

			npz_ch = Channel.fromPath(input_dir + "/**.npz")
			
			if (params.run_mode == "samestr_post_convert") {

				if (params.load_merge_tarball) {

					merge_input_ch = Channel.fromPath(params.load_merge_tarball)

				} else {

					// ss_converted = Channel.fromPath(input_dir + "/**.npz")
					merge_input_ch = npz_ch
						.map { file ->
							def species = file.name.replaceAll(/[.].*/, "")
							return [ species, file ]
						}
						.groupTuple(sort: true)

				}

				samestr_post_convert(npz_ch, mp4_tables)        

			} else if (params.run_mode == "samestr_post_merge") {

				if (params.load_filter_tarball) {

					post_merge_input_ch = Channel.fromPath(params.load_filter_tarball)

				} else {			

					// npz_ch = Channel.fromPath(input_dir + "/**.npz")
					npz_ch = npz_ch
						.map { file -> [ file.name.replaceAll(/.npz$/, ""), file ] }
					names_ch = Channel.fromPath(input_dir + "/**.names.txt")
						.map { file -> [ file.name.replaceAll(/.names.txt$/, ""), file ] }

					post_merge_input_ch = npz_ch.join(names_ch, by: 0)

				}

				samestr_post_merge(post_merge_input_ch, mp4_tables)

			}

		}

	}	

}
