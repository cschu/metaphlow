#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { fastq_input } from "./nevermore/workflows/input"
include { run_metaphlan4; combine_metaphlan4; collate_metaphlan4_tables } from "./nevermore/modules/profilers/metaphlan4"
include { run_metaphlan3; combine_metaphlan3; collate_metaphlan3_tables } from "./nevermore/modules/profilers/metaphlan3"
include { humann3 } from "./metaphlow/workflows/humann3"
include { samestr_full; samestr_post_convert } from "./metaphlow/workflows/samestr"


def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir

if (!params.fastq_input_pattern) {
	params.fastq_input_pattern = "**[._]{fastq.gz,fq.gz,fastq.bz2,fq.bz2}"
}
def fastq_input_pattern = input_dir + "/" + params.fastq_input_pattern

params.skip_alignment = true


workflow {

	if (!params.skip_profiling) {

		fastq_input(
			Channel.fromPath(input_dir + "/**"),
			Channel.of(null)
		)

		fastq_input_ch = fastq_input.out.fastqs

		if (params.ignore_samples) {
			ignore_samples = params.ignore_samples.split(",")
			print ignore_samples
			fastq_input_ch = fastq_input_ch
				.filter { !ignore_samples.contains(it[0].id) }
		}
		
		fastq_input_ch.dump(pretty: true, tag: "fastq_input_ch")
		nevermore_main(fastq_input_ch)

		fastq_ch = nevermore_main.out.fastqs	

		fastq_ch = fastq_ch
			.map { sample, fastqs ->
				sample_id = sample.id.replaceAll(/\.singles$/, "")
				return tuple(sample_id, fastqs)
			}
			.groupTuple()
			.map { sample_id, fastqs ->
				def meta = [:]
				meta.id = sample_id				
				return tuple(meta, [fastqs].flatten())
			}

		run_metaphlan4(fastq_ch, params.mp4_db)
		
		if (params.mp4_collate || params.run_humann3) {
			collate_metaphlan4_tables(
				run_metaphlan4.out.mp4_table
					.map { sample, table -> return table }
					.collect()
			)
		}

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

		if (params.run_samestr) {
			samestr_full(
				run_metaphlan4.out.mp4_sam,
				run_metaphlan4.out.mp4_table
			)
		}	
      
    } else if (params.run_samestr) {

		ss_converted = Channel.fromPath(input_dir + "/**.npz")
			.map { file ->
					def species = file.name.replaceAll(/[.].*/, "")
					return tuple(species, file)
			}
			.groupTuple(sort: true)

		mp4_tables = Channel.fromPath(input_dir + "/**.mp4.txt")
			.map { file ->
				def meta = [:]
				meta.id = file.name.replaceAll(/\.txt$/, "")
				return tuple(meta, file)
			}

		samestr_post_convert(ss_converted, mp4_tables)        

	}

}


