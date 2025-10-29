include { run_samestr_convert; run_samestr_merge; run_samestr_filter; run_samestr_stats; run_samestr_compare; run_samestr_summarize; collate_samestr_stats } from "../modules/profilers/samestr"
include { sstr_tarball as sstr_compare_tarball; sstr_tarball as sstr_convert_tarball; sstr_tarball as sstr_filter_tarball; sstr_tarball as sstr_merge_tarball } from "../modules/profilers/samestr"
include { samestr_buffer as sstr_convert_buffer; samestr_buffer as sstr_merge_buffer; samestr_buffer as sstr_filter_buffer } from "../modules/profilers/samestr"

include { ss_load_convert_tarball } from "../modules/profilers/samestr_load"

workflow samestr_post_merge {
	take:
		post_merge_input
		tax_profiles

	main:
		if (params.load_filter_tarball) {

			ss_load_manifest_tarball("sstr_filter", post_merge_input)

			// need to have batch_info here
			sstr_filter_buffer("filter", ss_load_manifest_tarball.out.batch_info, params.filter_batch_size)

			filter_output = sstr_filter_buffer.out.batches
				.splitCsv(header: ['batch_id', 'batch_size', 'file_path'], sep: '\t' )
				.map { item -> [item.batch_id, item.batch_size, item.file_path] }
				.groupTuple(by: 0)

		} else {

			filter_input = post_merge_input
				.splitCsv(header: ['batch_id', 'batch_size', 'file_path'], sep: '\t' )
				.map { item -> [item.batch_id, item.batch_size, item.file_path] }
				.groupTuple(by: 0)

			run_samestr_filter(filter_input, params.samestr_marker_db, params.samestr_sqlite)

			if (!params.skip_filter_tarball) {
				sstr_filter_tarball("sstr_filter", run_samestr_filter.out.sstr_npy.map { batch_id, batch_size, samples, sstr_npy -> [samples, sstr_npy].flatten() }.collect())
			}

			filter_output = run_samestr_filter.out.sstr_npy
				.join(by: 0, run_samestr_filter.out.filter_sentinel)
				.map { batch_id, batch_size, samples, sstr_npy, sentinel -> [ batch_id, batch_size, samples, sstr_npy ] }

		}

		run_samestr_stats(filter_output, params.samestr_marker_db)
		collate_samestr_stats(run_samestr_stats.out.sstr_stats.collect())

		compare_input = filter_output
			.map { batch_id, batch_size, samples, sstr_npy -> [samples, sstr_npy] }
			.flatten()
			.map { file -> [ file.name.replaceAll(/\.(npz|names\.txt)$/, ""), file ] }
			.groupTuple(size: 2, sort: true)
			.map { clade, files -> [ clade, files[1], files[0] ] }

		run_samestr_compare(compare_input, params.samestr_marker_db)
		
		if (!params.skip_compare_tarball) {
			sstr_compare_tarball("sstr_compare", run_samestr_compare.out.sstr_compare.collect())
		}

		run_samestr_summarize(
			run_samestr_compare.out.sstr_compare.collect(),
			tax_profiles.map { sample, table -> return table }.collect(),
			params.samestr_marker_db
		)
}


workflow samestr_post_convert {
	take:
		input_ch
		tax_profiles

	main:
		if (params.load_merge_tarball) {

			ss_load_manifest_tarball("sstr_merge", input_ch)

			merge_info = ss_load_manifest_tarball.out.merge_info

		} else {

			run_samestr_merge(input_ch, params.samestr_marker_db)

			if (!params.skip_merge_tarball) {
				sstr_merge_tarball("sstr_merge", run_samestr_merge.out.sstr_npy.collect())
			}

			merge_info = run_samestr_merge.out.merge_info.collect()

		}

		sstr_merge_buffer("merge", merge_info, params.merge_batch_size)

		merge_output = sstr_merge_buffer.out.batches
			// .splitCsv(header: ['batch_id', 'batch_size', 'file_path'], sep: '\t' )
			// .map { item -> [item.batch_id, item.batch_size, item.file_path] }
			// .groupTuple(by: 0)

		samestr_post_merge(merge_output, tax_profiles)
}


workflow samestr_convert {
	take:
		alignments
		tax_profiles

	main:
		run_samestr_convert(
			alignments.join(tax_profiles),
			params.samestr_marker_db,
			params.samestr_sqlite
		)

		if (!params.skip_convert_tarball) {
			sstr_convert_tarball("sstr_convert", run_samestr_convert.out.sstr_npy.map { sample, data -> [data].flatten() }.collect())
		}

		convert_info = run_samestr_convert.out.convert_info
				.join(run_samestr_convert.out.convert_sentinel, by: 0)
				.map { sample, data, sentinel -> return data }
				.collect()
	
	emit:
		batch_info = convert_info
}



workflow samestr_full {
	take:
		input_ch
		tax_profiles

	main:
		if (params.load_convert_tarball) {

			ss_load_convert_tarball("sstr_convert", input_ch)

			convert_info = ss_load_convert_tarball.out.batch_info

		} else {
			samestr_convert(
				input_ch,
				tax_profiles,
				params.samestr_marker_db,
				params.samestr_sqlite
			)

			convert_info = samestr_convert.batch_info
		}

		if (!params.stop_after_convert) {
			sstr_convert_buffer("convert", convert_info, params.convert_batch_size)

			grouped_npy_ch = sstr_convert_buffer.out.batches
				.splitCsv(header: ['batch_id', 'batch_size', 'file_path'], sep: '\t' )
				.map { item -> [item.batch_id, item.batch_size, item.file_path] }
				.groupTuple(by: 0)

			samestr_post_convert(grouped_npy_ch, tax_profiles)
		}
}
