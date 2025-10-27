include { run_samestr_convert; run_samestr_merge; run_samestr_filter; run_samestr_stats; run_samestr_compare; run_samestr_summarize; collate_samestr_stats } from "../modules/profilers/samestr"
include { sstr_tarball as sstr_compare_tarball; sstr_tarball as sstr_convert_tarball; sstr_tarball as sstr_filter_tarball; sstr_tarball as sstr_merge_tarball } from "../modules/profilers/samestr"
include { samestr_buffer as sstr_convert_buffer; samestr_buffer as sstr_merge_buffer } from "../modules/profilers/samestr"


workflow samestr_post_merge {
	take:
		ss_merged
		tax_profiles
	main:

		run_samestr_filter(ss_merged, params.samestr_marker_db, params.samestr_sqlite)

		if (!params.skip_filter_tarball) {
			sstr_filter_tarball("sstr_filter", run_samestr_filter.out.sstr_npy.map { batch_id, batch_size, samples, sstr_npy -> [samples, sstr_npy].flatten() }.collect())
		}

		filter_output = run_samestr_filter.out.sstr_npy
			.join(by: 0, run_samestr_filter.out.filter_sentinel)
			.map { batch_id, batch_size, samples, sstr_npy, sentinel -> [ batch_id, batch_size, samples, sstr_npy ] }

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
		ss_converted
		tax_profiles
	main:

		run_samestr_merge(ss_converted, params.samestr_marker_db)
		if (!params.skip_merge_tarball) {
			sstr_merge_tarball("sstr_merge", run_samestr_merge.out.sstr_npy.collect())
		}

		merge_info = run_samestr_merge.out.merge_info
			.collect()

		sstr_merge_buffer("merge", merge_info, params.merge_batch_size)

		merge_output = sstr_merge_buffer.out.batches
			.splitCsv(header: ['batch_id', 'batch_size', 'file_path'], sep: '\t' )
			.map { item -> [item.batch_id, item.batch_size, item.file_path] }
			.groupTuple(by: 0)

		samestr_post_merge(merge_output, tax_profiles)
}




workflow samestr_full {

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

		grouped_npy_ch = run_samestr_convert.out.sstr_npy
			.join(run_samestr_convert.out.convert_sentinel, by: 0)
			.map { sample, data, sentinel -> return data }
			.flatten()
			.map { file ->
					def species = file.name.replaceAll(/[.].*/, "")
					return [ species, file ]
			}
			.groupTuple(sort: true)
		

		if (!params.stop_after_convert) {

			convert_info = run_samestr_convert.out.convert_info
				.join(run_samestr_convert.out.convert_sentinel, by: 0)
					.map { sample, data, sentinel -> return data }
					.collect()

			sstr_convert_buffer("convert", convert_info, params.convert_batch_size)

			grouped_npy_ch = sstr_convert_buffer.out.batches
				.splitCsv(header: ['batch_id', 'batch_size', 'file_path'], sep: '\t' )
				.map { item -> [item.batch_id, item.batch_size, item.file_path] }
				.groupTuple(by: 0)
			samestr_post_convert(grouped_npy_ch, tax_profiles)
		}
}
