include { run_samestr_convert; run_samestr_merge; run_samestr_filter; run_samestr_stats; run_samestr_compare; run_samestr_summarize; collate_samestr_stats } from "../modules/profilers/samestr"

include { failure_guard as convert_failure_guard; failure_guard as filter_failure_guard } from "../modules/profilers/samestr"


workflow samestr_post_merge {
	take:
		ss_merged
		tax_profiles
	main:
		run_samestr_filter(ss_merged, params.samestr_marker_db)

		// for each pre-filtered clade "snp matrix"
		ss_merged
			// join 1 matrix to 1 filter sentinel
			.join(run_samestr_filter.out.sentinel, by: 0, remainder: true)
			.branch {
				// if sentinel is not present (i.e., filter failed) -> FAILURE
				failure: it[3] == null
				// else -> SUCCESS
				success: true
			}
			.set { filter_status_ch }

		// this is just a log operation
		ss_merged
			.join(run_samestr_filter.out.sentinel, by: 0, remainder: true)
			.dump(pretty: true, tag: "filter_status_ch")

		// send filter-failures to filter guard -> if any, tear down pipeline
		filter_failure_guard(filter_status_ch.failure.map { species, data, names, sentinel -> species }, "filter")

		// join filtered matrices to success sentinels
		filtered_ch = run_samestr_filter.out.sstr_npy
			.join(filter_status_ch.success, by: 0)
			// strip off unfiltered data and sentinel
			.map { species, data, names, input_data, input_names, sentinel -> [ species, data, names ] }

		filtered_ch.dump(pretty: true, tag: "filtered_ch")

		// route filtered data to downstream processes
		run_samestr_stats(filtered_ch, params.samestr_marker_db)
		collate_samestr_stats(run_samestr_stats.out.sstr_stats.collect())

		run_samestr_compare(filtered_ch, params.samestr_marker_db)

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

		samestr_post_merge(run_samestr_merge.out.sstr_npy, tax_profiles)
}

workflow samestr_full {

	take:
		alignments
		tax_profiles

	main:
		run_samestr_convert(
			alignments.join(tax_profiles),
			params.samestr_marker_db
		)

		// for each tax profile, check if the convert run was successful
		tax_profiles
			// join 1 tax profile to 1 sentinel
			.join(run_samestr_convert.out.sentinel, by: 0, remainder: true)
			.branch {
				// if sentinel is not present (i.e., convert failed) -> FAILURE
				failure: it[2] == null
				// else -> SUCCESS
				success: true
			}
			.set { convert_status_ch }

		// send all FAILURES to the guard, which will do nothing but report, fail, and tear down the pipeline
		convert_failure_guard(convert_status_ch.failure.map { sample, data, sentinel -> sample.id }, "convert")

		// for each clade-sample collection
		grouped_npy_ch = run_samestr_convert.out.sstr_npy
			// join 1 collection to 1 sentinel
			.join(convert_status_ch.success, by: 0)
			// strip off metadata, tax profiles, and sentinel
			.map { sample, data, tax_profiles, sentinel -> data }
			// flatten collections into stream
			.flatten()
			// add species from filename
			.map { file ->
				def species = file.name.replaceAll(/[.].*/, "")
				return [species, file]
			}
			// group by species
			.groupTuple(sort: true)

		if (!params.stop_after_convert) {
			samestr_post_convert(grouped_npy_ch, tax_profiles)
		}
}
