include { run_samestr_convert } from "../../../modules/profilers/samestr"
include { sstr_tarball as sstr_convert_tarball } from "../../../modules/profilers/samestr"


workflow samestr_convert {
	take:
		profiler_data  // alignments.join(tax_profiles),

	main:
		run_samestr_convert(
			profiler_data,
			params.samestr_marker_db,
			params.samestr_sqlite
		)

		if (!params.skip_convert_tarball) {
			sstr_convert_tarball("sstr_convert", run_samestr_convert.out.sstr_npy.map { sample, data -> [data].flatten() }.collect())
		}

		// not needed in batched flow
		converted_ch = run_samestr_convert.out.sstr_npy
			.join(run_samestr_convert.out.convert_sentinel, by: 0)
			.map { sample, data, sentinel -> return data }
			.flatten()

		convert_info = run_samestr_convert.out.convert_info
				.join(run_samestr_convert.out.convert_sentinel, by: 0)
				.map { sample, data, sentinel -> return data }
				.collect()
	
	emit:
		data = converted_ch
		batch_info = convert_info
}