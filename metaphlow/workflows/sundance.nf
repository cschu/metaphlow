include { calc_position_stats; sfacts_metagenotype } from "../modules/profilers/sundance"


workflow sundance {

	take:
	sstr_filter_ch

	main:

	sstr_filter_ch = sstr_filter_ch
		.map { clade, npz, names -> [ params.study, clade, npz ] }

	calc_position_stats(sstr_filter_ch)

	sfacts_metagenotype(sstr_filter_ch, calc_position_stats.out.position_stats, params.sundance_global_positions)

}