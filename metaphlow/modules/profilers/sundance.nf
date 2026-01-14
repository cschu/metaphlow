process calc_position_stats {
	container "registry.git.embl.org/schudoma/sundance-docker:latest"
	publishDir params.output_dir, mode: "copy"

	input:
	tuple val(study), val(clade), path(ss_filter_npz)

	output:
	tuple val(study), val(clade), path("sundance/${clade}/${study}.nc"), emit: position_stats

	script:
	"""
	mkdir -p sundance/${clade}/
	calculate_position_stats.py -i ${ss_filter_npz} -o sundance/${clade}/${study}.nc
	"""

}


process sfacts_metagenotype {
	container "registry.git.embl.org/schudoma/sundance-docker:latest"
	publishDir params.output_dir, mode: "copy"

	input:
	tuple val(study), val(clade), path(ss_filter_npz)
	tuple val(study), val(clade), path(position_stats)
	path(global_positions)

	output:
	tuple val(study), val(clade), path("metagenotype/${study}/${clade}.nc")

	script:
	"""
	mkdir -p metagenotype/${study}
	
	run_metagenotype.py -i ${ss_filter_npz} -o metagenotype/${study}/${clade}.nc -s ${position_stats} -p ${global_positions}/${clade}.tsv
	"""



}


-i /g/bork6/schudoma/data_processing/samestr/issue_2867/studies/790/output/sstr_filter/t__SGB17248.npz \
-o /g/bork6/salim/projects/strain_sandbox/metagenotype/study_790/t__SGB17248.nc \
-s /g/bork6/salim/projects/strain_sandbox/position_stats/t__SGB17248.tar.xz \
-m g/bork6/salim/projects/strain_sandbox/position_stats/t__SGB17248/study_790.nc \
-p /g/bork6/salim/projects/strain_sandbox/global_position/t__SGB17248.tsv