process ss_load_convert_tarball {
	input:
		val(procname)
		path(tarball)

	output:
		path("${procname}/*.npz"), emit: sstr_npy
        path("loaded.samestr_${procname}_clades.txt"), emit: batch_info

	script:

	"""
	mkdir -p ${procname}

	tar xf ${tarball} -C ${procname}

	find \$PWD/${procname} -type f -name '*.npz' >  loaded.samestr_${procname}_clades.txt
	"""
}