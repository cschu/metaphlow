process ss_load_convert_tarball {
	input:
		val(procname)
		path(tarball)

	output:
		path("${procname}/*.npz"), emit: sstr_npy
        path("loaded.${procname}_clades.txt"), emit: batch_info

	script:

	"""
	mkdir -p ${procname}

	tar xf ${tarball} -C ${procname}

	find \$PWD/${procname} -type f -name '*.npz' > loaded.${procname}_clades.txt
	"""
}

process ss_load_manifest_tarball {
	input:
		val(procname)
		path(tarball)
	
	output:
		path("${procname}/*.npz"), emit: sstr_npy
        path("loaded.${procname}_clades.txt"), emit: batch_info

	script:

	"""
	mkdir -p ${procname} 

	tar xf ${tarball} -C ${procname}

	mv -v ${procname}/${procname}_manifest.txt .

	awk -F '\\t' 'BEGIN {fout = "";} {if (\$1 != fout) {if (fout != "") close(fout); fout = \$1;} print \$2 >> \$1 } END {if (fout != "") close(fout);}' ${procname}_manifest.txt

	find \$(pwd)/${procname}/ -type f -name '*.npz' > loaded.${procname}_clades.txt
    find \$(pwd)/${procname}/ -type f -name '*.names.txt' >> loaded.${procname}_clades.txt
	"""

}