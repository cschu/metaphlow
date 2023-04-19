process reduce_metaphlan_profiles {

	input:
		path(mp_collated_profiles)
		value(reduce_function)
	
	output:
		path("mp_${reduce_function}_reduced_profiles.txt"), emit: mp_reduced_profiles
	
	script:
		"""
		humann_reduce_table --input ${mp_collated_profiles} --output mp_${reduce_function}_reduced_profiles.txt.tmp --function ${reduce_function} --sort-by level
		"""
		// sed -i "2 s/^/xxx/" mp_${reduce_function}_reduced_profiles.txt.tmp
		// mv mp_${reduce_function}_reduced_profiles.txt.tmp mp_${reduce_function}_reduced_profiles.txt

}


process generate_humann_joint_index {
	
	input:
		path(mp_reduced_profiles)
		path(nuc_db)

	output:
		path("joint_bt2_index/*"), emit: joint_bt2_index

	script:
		"""
		mkdir -p joint_bt2_index/
		echo -e "@dummy\nA\n+\n5" > joint.fastq

		humann --threads ${task.cpus} --input joint.fastq --input-format fastq --output joint_bt2_index/ --bypass-translated-search --taxonomic-profile ${mp_reduced_profiles} --nucleotide-database ${nuc_db}
		"""

}
