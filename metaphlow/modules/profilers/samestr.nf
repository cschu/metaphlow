process run_samestr_convert {
    container "ghcr.io/danielpodlesny/samestr:v1.2025.102"
    tag "${sample.id}"
    label "large"
    label "samestr"

    
    input:
		tuple val(sample), path(mp_sam), path(mp_profile)
		path(marker_db)
        path(marker_sqlite)

    output:
        tuple val(sample), path("sstr_convert/*/*.npz"), emit: sstr_npy, optional: true
        tuple val(sample), path("samestr_convert_DONE"), emit: convert_sentinel
        tuple val(sample), path("${sample.id}.samestr_convert_clades.txt"), emit: convert_info

    script:
    def keep_intermediates = (params.debug_convert) ? "--keep-tmp-files" : ""
    def sqlite_db = (params.db_dev) ? "--sqlitedb samestr.db.copy" : ""
    def prep_db = (params.db_dev && params.copy_sqlite) ? "cp -v ${marker_sqlite} samestr.db.copy" : "ln -s ${marker_sqlite} samestr.db.copy"
    def delete_db = (params.db_dev) ? "rm -fv samestr.db.copy" : ""

    """
    set -e -o pipefail

    ${prep_db}

    samestr --verbosity DEBUG \
    convert ${keep_intermediates} ${sqlite_db} \
        --input-files ${mp_sam} \
        --min-vcov 1 \
        --min-aln-qual 0 \
        --marker-dir ${marker_db} \
        --output-dir sstr_convert/ \
        --nprocs ${task.cpus} \
        --tax-profiles-extension .txt
    
    ${delete_db}

    find \$(pwd)/sstr_convert/ -type f -name '*.npz' > ${sample.id}.samestr_convert_clades.txt

    touch samestr_convert_DONE
    """
}

process run_samestr_merge {
    publishDir params.output_dir, mode: "copy"
    container "ghcr.io/danielpodlesny/samestr:v1.2025.102"
    // tag "${species}"
    tag "batch_${batch_id}"
    label "large"
    label "samestr"
    
    input:
        tuple val(batch_id), path(sstr_npy)
        // tuple val(species), path(sstr_npy)
	    path(marker_db)

    output:
        // tuple \
        //     val(species), \
        //     path("sstr_merge/${species}.npz"), \
        //     path("sstr_merge/${species}.names.txt"), \
        path("sstr_merge/*.{npz,names.txt}"), emit: sstr_npy

    script:
    """
    samestr --verbosity DEBUG \
    merge \
        --input-files ${sstr_npy} \
        --output-dir sstr_merge/ \
	    --marker-dir ${marker_db} \
        --nprocs ${task.cpus}
    """
        // --clade ${species} \
}

process run_samestr_filter {
    container "ghcr.io/danielpodlesny/samestr:v1.2025.102"
    tag "batch_${batch_id}"
    label "large"
    label "samestr"
    
    input:
        // tuple val(species), path(sstr_npy), path(sstr_names)
        tuple val(batch_id), path(files)
	    path(marker_db)
        path(marker_sqlite)

    output:
            // val(species), \
        tuple \
            val(batch_id),
            path("sstr_filter/*.npz"), \
            path("sstr_filter/*.names.txt"), \
        emit: sstr_npy, optional: true

    script:
    // #    --global-pos-min-n-vcov 10 \
    // #    --sample-pos-min-n-vcov 2 \
    def sqlite_db = (params.db_dev) ? "--sqlitedb samestr.db.copy" : ""
    def prep_db = (params.db_dev && params.copy_sqlite) ? "cp -v ${marker_sqlite} samestr.db.copy" : "ln -s ${marker_sqlite} samestr.db.copy"
    def delete_db = (params.db_dev) ? "rm -fv samestr.db.copy" : ""

    """
    ${prep_db}

    samestr --verbosity DEBUG \
    filter ${sqlite_db} \
        --input-files *.npz \
        --input-names *.names.txt \
        --output-dir sstr_filter/ \
        --marker-dir ${marker_db} \
        --marker-trunc-len 50 \
        --global-pos-min-f-vcov 0.25 \
        --sample-pos-min-sd-vcov 3 \
        --samples-min-n-hcov 5000 \
        --sample-var-min-n-vcov 2 \
        --sample-var-min-f-vcov 0.025 \
        --clade-min-samples 1 \
        --nprocs ${task.cpus}

    ${delete_db}
    """
}

process run_samestr_stats {
    publishDir params.output_dir, mode: "copy"
    container "ghcr.io/danielpodlesny/samestr:v1.2025.102"
    tag "batch_${batch_id}"
    label "large"
    label "samestr"
    
    input:
        // tuple val(species), path(sstr_npy), path(sstr_names)
        tuple val(batch_id), path(ssty_npy), path(sstr_names)
	    path(marker_db)

    output:
        // path "sstr_stats/${species}.aln_stats.txt", emit: sstr_stats
        path "sstr_stats/*.aln_stats.txt", emit: sstr_stats

    script:
    """
    samestr --verbosity DEBUG \
    stats \
    --input-files *.npz \
    --input-names *.names.txt \
    --marker-dir ${marker_db} \
    --nprocs ${task.cpus} \
    --output-dir sstr_stats/
    """
}

process collate_samestr_stats {
    publishDir params.output_dir, mode: "copy"
    tag "ohm..nom..nom!"
    label "large"

    input:
    path(stats_files)

    output:
    path("sstr_aln_stats.collated.tsv.gz"), emit: sstr_stats

    script:
    // find . -maxdepth 1 -mindepth 1 -type f -name '*.aln_stats.txt' | sort | xargs -I {} awk -F '\t' -v OFS='\t' -v clade={} 'NR>1 { print gensub(/\.aln_stats.txt/, "", "g", gensub(/.+\//, "", "g", clade)),$0 }' {} | head
    """
    head -n 1 ${stats_files[0]} | sed "s/^/Clade\\t/" > sstr_aln_stats.collated.tsv
    find . -maxdepth 1 -mindepth 1 -name '*.aln_stats.txt' | sort | xargs -I {} awk -F '\\t' -v OFS='\\t' -v clade={} 'NR>1 { print gensub(/\\.aln_stats.txt/, "", "g", gensub(/.+\\//, "", "g", clade)),\$0; }' {} >> sstr_aln_stats.collated.tsv
    gzip sstr_aln_stats.collated.tsv
    """
    
}

process run_samestr_compare {
    publishDir params.output_dir, mode: "copy"
    container "ghcr.io/danielpodlesny/samestr:v1.2025.102"
    tag "${species}"
    label "large"
    label "samestr"
    
    input:
        tuple val(species), path(sstr_npy), path(sstr_names)
	    path(marker_db)

    output:
        tuple \
            path("sstr_compare/${species}.closest.txt"), \
            path("sstr_compare/${species}.fraction.txt"), \
            path("sstr_compare/${species}.overlap.txt"), \
        emit: sstr_compare

    script:
    """
    samestr --verbosity DEBUG \
    compare \
        --input-files ${sstr_npy} \
        --input-names ${sstr_names} \
        --marker-dir ${marker_db} \
        --output-dir sstr_compare/ \
        --nprocs ${task.cpus}
    """
}

process run_samestr_summarize {
    publishDir params.output_dir, mode: "copy"
    container "ghcr.io/danielpodlesny/samestr:v1.2025.102"
    label "large"
    label "samestr"
    
    input:
        path(sstr_data)
        path(mp_profiles)
	path(marker_db)

    output:
        tuple \
            path("sstr_summarize/taxon_counts.tsv"), \
            path("sstr_summarize/sstr_cooccurrences.tsv"), \
            path("sstr_summarize/sstr_strain_events.tsv"), \
        emit: sstr_summarize

    script:
    """
    mkdir profiles/
    # TODO: this will only work with metaphlan, not motus
    find . -maxdepth 1  \\( -name '*.mp4.txt' -o -name '*.motus.txt' \\) -exec mv -v {} profiles/ \\;

    samestr --verbosity DEBUG \
    summarize \
        --input-dir ./ \
        --marker-dir ${marker_db} \
        --tax-profiles-dir ./profiles/ \
        --tax-profiles-extension .txt \
        --output-dir sstr_summarize/
    """
}


process sstr_tarball {
    label "samestr_tarball"
    tag "${procname}"
    publishDir params.output_dir, mode: "copy"

    input:
    val(procname)
    path(files), name: "input/*"

    output:
    path("*.tar.gz")

    script:
    """
    mv -v input ${procname}
    tar chvzf ${procname}.tar.gz ${procname}
    """
}