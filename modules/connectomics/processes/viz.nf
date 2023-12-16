#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process VISUALIZE_CONNECTIVITY {
    cpus 1
    memory "2 GB"

    input:
        tuple val(sid), path(npy)
    output:
        tuple val(sid), path("*.png")
    script:
    String npy_list = npy.join(", ").replace(',', '')
    """
    for matrix in $npy_list; do
        scil_visualize_connectivity.py \$matrix \${matrix/.npy/_matrix.png} \
            --name_axis --display_legend --histogram \${matrix/.npy/_hist.png} --nb_bins 50 \
            --exclude_zeros --axis_text_size 5 5
    done
    """
}