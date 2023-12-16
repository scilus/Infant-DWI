#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process DECOMPOSE_CONNECTIVITY {
    cpus 1
    memory { 7.B * trk.size() }

    input:
        tuple val(sid), path(trk), path(labels)
    output:
        tuple val(sid), path("${sid}__decompose.h5"), emit: decompose
    script:
    no_pruning_arg = ""
    if ( params.no_pruning ) {
        no_pruning_arg = "--no_pruning"
    }
    
    no_remove_loops_arg = ""
    if ( params.no_remove_loops ) {
        no_remove_loops_arg = "--no_remove_loops"
    }

    no_remove_outliers_arg = ""
    if ( params.no_pruning ) {
        no_remove_outliers_arg = "--no_pruning"
    }

    no_remove_outliers_arg = ""
    if ( params.no_remove_outliers ) {
        no_remove_outliers_arg = "--no_remove_outliers"
    }
    """
    scil_decompose_connectivity.py $trk $labels ${sid}__decompose.h5 --no_remove_curv_dev \
        $no_pruning_arg $no_remove_loops_arg $no_remove_outliers_arg --min_length $params.min_length --max_length $params.max_length \
        --loop_max_angle $params.loop_max_angle --outlier_threshold $params.outlier_threshold -v
    """
}