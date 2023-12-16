#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process COMPUTE_AFD_FIXEL {
    cpus params.processes_afd_fixel
    memory '2 GB'

    input:
        tuple val(sid), path(h5), path(fodf)
    output:
        tuple val(sid), path("${sid}__decompose_afd_fixel.h5"), emit: decompose_afd
    script:
    """
    scil_compute_mean_fixel_afd_from_hdf5.py $h5 $fodf ${sid}__decompose_afd_fixel.h5 \
        --sh_basis "descoteaux07" --processes $params.processes_afd_fixel
    """
}

process COMPUTE_CONNECTIVITY {
    cpus params.processes_connectivity
    memory '2 GB'

    input:
        tuple val(sid), path(h5), path(labels), path(metrics)
    output:
        tuple val(sid), path("*.npy"), emit: metrics
    script:
    String metrics_list = metrics.join(", ").replace(',', '')
    """
    metrics_args=""
    for metric in $metrics; do
        base_name=\$(basename \${metric})
        metrics_args="\${metrics_args} --metrics \${metric} \$(basename \$base_name .nii.gz).npy"
    done

    scil_compute_connectivity.py $h5 $labels \
        --volume ${sid}__vol.npy --streamline_count ${sid}__sc.npy \
        --length ${sid}__len.npy \$metrics_args --density_weighting \
        --no_self_connection --include_dps ./ \
        --processes $params.processes_connectivity
    """
}