#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TRANSFORM_LABELS {
    cpus 1
    memory '2 GB'

    input:
        tuple val(sid), path(labels), path(t2), path(mat), path(syn)
    output:
        tuple val(sid), path("${sid}__labels_warped.nii.gz"), emit: labels_warped
    script:
    """
    antsApplyTransforms -d 3 -i $labels -r $t2 -o ${sid}__labels_warped.nii.gz \
        -t $syn $mat -n NearestNeighbor
    scil_image_math.py convert ${sid}__labels_warped.nii.gz ${sid}__labels_warped.nii.gz \
        --data_type int16 -f
    """
}

process TRANSFORM_T1 {
    cpus 1
    memory '2 GB'

    input:
        tuple val(sid), path(t1), path(dwi), path(bval), path(bvec), path(mat), path(syn)
    output:
        tuple val(sid), path("${sid}__t1_warped.nii.gz"), emit: t1_warped
    script:
    """
    scil_extract_b0.py $dwi $bval $bvec b0.nii.gz --mean\
        --b0_thr $params.b0_thr --force_b0_threshold
    antsApplyTransforms -d 3 -i $t1 -r b0.nii.gz -o ${sid}__t1_warped.nii.gz \
        -t $syn $mat -n Linear
    """
}