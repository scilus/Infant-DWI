#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process EXTRACT_DTI_SHELL {
    cpus 3

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
    output:
        tuple val(sid), path("${sid}__dwi_dti.nii.gz"), path("${sid}__bval_dti"),
        path("${sid}__bvec_dti"), emit: dti_files
    script:
    if (params.dti_shells)
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_extract_dwi_shell.py $dwi \
            $bval $bvec $params.dti_shells ${sid}__dwi_dti.nii.gz \
            ${sid}__bval_dti ${sid}__bvec_dti -t $params.dwi_shell_tolerance -f
    """
    else
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1

        shells=\$(awk -v max="$params.max_dti_shell_value" '{for (i = 1; i <= NF; i++) {v = int(\$i);if (v <= max) shells[v] = 1;}}END {for (v in shells) print v;}' "$bval" |\
                sort -n | tr '\n' ' ')

        scil_extract_dwi_shell.py $dwi \
            $bval $bvec \$shells ${sid}__dwi_dti.nii.gz \
            ${sid}__bval_dti ${sid}__bvec_dti -t $params.dwi_shell_tolerance -f
    """
}

process DTI_METRICS {
    cpus 3

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(b0_mask)
    output:
        tuple val(sid), path("${sid}__fa.nii.gz"), path("${sid}__md.nii.gz"), emit: fa_and_md
        tuple val(sid),
        path("${sid}__ad.nii.gz"),
        path("${sid}__rd.nii.gz"), emit: ad_and_rd
        tuple val(sid),
        path("${sid}__evecs.nii.gz"),
        path("${sid}__evecs_v1.nii.gz"),
        path("${sid}__evecs_v2.nii.gz"),
        path("${sid}__evecs_v3.nii.gz"),
        path("${sid}__evals.nii.gz"),
        path("${sid}__evals_e1.nii.gz"),
        path("${sid}__evals_e2.nii.gz"),
        path("${sid}__evals_e3.nii.gz"),
        path("${sid}__ga.nii.gz"),
        path("${sid}__rgb.nii.gz"),
        path("${sid}__mode.nii.gz"),
        path("${sid}__norm.nii.gz"),
        path("${sid}__tensor.nii.gz"),
        path("${sid}__nonphysical.nii.gz"),
        path("${sid}__pulsation_std_dwi.nii.gz"),
        path("${sid}__residual.nii.gz"), 
        path("${sid}__residual_iqr_residuals.npy"),
        path("${sid}__residual_mean_residuals.npy"),
        path("${sid}__residual_q1_residuals.npy"),
        path("${sid}__residual_q3_residuals.npy"),
        path("${sid}__residual_residuals_stats.png"),
        path("${sid}__residual_std_residuals.npy"), emit: dti_metrics
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_dti_metrics.py $dwi $bval $bvec --mask $b0_mask\
        --ad ${sid}__ad.nii.gz --evecs ${sid}__evecs.nii.gz\
        --evals ${sid}__evals.nii.gz --fa ${sid}__fa.nii.gz\
        --ga ${sid}__ga.nii.gz --rgb ${sid}__rgb.nii.gz\
        --md ${sid}__md.nii.gz --mode ${sid}__mode.nii.gz\
        --norm ${sid}__norm.nii.gz --rd ${sid}__rd.nii.gz\
        --tensor ${sid}__tensor.nii.gz\
        --non-physical ${sid}__nonphysical.nii.gz\
        --pulsation ${sid}__pulsation.nii.gz\
        --residual ${sid}__residual.nii.gz\
        -f --force_b0_threshold
    """
}