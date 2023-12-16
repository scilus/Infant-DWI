#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SH_FITTING_SHELL {
    cpus 3

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
    output:
        tuple val(sid), path("${sid}__dwi_sh_fitting.nii.gz"), 
        path("${sid}__bval_sh_fitting"), path("${sid}__bvec_sh_fitting"), emit: dwi_sh
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_dwi_shell.py $dwi $bval $bvec $params.sh_fitting_shell\
        ${sid}__dwi_sh_fitting.nii.gz ${sid}__bval_sh_fitting ${sid}__bvec_sh_fitting\
        -t $params.dwi_shell_tolerance -f
    """
}

process SH_FITTING {
    cpus 1
    
    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
    output:
        tuple val(sid), path("${sid}__dwi_sh.nii.gz")
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_sh_from_signal.py --sh_order $params.sh_fitting_order --sh_basis $params.sh_fitting_basis\
        $dwi $bval $bvec ${sid}__dwi_sh.nii.gz
    """
}