#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process REGISTER_T2 {
    cpus params.processes_registration

    input:
        tuple val(sid), path(md), path(t2w), path(wm_mask)
    output:
        tuple val(sid), 
        path("${sid}__t2w_warped.nii.gz"),
        path("${sid}__wm_mask_warped.nii.gz"), emit: warped_anat
        tuple val(sid),
        path("output0GenericAffine.mat"),
        path("output1Warp.nii.gz"), emit: transfos
        tuple val(sid),
        path("outputInverseWarped.nii.gz"), emit: inverse_transfo
    when:
        params.infant_config

    script:
    // ** Registration from t2w to diff space in infant returns better result when using the MD map due ** //
    // ** to similar intensities (white = CSF in both volumes). See: Uus A, Pietsch M, Grigorescu I, ** //
    // ** Christiaens D, Tournier JD, Grande LC, Hutter J, Edwards D, Hajnal J, Deprez M. Multi-channel ** //
    // ** Registration for Diffusion MRI: Longitudinal Analysis for the Neonatal Brain. Biomedical Image ** //
    // ** Registration. 2020 May 13;12120:111–21. doi: 10.1007/978-3-030-50120-4_11. PMCID: PMC7279935. ** //
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export ANTS_RANDOM_SEED=1234
    antsRegistration --dimensionality 3 --float 0 \
        --collapse-output-transforms 1 \
        --output [ output,outputWarped.nii.gz,outputInverseWarped.nii.gz ] \
        --interpolation Linear --use-histogram-matching 0 \
        --winsorize-image-intensities [ 0.005,0.995 ] \
        --initial-moving-transform [ $md,$t2w,1 ] \
        --transform Rigid[ 0.1 ] \
        --metric MI[ $md,$t2w,1,32,Regular,0.25 ] \
        --convergence [ 1000x500x250x100,1e-6,10 ] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
        --transform Affine[ 0.1 ] --metric MI[ $md,$t2w,1,32,Regular,0.25 ] \
        --convergence [ 1000x500x250x100,1e-6,10 ] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
        --transform SyN[ 0.1,3,0 ] \
        --metric CC[ $md,$t2w,1,4 ] \
        --convergence [ 200x150x200x200,1e-6,10 ] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
    mv outputWarped.nii.gz ${sid}__t2w_warped.nii.gz
    antsApplyTransforms -d 3 -i $wm_mask -r ${sid}__t2w_warped.nii.gz \
        -o ${sid}__wm_mask_warped.nii.gz -n NearestNeighbor \
        -t output1Warp.nii.gz output0GenericAffine.mat -u int
    """
}

process REGISTER_T1 {
    cpus params.processes_registration

    input:
        tuple val(sid), path(fa), path(t1), path(t1_mask), path(b0)
    output:
        tuple val(sid), path("${sid}__t1_warped.nii.gz"), emit: t1_warped
        tuple val(sid), path("${sid}__output0GenericAffine.mat"),
        path("${sid}__output1Warp.nii.gz"), emit: transfos
        tuple val(sid), path("${sid}__output1InverseWarp.nii.gz"), emit: inverse_transfo
        tuple val(sid), path("${sid}__t1_mask_warped.nii.gz"), emit: t1_mask_warped
    when:
        !params.infant_config
    
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export ANTS_RANDOM_SEED=1234
    antsRegistration --dimensionality 3 --float 0\
        --output [output,outputWarped.nii.gz,outputInverseWarped.nii.gz]\
        --interpolation Linear --use-histogram-matching 0\
        --winsorize-image-intensities [0.005,0.995]\
        --initial-moving-transform [$b0,$t1,1]\
        --transform Rigid['0.2']\
        --metric MI[$b0,$t1,1,32,Regular,0.25]\
        --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
        --smoothing-sigmas 3x2x1x0\
        --transform Affine['0.2']\
        --metric MI[$b0,$t1,1,32,Regular,0.25]\
        --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
        --smoothing-sigmas 3x2x1x0\
        --transform SyN[0.1,3,0]\
        --metric MI[$b0,$t1,1,32]\
        --metric CC[$fa,$t1,1,4]\
        --convergence [50x25x10,1e-6,10] --shrink-factors 4x2x1\
        --smoothing-sigmas 3x2x1
    mv outputWarped.nii.gz ${sid}__t1_warped.nii.gz
    mv output0GenericAffine.mat ${sid}__output0GenericAffine.mat
    mv output1InverseWarp.nii.gz ${sid}__output1InverseWarp.nii.gz
    mv output1Warp.nii.gz ${sid}__output1Warp.nii.gz
    antsApplyTransforms -d 3 -i $t1_mask -r ${sid}__t1_warped.nii.gz\
        -o ${sid}__t1_mask_warped.nii.gz -n NearestNeighbor\
        -t ${sid}__output1Warp.nii.gz ${sid}__output0GenericAffine.mat
    scil_image_math.py convert ${sid}__t1_mask_warped.nii.gz ${sid}__t1_mask_warped.nii.gz\
        --data_type uint8 -f
    """
}