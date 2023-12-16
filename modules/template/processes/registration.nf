#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process REGISTER_POP {
    label "REGISTRATION_POP"
    cpus params.processes_registration

    input:
        tuple val(sid), path(anat), path(ref)
    output:
        tuple val(sid), path("${sid}__t2w_warped.nii.gz"), emit: warped_anat
    script:
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
        --initial-moving-transform [ $ref,$anat,1 ] \
        --transform Rigid[ 0.1 ] \
        --metric MI[ $ref,$anat,1,32,Regular,0.25 ] \
        --convergence [ 1000x500x250x100,1e-6,10 ] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
        --transform Affine[ 0.1 ] --metric MI[ $ref,$anat,1,32,Regular,0.25 ] \
        --convergence [ 1000x500x250x100,1e-6,10 ] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
        --transform SyN[ 0.1,3,0 ] \
        --metric CC[ $ref,$anat,1,4 ] \
        --convergence [ 200x150x200x200,1e-6,10 ] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
    mv outputWarped.nii.gz ${sid}__t2w_warped.nii.gz
    """
}

process REGISTER_FA {
    label "REGISTRATION_FA"
    cpus params.processes_registration

    input:
        tuple val(sid), path(moving), path(ref)
    output:
        tuple val(sid), path("${sid}__fa_warped.nii.gz"), emit: fa_warped
        tuple val(sid),
        path("output0GenericAffine.mat"),
        path("output1Warp.nii.gz"), emit: transfos    
    
    script:
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
        --initial-moving-transform [ $ref,$moving,1 ] \
        --transform Rigid[ 0.1 ] \
        --metric MI[ $ref,$moving,1,32,Regular,0.25 ] \
        --convergence [ 1000x500x250x100,1e-6,10 ] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
        --transform Affine[ 0.1 ] --metric MI[ $ref,$moving,1,32,Regular,0.25 ] \
        --convergence [ 1000x500x250x100,1e-6,10 ] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
        --transform SyN[ 0.1,3,0 ] \
        --metric CC[ $ref,$moving,1,4 ] \
        --convergence [ 200x150x200x200,1e-6,10 ] \
        --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
    mv outputWarped.nii.gz ${sid}__fa_warped.nii.gz
    """
}

process APPLY_TRANSFORM_DWI_BVECS {
    label "APPLY_TRANSFORM"
    cpus 1

    input:
        tuple val(sid), path(warped_fa), path(dwi), path(bvec), path(mat), path(warp)
    output:
        tuple val(sid), path("${sid}__warped_dwi.nii.gz"), emit: dwi_warped
        tuple val(sid), path("${sid}__warped_dwi.bvec"), emit: bvec_warped
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export ANTS_RANDOM_SEED=1234
    antsApplyTransforms -d 3 -e 3 \
        -i $dwi -r $warped_fa \
        -n Linear \
        -t $warp $mat \
        -o ${sid}__warped_dwi.nii.gz
    scil_image_math.py convert ${sid}__warped_dwi.nii.gz ${sid}__warped_dwi.nii.gz --data_type float32 -f
    scil_apply_transform_to_bvecs.py $bvec $mat ${sid}__warped_dwi.bvec
    """
}