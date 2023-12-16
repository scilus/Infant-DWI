#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BET_DWI {
    cpus 2

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
    output:
        tuple val(sid), path("${sid}__dwi_bet.nii.gz"), emit: bet_dwi
    when:
        !params.skip_dwi_preprocessing

    script:
    // ** Using a combination of preliminary bet, powder average computation and then final bet. ** //
    // ** This might not be necessary for good quality data, but returns much more robust results on ** //
    // ** infant data. ** //
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0_mean.nii.gz\
        --b0_thr $params.b0_thr --force_b0_threshold --mean
    bet ${sid}__b0_mean.nii.gz ${sid}__b0_bet -f $params.initial_bet_f -m -R
    scil_image_math.py convert ${sid}__b0_bet_mask.nii.gz ${sid}__b0_bet_mask.nii.gz\
        --data_type uint8 -f
    mrcalc $dwi ${sid}__b0_bet_mask.nii.gz -mult ${sid}__dwi_bet_prelim.nii.gz\
        -quiet -force -nthreads 1
    scil_compute_powder_average.py ${sid}__dwi_bet_prelim.nii.gz $bval\
        ${sid}__powder_avg.nii.gz --b0_thr $params.b0_thr -f
    bet ${sid}__powder_avg.nii.gz ${sid}__powder_avg_bet -m -R -f $params.final_bet_f
    scil_image_math.py convert ${sid}__powder_avg_bet_mask.nii.gz ${sid}__powder_avg_bet_mask.nii.gz\
        --data_type uint8 -f
    mrcalc $dwi ${sid}__powder_avg_bet_mask.nii.gz -mult ${sid}__dwi_bet.nii.gz\
        -quiet -force -nthreads 1
    """
}

process BET_T2 {
    cpus 2

    input:
        tuple val(sid), path(anat)
    output:
        tuple val(sid), path("${sid}__t2w_bet.nii.gz"), emit: t2_bet
    when:
        params.infant_config
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    bet $anat ${sid}__t2w_bet.nii.gz -f $params.bet_anat_f -R
    """
}

process DENOISING {
    cpus params.processes_denoise_dwi

    input:
        tuple val(sid), path(dwi)
    output:
        tuple val(sid), path("${sid}__dwi_denoised.nii.gz"), emit: denoised_dwi
    when:
        !params.skip_dwi_preprocessing
    
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    dwidenoise $dwi ${sid}__dwi_denoised.nii.gz -extent 7 -nthreads 6
    fslmaths ${sid}__dwi_denoised.nii.gz -thr 0 ${sid}__dwi_denoised.nii.gz
    """
}

process TOPUP {
    cpus 4

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(revb0)
    output:
        tuple val(sid), path("${sid}__corrected_b0s.nii.gz"), path("${params.topup_prefix}_fieldcoef.nii.gz"),
        path("${params.topup_prefix}_movpar.txt"), emit: topup_result
    when:
        !params.skip_dwi_preprocessing

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=$task.cpus
    export OPENBLAS_NUM_THREADS=1
    export ANTS_RANDOM_SEED=1234
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0_mean.nii.gz --mean\
        --b0_thr $params.b0_thr --force_b0_threshold
    bet2 $revb0 ${sid}__revb0_bet -m -f $params.initial_bet_f
    antsRegistrationSyN.sh -d 3 -f ${sid}__b0_mean.nii.gz -m ${sid}__revb0_bet.nii.gz\
        -o output
    mv outputWarped.nii.gz ${sid}__rev_b0_warped.nii.gz
    scil_prepare_topup_command.py ${sid}__b0_mean.nii.gz ${sid}__rev_b0_warped.nii.gz\
        --config $params.topup_config\
        --encoding_direction $params.encoding_direction\
        --readout $params.readout --out_prefix $params.topup_prefix\
        --out_script
    bash topup.sh
    mv corrected_b0s.nii.gz ${sid}__corrected_b0s.nii.gz
    """
}

process EDDY_TOPUP {
    cpus params.processes_eddy
    memory { 5.GB * task.attempt }

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(b0s_corrected), path(field), path(movpar)
    output:
        tuple val(sid), path("${sid}__dwi_corrected.nii.gz"), path("${sid}__bval_eddy"), 
        path("${sid}__dwi_eddy_corrected.bvec"), emit: dwi_bval_bvec
        tuple val(sid), path("${sid}__b0_bet_mask.nii.gz"), emit: b0_mask
    when:
        !params.skip_dwi_preprocessing

    script:
    slice_drop_flag=""
    if (params.use_slice_drop_correction)
        slice_drop_flag="--slice_drop_correction"
    """
    export OMP_NUM_THREADS=$task.cpus
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OPENBLAS_NUM_THREADS=1
    mrconvert $b0s_corrected b0_corrected.nii.gz -coord 3 0 -axes 0,1,2 -nthreads 1 -quiet
    bet2 b0_corrected.nii.gz ${sid}__b0_bet -m\
        -f $params.topup_bet_f
    scil_prepare_eddy_command.py $dwi $bval $bvec ${sid}__b0_bet_mask.nii.gz\
        --topup $params.topup_prefix --eddy_cmd $params.eddy_cmd\
        --b0_thr $params.b0_thr\
        --encoding_direction $params.encoding_direction\
        --readout $params.readout --out_script --fix_seed\
        $slice_drop_flag
    echo "--very_verbose" >> eddy.sh
    bash eddy.sh
    fslmaths dwi_eddy_corrected.nii.gz -thr 0 ${sid}__dwi_corrected.nii.gz
    mv dwi_eddy_corrected.eddy_rotated_bvecs ${sid}__dwi_eddy_corrected.bvec
    mv $bval ${sid}__bval_eddy
    """
}

process N4 {
    cpus 1

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(b0_mask)
    output:
        tuple val(sid), path("${sid}__dwi_n4.nii.gz"), emit: dwi_n4
    when:
        !params.skip_dwi_preprocessing

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0.nii.gz --mean\
        --b0_thr $params.b0_thr --force_b0_threshold
    scil_image_math.py convert $b0_mask ${sid}__b0_mask.nii.gz --data_type uint8 -f
    N4BiasFieldCorrection -i ${sid}__b0.nii.gz\
        -o [${sid}__b0_n4.nii.gz, bias_field_b0.nii.gz]\
        -c [500x250x125x75, 1e-6] -v 1
    scil_apply_bias_field_on_dwi.py $dwi bias_field_b0.nii.gz\
        ${sid}__dwi_n4.nii.gz --mask ${sid}__b0_mask.nii.gz -f
    """
}

process CROP_DWI {
    cpus 1

    input:
        tuple val(sid), path(dwi), path(b0_mask)
    output:
        tuple val(sid), path("${sid}__dwi_cropped.nii.gz"), emit: dwi
        tuple val(sid), path("${sid}__b0_mask_cropped.nii.gz"), emit: mask
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    mrcalc $dwi $b0_mask -mult ${sid}__dwi_temp.nii.gz -quiet -force -nthreads 1
    scil_crop_volume.py ${sid}__dwi_temp.nii.gz ${sid}__dwi_cropped.nii.gz\
        --output_bbox dwi_boundingBox.pkl -f
    scil_crop_volume.py $b0_mask ${sid}__b0_mask_cropped.nii.gz\
        --input_bbox dwi_boundingBox.pkl -f
    scil_image_math.py convert ${sid}__b0_mask_cropped.nii.gz ${sid}__b0_mask_cropped.nii.gz\
        --data_type uint8 -f
    """
}

process DENOISE_T1 {
    cpus params.processes_denoise_t1

    input:
        tuple val(sid), path(t1)
    output:
        tuple val(sid), path("${sid}__t1_denoised.nii.gz"), emit: t1_denoised
    when:
        !params.infant_config
    
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_run_nlmeans.py $t1 ${sid}__t1_denoised.nii.gz 1\
        --processes $task.cpus -f
    """
}

process N4_T1 {
    cpus 1

    input:
        tuple val(sid), path(t1)
    output:
        tuple val(sid), path("${sid}__t1_n4.nii.gz"), emit: t1_n4
    when:
        !params.infant_config
    
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    N4BiasFieldCorrection -i $t1\
        -o [${sid}__t1_n4.nii.gz, bias_field_t1.nii.gz]\
        -c [300x150x75x50, 1e-6] -v 1
    """
}

process CROP_ANAT {
    cpus 1

    input:
        tuple val(sid), path(anat), path(mask)
    output:
        tuple val(sid), 
        path("${sid}__anat_cropped.nii.gz"), 
        path("${sid}__mask_cropped.nii.gz"), emit: cropped_anat_and_mask
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_crop_volume.py $anat ${sid}__anat_cropped.nii.gz\
        --output_bbox boundingBox.pkl -f
    scil_crop_volume.py $mask ${sid}__mask_cropped.nii.gz\
        --input_bbox boundingBox.pkl -f
    """
}

process RESAMPLE_T1 {
    cpus 1

    input:
        tuple val(sid), path(t1)
    output:
        tuple val(sid), path("${sid}__t1_resampled.nii.gz"), emit: t1_resampled
    when:
        !params.infant_config

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_resample_volume.py $t1 ${sid}__t1_resampled.nii.gz\
        --voxel_size $params.anat_resolution \
        --interp $params.anat_interpolation
    """
}

process BET_T1 {
    cpus params.processes_bet_t1

    input:
        tuple val(sid), path(t1)
    output:
        tuple val(sid), 
        path("${sid}__t1_bet.nii.gz"), 
        path("${sid}__t1_bet_mask.nii.gz"), emit: t1_and_mask_bet
    when:
        !params.infant_config
    
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export ANTS_RANDOM_SEED=1234
    antsBrainExtraction.sh -d 3 -a $t1 -e $params.template_t1/t1_template.nii.gz\
        -o bet/ -m $params.template_t1/t1_brain_probability_map.nii.gz -u 0
    scil_image_math.py convert bet/BrainExtractionMask.nii.gz ${sid}__t1_bet_mask.nii.gz\
        --data_type uint8
    mrcalc $t1 ${sid}__t1_bet_mask.nii.gz -mult ${sid}__t1_bet.nii.gz -nthreads 1
    """
}

process RESAMPLE_ANAT {
    cpus 1

    input:
        tuple val(sid), path(t2w), path(mask)
    output:
        tuple val(sid), path("${sid}__t2w_resampled.nii.gz"), path("${sid}__mask_resampled.nii.gz"), emit: t2w_and_mask
    when:
        params.infant_config
        
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_resample_volume.py $t2w ${sid}__t2w_resampled.nii.gz\
        --voxel_size $params.anat_resolution --interp $params.anat_interpolation -f
    scil_resample_volume.py $mask ${sid}__mask_resampled.nii.gz\
        --voxel_size $params.anat_resolution --interp $params.mask_interpolation\
        -f
    scil_image_math.py convert ${sid}__mask_resampled.nii.gz ${sid}__mask_resampled.nii.gz\
        --data_type uint8 -f
    """
}

process NORMALIZE {
    cpus 3

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(b0_mask)
    output:
        tuple val(sid), path("${sid}__dwi_normalized.nii.gz"), emit: dwi_normalized
    when:
        !params.skip_dwi_preprocessing

    script:
    if (params.dti_shells)
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_dwi_shell.py $dwi $bval $bvec $params.dti_shells\
        dwi_dti.nii.gz bval_dti bvec_dti -t $params.dwi_shell_tolerance
    scil_compute_dti_metrics.py dwi_dti.nii.gz bval_dti bvec_dti --mask $b0_mask\
        --not_all --fa fa.nii.gz --force_b0_threshold
    mrthreshold fa.nii.gz fa_wm_mask.nii.gz -abs $params.fa_mask_threshold -nthreads 1
    dwinormalise individual $dwi fa_wm_mask.nii.gz ${sid}__dwi_normalized.nii.gz\
        -fslgrad $bvec $bval -nthreads 1
    """
    else
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1

    shells=\$(awk -v max="$params.max_dti_shell_value" '{for (i = 1; i <= NF; i++) {v = int(\$i);if (v <= max) shells[v] = 1;}}END {for (v in shells) print v;}' "$bval" |\
                sort -n | tr '\n' ' ')
    
    scil_extract_dwi_shell.py $dwi $bval $bvec \$shells\
        dwi_dti.nii.gz bval_dti bvec_dti -t $params.dwi_shell_tolerance
    scil_compute_dti_metrics.py dwi_dti.nii.gz bval_dti bvec_dti --mask $b0_mask\
        --not_all --fa fa.nii.gz --force_b0_threshold
    mrthreshold fa.nii.gz fa_wm_mask.nii.gz -abs $params.fa_mask_threshold -nthreads 1
    dwinormalise individual $dwi fa_wm_mask.nii.gz ${sid}__dwi_normalized.nii.gz\
        -fslgrad $bvec $bval -nthreads 1
    """
}

process RESAMPLE_DWI {
    cpus 3

    input:
        tuple val(sid), path(dwi), path(mask)
    output:
        tuple val(sid), path("${sid}__dwi_resampled.nii.gz"), emit: dwi_resampled
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_resample_volume.py $dwi dwi_resample.nii.gz\
        --voxel_size $params.dwi_resolution --interp $params.dwi_interpolation -f
    fslmaths dwi_resample.nii.gz -thr 0 dwi_resample_clipped.nii.gz
    scil_resample_volume.py $mask mask_resample.nii.gz\
        --ref dwi_resample.nii.gz\
        --enforce_dimensions\
        --interp $params.mask_dwi_interpolation -f
    mrcalc dwi_resample_clipped.nii.gz mask_resample.nii.gz\
        -mult ${sid}__dwi_resampled.nii.gz -quiet -nthreads 1
    """
}

process EXTRACT_B0 {
    cpus 3

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
    output:
        tuple val(sid), path("${sid}__b0_resampled.nii.gz"), path("${sid}__b0_mask_resampled.nii.gz"), emit: b0_and_mask
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0_resampled.nii.gz --mean\
        --b0_thr $params.b0_thr --force_b0_threshold
    mrthreshold ${sid}__b0_resampled.nii.gz ${sid}__b0_mask_resampled.nii.gz\
        --abs 0.00001 -nthreads 1
    """
}

process DWI_MASK {
    cpus 1

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
    output:
        tuple val(sid), path("${sid}__b0_mask.nii.gz"), emit: dwi_mask
    
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_b0.py $dwi $bval $bvec b0.nii.gz --mean\
        --b0_thr $params.b0_thr --force_b0_threshold
    mrthreshold b0.nii.gz ${sid}__b0_mask.nii.gz\
        --abs 0.00001 -nthreads 1
    """
}