#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SEGMENT_TISSUES {
    cpus 1

    input:
        tuple val(sid), path(anat)
    output:
        tuple val(sid), path("${sid}__map_wm.nii.gz"), path("${sid}__map_gm.nii.gz"),
        path("${sid}__map_csf.nii.gz"), emit: maps
        tuple val(sid), path("${sid}__mask_wm.nii.gz"), path("${sid}__mask_gm.nii.gz"),
        path("${sid}__mask_csf.nii.gz"), emit: masks

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    fast -t 1 -n $params.number_of_tissues\
        -H 0.1 -I 4 -l 20.0 -g -o anat.nii.gz $anat
    scil_image_math.py convert anat_seg_2.nii.gz ${sid}__mask_wm.nii.gz --data_type uint8
    scil_image_math.py convert anat_seg_1.nii.gz ${sid}__mask_gm.nii.gz --data_type uint8
    scil_image_math.py convert anat_seg_0.nii.gz ${sid}__mask_csf.nii.gz --data_type uint8
    mv anat_pve_2.nii.gz ${sid}__map_wm.nii.gz
    mv anat_pve_1.nii.gz ${sid}__map_gm.nii.gz
    mv anat_pve_0.nii.gz ${sid}__map_csf.nii.gz
    """    
}

process ATROPOS_SEG {
    cpus 1

    input:
        tuple val(sid), path(anat), path(mask)
    output:
        tuple val(sid), path("${sid}__map_wm.nii.gz"), emit: wm_map
        tuple val(sid), path("${sid}__map_gm.nii.gz"), emit: gm_map
        tuple val(sid), path("${sid}__map_csf.nii.gz"), emit: csf_map
    
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    antsAtroposN4.sh -d 3 -a $anat -x $mask -c 3 \
        -o ${sid}__ -m $params.atropos_m -n $params.atropos_n \
        -b $params.atropos_formulation
    mv ${sid}__SegmentationPosteriors3.nii.gz ${sid}__map_wm.nii.gz
    mv ${sid}__SegmentationPosteriors2.nii.gz ${sid}__map_gm.nii.gz
    mv ${sid}__SegmentationPosteriors1.nii.gz ${sid}__map_csf.nii.gz
    """
}

process GENERATE_MASKS {
    cpus 1

    input:
        tuple val(sid), path(wm_mask), path(fa)
    output:
        tuple val(sid), path("${sid}__seeding_mask.nii.gz"),
        path("${sid}__tracking_mask.nii.gz"), emit: masks
        tuple val(sid), path("${sid}__fa_mask.nii.gz")
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    bet2 $fa fa_bet -m -f 0.16
    scil_image_math.py erosion fa_bet_mask.nii.gz $params.erosion fa_bet_mask.nii.gz -f
    mrcalc fa_bet.nii.gz fa_bet_mask.nii.gz -mult fa_eroded.nii.gz
    mrthreshold fa_eroded.nii.gz ${sid}__fa_mask.nii.gz -abs $params.local_fa_seeding_mask_thr -nthreads 1 -force
    scil_image_math.py union ${sid}__fa_mask.nii.gz $wm_mask\
        ${sid}__seeding_mask.nii.gz --data_type uint8 -f
    cp ${sid}__seeding_mask.nii.gz ${sid}__tracking_mask.nii.gz
    """
}

process LOCAL_TRACKING_MASK {
    cpus 1

    input:
        tuple val(sid), path(wm), path(fa)
    output:
        tuple val(sid), path("${sid}__local_tracking_mask.nii.gz"), emit: tracking_mask
    when:
        params.run_local_tracking
    
    script:
    if (params.local_tracking_mask_type == "wm")
        """
        mv $wm ${sid}__local_tracking_mask.nii.gz
        """
    else if (params.local_tracking_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.local_fa_tracking_mask_thr -ge ${sid}__local_tracking_mask.nii.gz\
            -datatype uint8
        """
}

process LOCAL_SEEDING_MASK {
    cpus 1

    input:
        tuple val(sid), path(wm), path(fa)
    output:
        tuple val(sid), path("${sid}__local_seeding_mask.nii.gz"), emit: seeding_mask
    when:
        params.run_local_tracking
    
    script:
    if (params.local_seeding_mask_type == "wm")
        """
        mv $wm ${sid}__local_seeding_mask.nii.gz
        """
    else if (params.local_seeding_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.local_fa_seeding_mask_thr -ge ${sid}__local_seeding_mask.nii.gz\
            -datatype uint8
        """
}

process LOCAL_TRACKING {
    cpus 2

    input:
        tuple val(sid), path(fodf), path(seeding_mask), path(tracking_mask)
    output:
        tuple val(sid), path("${sid}__local_tracking.trk"), emit: tractogram
    when:
        params.run_local_tracking

    script:
    compress = params.local_compress_streamlines ? '--compress ' + params.local_compress_value : ''
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_local_tracking.py $fodf $seeding_mask $tracking_mask\
        tmp.trk --algo $params.local_algo --$params.local_seeding $params.local_nbr_seeds\
        --seed $params.local_tracking_seed --step $params.local_step_size --theta $params.local_theta\
        --sfthres $params.local_sfthres --min_length $params.local_min_len\
        --max_length $params.local_max_len $compress --sh_basis $params.basis
    scil_remove_invalid_streamlines.py tmp.trk\
        ${sid}__local_tracking.trk --remove_single_point
    """
}

process PFT_SEEDING_MASK {
    cpus 1

    input:
        tuple val(sid), path(wm), path(fa), path(interface_mask)
    output:
        tuple val(sid), path("${sid}__pft_seeding_mask.nii.gz"), emit: seeding_mask
    when:
        params.run_pft_tracking
    
    script:
    if (params.pft_seeding_mask_type == "wm")
        """    
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_image_math.py union $wm $interface_mask ${sid}__pft_seeding_mask.nii.gz\
            --data_type uint8
        """
    else if (params.pft_seeding_mask_type == "interface")
        """
        mv $interface_mask ${sid}__pft_seeding_mask.nii.gz
        """
    else if (params.pft_seeding_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.pft_fa_seeding_mask_thr -ge ${sid}__pft_seeding_mask.nii.gz\
            -datatype uint8
        """
}

process PFT_TRACKING_MASK {
    cpus 1

    input:
        tuple val(sid), path(wm), path(gm), path(csf)
    output:
        tuple val(sid), path("${sid}__map_include.nii.gz"), path("${sid}__map_exclude.nii.gz"), emit: tracking_maps
        tuple val(sid), path("${sid}__interface.nii.gz"), emit: interface_map
    when:
        params.run_pft_tracking

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_maps_for_particle_filter_tracking.py $wm $gm $csf\
        --include ${sid}__map_include.nii.gz\
        --exclude ${sid}__map_exclude.nii.gz\
        --interface ${sid}__interface.nii.gz -f
    """
}

process PFT_TRACKING {
    cpus 2

    input:
        tuple val(sid), path(fodf), path(include), path(exclude), path(seed)

    output:
        tuple val(sid), path("${sid}__pft_tracking.trk"), emit: tractogram
    when:
        params.run_pft_tracking

    script:
    compress = params.pft_compress_streamlines ? '--compress ' + params.pft_compress_value : ''
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_pft.py $fodf $seed $include $exclude\
        tmp.trk\
        --algo $params.pft_algo --$params.pft_seeding $params.pft_nbr_seeds\
        --seed $params.pft_random_seed --step $params.pft_step_size --theta $params.pft_theta\
        --sfthres $params.pft_sfthres --sfthres_init $params.pft_sfthres_init\
        --min_length $params.pft_min_len --max_length $params.pft_max_len\
        --particles $params.pft_particles --back $params.pft_back\
        --forward $params.pft_front $compress --sh_basis $params.basis
    scil_remove_invalid_streamlines.py tmp.trk\
        ${sid}__pft_tracking.trk --remove_single_point
    """
}