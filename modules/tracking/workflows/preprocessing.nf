#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    BET_DWI;
    BET_T2;
    BET_T1;
    DENOISING;
    DENOISE_T1;
    TOPUP;
    EDDY_TOPUP;
    N4;
    N4_T1;
    CROP_DWI;
    CROP_ANAT;
    RESAMPLE_ANAT;
    RESAMPLE_T1;
    NORMALIZE;
    RESAMPLE_DWI;
    EXTRACT_B0;
    DWI_MASK
} from '../processes/preprocess.nf'

workflow DWI {
    take:
        dwi_channel
        rev_channel

    main:

        // ** Bet ** //
        BET_DWI(dwi_channel)

        // ** Denoising ** //
        DENOISING(BET_DWI.out)
        
        // ** Topup ** //
        topup_channel = dwi_channel
            .map{[it[0], it[2], it[3]]}
            .combine(DENOISING.out, by: 0)
            .combine(rev_channel, by: 0)
            .map{ sid, bvals, bvecs, dwi, rev -> tuple(sid, dwi, bvals, bvecs, rev)}
        TOPUP(topup_channel)

        // ** Eddy ** //
        eddy_channel = dwi_channel
            .map{[it[0], it[2], it[3]]}
            .combine(DENOISING.out, by: 0)
            .combine(TOPUP.out.topup_result, by: 0)
            .map{ sid, bvals, bvecs, dwi, corrected_b0s, field, movpar -> tuple(sid, dwi, bvals, bvecs, corrected_b0s, 
                                                                                field, movpar)}
        EDDY_TOPUP(eddy_channel)
        
        // ** N4 ** //
        n4_channel = EDDY_TOPUP.out.dwi_bval_bvec
            .combine(EDDY_TOPUP.out.b0_mask, by: 0)
        N4(n4_channel)

        // ** Crop ** //
        if ( params.skip_dwi_preprocessing ) {
            DWI_MASK(dwi_channel)
            crop_channel = dwi_channel.map{ [it[0], it[1]] }
                                .combine(DWI_MASK.out.dwi_mask, by: 0)
            CROP_DWI(crop_channel)
        } else {
        dwi_crop_channel = N4.out
            .combine(EDDY_TOPUP.out.b0_mask, by: 0)
        CROP_DWI(dwi_crop_channel)
        }

        // ** Normalization ** //
        normalize_channel = CROP_DWI.out.dwi
            .combine(EDDY_TOPUP.out.dwi_bval_bvec.map{[it[0], it[2], it[3]]}, by: 0)
            .combine(CROP_DWI.out.mask, by: 0)
        NORMALIZE(normalize_channel)

        // ** Resampling ** //
        if ( params.skip_dwi_preprocessing ) {
            resample_channel = CROP_DWI.out.dwi
                                    .combine(CROP_DWI.out.mask, by: 0)
            RESAMPLE_DWI(resample_channel)
        } else {
        resample_dwi_channel = NORMALIZE.out.dwi_normalized
            .combine(CROP_DWI.out.mask, by: 0)
        RESAMPLE_DWI(resample_dwi_channel)
        }

        // ** Extracting b0 ** //
        if ( params.skip_dwi_preprocessing ) {
            extract_b0_channel = RESAMPLE_DWI.out.dwi_resampled
                            .combine(dwi_channel.map{ [it[0], it[2], it[3]] }, by: 0)
            EXTRACT_B0(extract_b0_channel)
        } else {
        extract_b0_channel = EDDY_TOPUP.out.dwi_bval_bvec
            .map{[it[0], it[2], it[3]]}
            .combine(RESAMPLE_DWI.out.dwi_resampled, by: 0)
            .map{ sid, bval, bvec, dwi -> tuple(sid, dwi, bval, bvec)}
        EXTRACT_B0(extract_b0_channel)
        }

    emit:
        dwi_bval_bvec = extract_b0_channel
        b0_and_mask = EXTRACT_B0.out.b0_and_mask
}

workflow ANAT {
    take:
        anat_channel

    main: 
        if ( ! params.infant_config ) {
            // ** Denoising ** //
            DENOISE_T1(anat_channel)

            // ** N4 ** //
            N4_T1(DENOISE_T1.out.t1_denoised)  
        }

        // ** Resampling ** //
        if ( params.infant_config ) {
            // ** Resample if -profile infant ** //
            RESAMPLE_ANAT(anat_channel)
        } else {
            RESAMPLE_T1(N4_T1.out.t1_n4)
        }

        // ** Bet ** //
        if ( params.infant_config ) {
            // ** Bet if -profile infant ** //
            BET_T2(RESAMPLE_ANAT.out.t2w_and_mask.map{ [it[0], it[1]] })
        } else {
            BET_T1(RESAMPLE_T1.out.t1_resampled)
        }

        // ** Crop ** //
        if ( params.infant_config ) {
            crop_channel = BET_T2.out.t2_bet
                                .combine(RESAMPLE_ANAT.out.t2w_and_mask.map{ [it[0], it[2]] }, by: 0)
            CROP_ANAT(crop_channel)
        } 
        else {
            CROP_ANAT(BET_T1.out.t1_and_mask_bet)
        }

    emit:
        anat_and_mask = CROP_ANAT.out.cropped_anat_and_mask
}
