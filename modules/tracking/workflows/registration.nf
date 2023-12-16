#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    REGISTER_T1;
    REGISTER_T2
} from '../processes/registration_processes.nf'

workflow REGISTRATION {
    take:
        fa_md_channel
        anat_and_mask
        b0_channel
    main:
        
        // ** If -profile infant is selected, will do registration from t2w on MD. ** //
        t2_reg_channel = fa_md_channel
                            .map{ [it[0], it[2]] }
                            .combine(anat_and_mask, by: 0)
        REGISTER_T2(t2_reg_channel)

        // ** Classical registration from t1w to b0/FA. ** //
        t1_reg_channel = fa_md_channel
                            .map{ [it[0], it[1]] }
                            .combine(anat_and_mask, by: 0)
                            .combine(b0_channel, by: 0)
        REGISTER_T1(t1_reg_channel)

        // ** Organising channel for output. ** //
        if ( params.infant_config ) {
            warped_anat = REGISTER_T2.out.warped_anat
            transfos = REGISTER_T2.out.transfos
        } else {
            warped_anat = REGISTER_T1.out.t1_warped
            transfos = REGISTER_T1.out.transfos
        }

    emit:
        warped_anat = warped_anat
        transfos = transfos
}