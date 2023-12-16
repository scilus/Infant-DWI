#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    BET_T2
} from '../../tracking/processes/preprocess.nf'

include {
    REGISTER_POP;
    REGISTER_FA;
    APPLY_TRANSFORM_DWI_BVECS
} from '../processes/registration.nf'

include {
    AVERAGE_VOLUMES_ANAT;
    AVERAGE_DWI
} from '../processes/average.nf'

workflow POPULATION_TEMPLATE {
    take:
        anat_channel
        dwi_channel
        fa_channel
        anat_ref_channel
        fa_ref_channel
    main:

        BET_T2(anat_channel)

        reg_channel = BET_T2.out.bet_t2
                        .combine(anat_ref_channel)
        
        REGISTER_POP(reg_channel)

        all_anats = REGISTER_POP.out.warped_anat
                        .map{ [it[1]] }
                        .collect()
        
        AVERAGE_VOLUMES_ANAT(all_anats)

        reg_fa_channel = fa_channel
                            .combine(fa_ref_channel)
        
        REGISTER_FA(reg_fa_channel)

        apply_transfo_channel = REGISTER_FA.out.fa_warped
                                .combine(dwi_channel, by: 0)
                                .combine(REGISTER_FA.out.transfos, by: 0)

        APPLY_TRANSFORM_DWI_BVECS(apply_transfo_channel)

        all_dwis = APPLY_TRANSFORM_DWI_BVECS.out.dwi_warped
                        .map{ [it[1]] }
                        .collect()
        
        AVERAGE_DWI(all_dwis)
        
}