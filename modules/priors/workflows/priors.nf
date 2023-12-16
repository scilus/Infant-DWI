#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DWI_MASK } from "../../tracking/processes/preprocess.nf"
include { EXTRACT_DTI_SHELL;
          DTI_METRICS } from "../../tracking/processes/DTI_processes.nf"
include { COMPUTE_PRIORS;
          AVERAGE_PRIORS } from "../../connectomics/processes/commit.nf"
include { COMPUTE_FRF;
          MEAN_FRF } from "../../tracking/processes/FODF_processes.nf"

workflow PRIORS {
    take:
        dwi_channel
    
    main:

    // ** Generate the mask ** //
    DWI_MASK(dwi_channel)

    // ** Extract the DTI shell ** //
    EXTRACT_DTI_SHELL(dwi_channel)

    // ** Compute the DTI metrics ** //
    dti_channel = EXTRACT_DTI_SHELL.out.dti_files
                    .combine(DWI_MASK.out.dwi_mask, by: 0)
    DTI_METRICS(dti_channel)

    // ** Compute the priors ** //
    priors_channel = DTI_METRICS.out.fa_and_md
                        .combine(DTI_METRICS.out.ad_and_rd, by: 0)
    COMPUTE_PRIORS(priors_channel)

    // ** AVERAGE THE PRIORS ** //
    avg_priors_channel = COMPUTE_PRIORS.out.para_diff
                            .join(COMPUTE_PRIORS.out.iso_diff)
                            .join(COMPUTE_PRIORS.out.perp_diff)
                            .groupTuple()
    
    AVERAGE_PRIORS(avg_priors_channel)

    // ** Compute the FRF ** //
    frf_channel = dwi_channel
                    .combine(DWI_MASK.out.dwi_mask, by: 0)
    COMPUTE_FRF(frf_channel)

    // ** Compute the mean FRF ** //
    all_frf = COMPUTE_FRF.out.frf
                .map{[it[1]]}
                .collect()
    MEAN_FRF(all_frf)

}
