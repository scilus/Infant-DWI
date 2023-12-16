#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRANSFORM_LABELS;
          TRANSFORM_T1 } from "../processes/transform.nf"
include { DECOMPOSE_CONNECTIVITY as INITIAL_DECOMPOSE } from "../processes/decompose.nf"
include { DECOMPOSE_CONNECTIVITY as FINAL_DECOMPOSE } from "../processes/decompose.nf"
include { COMMIT;
          COMMIT_ON_TRK;
          COMPUTE_PRIORS } from "../processes/commit.nf"
include { COMPUTE_AFD_FIXEL;
          COMPUTE_CONNECTIVITY } from "../processes/compute_metrics.nf"
include { VISUALIZE_CONNECTIVITY } from "../processes/viz.nf"

workflow CONNECTOMICS {
    take:
        tracking_channel
        labels_channel
        dwi_peaks_channel
        fodf_channel
        metrics_channel
        anat_channel
        transfos_channel
        fa_md_ad_rd_channel
        
    main:
        // ** Computing priors for COMMIT ** //
        COMPUTE_PRIORS(fa_md_ad_rd_channel)

        // ** If -profile freesurfer, transform t1 to diff space. ** //
        if ( params.run_freesurfer && !params.run_tracking ) {
            t1_for_transfo = anat_channel
                                .combine(dwi_peaks_channel.map{ [it[0], it[1], it[2], it[3]] }, by: 0)
                                .combine(transfos_channel, by: 0)   
            TRANSFORM_T1(t1_for_transfo)
            channel_for_transfo = labels_channel
                                    .combine(TRANSFORM_T1.out.t1_warped, by: 0)
                                    .combine(transfos_channel, by: 0)
        } else {
            channel_for_transfo = labels_channel
                .combine(anat_channel, by: 0)
                .combine(transfos_channel, by: 0)
        }

        // ** Transforming labels to diff space ** //
        TRANSFORM_LABELS(channel_for_transfo)

        // ** If -profile infant is used, first part will be run. COMMIT1 is the only supported ** //
        // ** method as of now, since running commit2 requires a decomposition first, which is not an ** //
        // ** easy task on infant data. This will be improved in the future. ** //
        if ( params.infant_config ) {

            // ** COMMIT1 processing on trk ** //
            commit_channel = tracking_channel
                                .combine(dwi_peaks_channel, by: 0)
            COMMIT_ON_TRK(commit_channel)

            // ** Decomposing tractogram ** //
            decompose_channel = COMMIT_ON_TRK.out.trk_commit
                                .combine(TRANSFORM_LABELS.out.labels_warped, by: 0)
            INITIAL_DECOMPOSE(decompose_channel)

            // ** Setting output channel ** //
            afd_fixel_channel = INITIAL_DECOMPOSE.out.decompose
                                    .combine(fodf_channel, by: 0)
        }
        else {
            // ** Decomposing tractogram ** //
            decompose_channel = tracking_channel
                                    .combine(TRANSFORM_LABELS.out.labels_warped, by: 0)
            INITIAL_DECOMPOSE(decompose_channel)

            // ** Running COMMIT1 or COMMIT2 ** //
            if ( params.use_both ) {

                if ( params.compute_priors ) {
                    commit_channel = INITIAL_DECOMPOSE.out.decompose
                                    .combine(dwi_peaks_channel, by: 0)
                                    .combine(COMPUTE_PRIORS.out.para_diff, by: 0)
                                    .combine(COMPUTE_PRIORS.out.iso_diff, by: 0)
                                    .combine(COMPUTE_PRIORS.out.perp_diff, by: 0)
                    COMMIT(commit_channel)
                } else {
                    commit_channel = INITIAL_DECOMPOSE.out.decompose
                                    .combine(dwi_peaks_channel, by: 0)
                                    .combine(Channel.of([[]]))
                                    .combine(Channel.of([[]]))
                                    .combine(Channel.of([[]]))
                    COMMIT(commit_channel)
                }

            decompose_channel = COMMIT.out.trk_commit
                                    .combine(TRANSFORM_LABELS.out.labels_warped, by: 0)
            FINAL_DECOMPOSE(decompose_channel)

            // ** Setting output channel ** //
            afd_fixel_channel = FINAL_DECOMPOSE.out.decompose
                                    .combine(fodf_channel, by: 0)
            }
            else {
                
                if ( params.compute_priors ) {
                    commit_channel = INITIAL_DECOMPOSE.out.decompose
                                    .combine(dwi_peaks_channel, by: 0)
                                    .combine(COMPUTE_PRIORS.out.para_diff, by: 0)
                                    .combine(COMPUTE_PRIORS.out.iso_diff, by: 0)
                                    .combine(COMPUTE_PRIORS.out.perp_diff, by: 0)
                    COMMIT(commit_channel)
                } else {
                    commit_channel = INITIAL_DECOMPOSE.out.decompose
                                    .combine(dwi_peaks_channel, by: 0)
                                    .combine(Channel.of([[]]))
                                    .combine(Channel.of([[]]))
                                    .combine(Channel.of([[]]))
                    COMMIT(commit_channel)
                }
            // ** Setting output channel ** //
            afd_fixel_channel = COMMIT.out.h5_commit
                                    .combine(fodf_channel, by: 0)
            }
        }

        // ** Computing AFD fixel ** //
        COMPUTE_AFD_FIXEL(afd_fixel_channel)

        // ** Computing Connectivity ** //
        compute_metrics_channel = COMPUTE_AFD_FIXEL.out.decompose_afd
            .combine(TRANSFORM_LABELS.out.labels_warped, by: 0)
            .combine(metrics_channel, by: 0)
        COMPUTE_CONNECTIVITY(compute_metrics_channel)
        
        // ** Visualizing Connectivity ** //
        VISUALIZE_CONNECTIVITY(COMPUTE_CONNECTIVITY.out.metrics)
}