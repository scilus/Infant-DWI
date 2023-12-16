#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    SEGMENT_TISSUES;
    GENERATE_MASKS;
    LOCAL_TRACKING_MASK;
    LOCAL_SEEDING_MASK;
    LOCAL_TRACKING;
    PFT_SEEDING_MASK;
    PFT_TRACKING_MASK;
    PFT_TRACKING
} from '../processes/tracking_processes.nf'

workflow TRACKING {
    take:
        anat_and_mask_channel
        fodf_channel
        fa_channel
    main:
        
        if ( params.infant_config ) {
            // ** Creating masks channel and generating seeding and tracking masks. ** //
            masks_channel = anat_and_mask_channel
                                .map{ [it[0], it[2]] }
                                .combine(fa_channel, by: 0)
            GENERATE_MASKS(masks_channel)

            // ** Performing local tracking. ** //
            tracking_channel = fodf_channel
                                .combine(GENERATE_MASKS.out.masks, by: 0)
            LOCAL_TRACKING(tracking_channel)
        } else {
            // ** Segmenting tissues using fslfast ** //
            anat_channel = anat_and_mask_channel.map{ [it[0], it[1]] }
            SEGMENT_TISSUES(anat_channel)

            // ** If --run_local_tracking is set, this will be run ** //
            // ** Generating seeding and tracking mask ** //
            local_masks_channel = SEGMENT_TISSUES.out.masks
                                    .map{ [it[0], it[1]] }
                                    .combine(fa_channel, by: 0)
            LOCAL_TRACKING_MASK(local_masks_channel)
            LOCAL_SEEDING_MASK(local_masks_channel)

            // ** Performing local tracking ** //
            local_tracking_channel = fodf_channel
                                        .combine(LOCAL_SEEDING_MASK.out.seeding_mask, by: 0)
                                        .combine(LOCAL_TRACKING_MASK.out.tracking_mask, by: 0)
            LOCAL_TRACKING(local_tracking_channel)

            // ** If --run_pft_tracking is set, this will be run ** //
            // ** Creating PFT masks. ** //
            PFT_TRACKING_MASK(SEGMENT_TISSUES.out.maps)
            pft_masks_channel = SEGMENT_TISSUES.out.masks
                                    .map{ [it[0], it[1]] }
                                    .combine(fa_channel, by: 0)
                                    .combine(PFT_TRACKING_MASK.out.interface_map, by: 0)
            PFT_SEEDING_MASK(pft_masks_channel)

            // ** Performing PFT tracking ** //
            pft_tracking_channel = fodf_channel
                                    .combine(PFT_TRACKING_MASK.out.tracking_maps, by: 0)
                                    .combine(PFT_SEEDING_MASK.out.seeding_mask, by: 0)
            PFT_TRACKING(pft_tracking_channel)
        }   

        if ( params.run_local_tracking ) {
            out_channel = LOCAL_TRACKING.out.tractogram
        } else {
            out_channel = PFT_TRACKING.out.tractogram
        }

    emit:
        trk = out_channel
}