#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    FREESURFER;
    RECON_SURF
} from '../processes/freesurfer.nf'
include {
    FS_BN_GL_SF;
    BN_CHILD;
    LOBES;
    LAUSANNE
} from '../processes/atlases.nf'

workflow FREESURFERFLOW {
    take:
        anat

    main:

        if ( params.recon_all ) {
        // ** Lauching FreeSurfer Recon-all ** //
        FREESURFER(anat)
        folder_channel = FREESURFER.out.folders
        t1 = FREESURFER.out.final_t1
        } else if ( params.recon_surf ) {
        // ** Launching FastSurfer ** //
        RECON_SURF(anat)
        folder_channel = RECON_SURF.out.folders
        t1 = RECON_SURF.out.final_t1
        }

        // ** Computing FS_BN_GL_SF atlases ** //
        FS_BN_GL_SF(folder_channel)

        // ** Computing BN_CHILD Atlas ** //
        BN_CHILD(folder_channel)

        // ** Computing lobes atlases ** //
        LOBES(folder_channel)

        // ** Computing lausanne atlas ** //
        scales = Channel.from(1,2,3,4,5)
        LAUSANNE(folder_channel,
                scales)

        // ** Reorganizing Lausanne multiscale atlas channel ** //
        lausanne1 = LAUSANNE.out.lausanne_1.map{ [it[0]] }
                        .merge(LAUSANNE.out.lausanne_1.map{ [it[1]] }.flatMap())
        lausanne2 = LAUSANNE.out.lausanne_2.map{ [it[0]] }
                        .merge(LAUSANNE.out.lausanne_2.map{ [it[1]] }.flatMap())
        lausanne3 = LAUSANNE.out.lausanne_3.map{ [it[0]] }
                        .merge(LAUSANNE.out.lausanne_3.map{ [it[1]] }.flatMap())
        lausanne4 = LAUSANNE.out.lausanne_4.map{ [it[0]] }
                        .merge(LAUSANNE.out.lausanne_4.map{ [it[1]] }.flatMap())
        lausanne5 = LAUSANNE.out.lausanne_5.map{ [it[0]] }
                        .merge(LAUSANNE.out.lausanne_5.map{ [it[1]] }.flatMap())

        // ** Work out a way for the user to select which atlas to use. ** //
        // ** Could be cleaner than a bunch of if statements in the future. ** //
        if ( params.use_freesurfer_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = FS_BN_GL_SF.out.freesurfer_dilated
            } else {
                labels = FS_BN_GL_SF.out.freesurfer
            }
        } else if ( params.use_brainnetome_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = FS_BN_GL_SF.out.brainnetome_dilated
            } else {
                labels = FS_BN_GL_SF.out.brainnetome
            }
        } else if ( params.use_brainnetome_child_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = BN_CHILD.out.brainnetome_child_dilated
            } else {
                labels = BN_CHILD.out.brainnetome_child
            }
        } else if ( params.use_glasser_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = FS_BN_GL_SF.out.glasser_dilated
            } else {
                labels = FS_BN_GL_SF.out.glasser
            }
        } else if ( params.use_schaefer_100_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = FS_BN_GL_SF.out.schaefer_100_dilated
            } else {
                labels = FS_BN_GL_SF.out.schaefer_100
            }
        } else if ( params.use_schaefer_200_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = FS_BN_GL_SF.out.schaefer_200_dilated
            } else {
                labels = FS_BN_GL_SF.out.schaefer_200
            }
        } else if ( params.use_schaefer_400_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = FS_BN_GL_SF.out.schaefer_400_dilated    
            } else {
                labels = FS_BN_GL_SF.out.schaefer_400
            }
        } else if ( params.use_lausanne_1_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = lausanne1.map{ [it[0], it[2]] }          
            } else {
                labels = lausanne1.map{ [it[0], it[1]] }
            }
        } else if ( params.use_lausanne_2_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = lausanne2.map{ [it[0], it[2]] }
            } else {
                labels = lausanne2.map{ [it[0], it[1]] }
            }
        } else if ( params.use_lausanne_3_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = lausanne3.map{ [it[0], it[2]] }
            } else {
                labels = lausanne3.map{ [it[0], it[1]] }
            }
        } else if ( params.use_lausanne_4_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = lausanne4.map{ [it[0], it[2]] }                
            } else {
                labels = lausanne4.map{ [it[0], it[1]] }
            }
        } else if ( params.use_lausanne_5_atlas ) {
            if ( params.use_dilated_labels ) {
                labels = lausanne5.map{ [it[0], it[2]] }
            } else {
                labels = lausanne5.map{ [it[0], it[1]] }
            }
        }

    emit:
        labels
        t1
}