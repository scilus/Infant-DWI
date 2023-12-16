#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = false
params.references = false

def fetch_id ( dir, dir_base ) {
    return dir_base.relativize(dir)
        .collect{ it.name }
        .join("_")
}

// ** Getting data for the -profile freesurfer ** //
workflow get_data_freesurfer {
    main:
        if (! params.input ) {
            log.info "You must provide an input folder containing all images required for FreesurferFlow :"
            log.info "        --input=/path/to/[input_folder]           Input folder containing your subjects."
            log.info "                              [input]"
            log.info "                               ├-- S1"
            log.info "                               |   └-- *t1.nii.gz"
            log.info "                               └-- S2"
            log.info "                                   └-- *t1.nii.gz"
            error "Please resubmit your command with the previous file structure."
        }

        input = file(params.input)

        // ** Loading files ** //
        anat_channel = Channel.fromPath("$input/**/*t1.nii.gz")
                        .map{ch1 -> [ch1.parent.name, ch1]}

    emit:
        anat = anat_channel
}

// ** Decided to split the data fetching steps for different profiles in different functions ** //
// ** for easier code-reading. ** //
workflow get_data_tracking {
    main:
        if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "        --input=/path/to/[input_folder]             Input folder containing multiple subjects for tracking"
            log.info ""
            log.info "                               [Input]"
            log.info "                               ├-- S1"
            log.info "                               |   ├-- *dwi.nii.gz"
            log.info "                               |   ├-- *dwi.bval"
            log.info "                               |   ├-- *dwi.bvec"
            log.info "                               |   ├-- *revb0.nii.gz"
            log.info "                               |   └-- *t1.nii.gz"
            log.info "                               └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"
            log.info "                                    ├-- *bval"
            log.info "                                    ├-- *bvec"
            log.info "                                    ├-- *revb0.nii.gz"
            log.info "                                    └-- *t1.nii.gz"
            error "Please resubmit your command with the previous file structure."
        }
        
        input = file(params.input)

        // ** Loading all files. ** //
        dwi_channel = Channel.fromFilePairs("$input/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true)
            { fetch_id(it.parent, input) }
        rev_channel = Channel.fromFilePairs("$input/**/*revb0.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        anat_channel = Channel.fromFilePairs("$input/**/*t1.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }

        // ** Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility. ** //
        dwi_channel = dwi_channel.map{sid, bvals, bvecs, dwi -> tuple(sid, dwi, bvals, bvecs)}
        
    emit:
        dwi = dwi_channel
        rev = rev_channel
        anat = anat_channel
}

// ** Getting data for the -profile tracking,infant ** //
workflow get_data_tracking_infant {
    main:
        if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "        --input=/path/to/[input_folder]             Input folder containing multiple subjects for tracking"
            log.info ""
            log.info "                               [Input]"
            log.info "                               ├-- S1"
            log.info "                               |   ├-- *dwi.nii.gz"
            log.info "                               |   ├-- *dwi.bval"
            log.info "                               |   ├-- *dwi.bvec"
            log.info "                               |   ├-- *revb0.nii.gz"
            log.info "                               |   ├-- *t2w.nii.gz"
            log.info "                               |   └-- *wm_mask.nii.gz"
            log.info "                               └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"
            log.info "                                    ├-- *bval"
            log.info "                                    ├-- *bvec"
            log.info "                                    ├-- *revb0.nii.gz"
            log.info "                                    ├-- *t2w.nii.gz"
            log.info "                                    └-- *wm_mask.nii.gz"
            error "Please resubmit your command with the previous file structure."
        }
        
        input = file(params.input)

        // ** Loading all files. ** //
        dwi_channel = Channel.fromFilePairs("$input/**/*dwi.{nii.gz,bval,bvec}", size: 3, flat: true)
            { fetch_id(it.parent, input) }
        rev_channel = Channel.fromFilePairs("$input/**/*revb0.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        anat_channel = Channel.fromFilePairs("$input/**/*t2w.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        wm_mask_channel = Channel.fromFilePairs("$input/**/*wm_mask.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }

        // ** Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility. ** //
        dwi_channel = dwi_channel.map{sid, bvals, bvecs, dwi -> tuple(sid, dwi, bvals, bvecs)}
        
    emit:
        dwi = dwi_channel
        rev = rev_channel
        anat = anat_channel
        wm_mask = wm_mask_channel
}

// ** Fetching data for -profile connectomics ** //
workflow get_data_connectomics {
    main:
        if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "     --input=/path/to/[input_folder]             Input folder containing multiple subjects"
            log.info ""
            log.info "                                [Input]"
            log.info "                                ├-- S1"
            log.info "                                |   ├-- *dwi.nii.gz"            
            log.info "                                |   ├-- *dwi.bval"            
            log.info "                                |   ├-- *dwi.bvec"                
            log.info "                                |   ├-- *t1.nii.gz"             
            log.info "                                |   ├-- *.trk"                  
            log.info "                                |   ├-- *labels.nii.gz"          
            log.info "                                |   ├-- *peaks.nii.gz"          
            log.info "                                |   ├-- *fodf.nii.gz"            
            log.info "                                |   ├-- OGenericAffine.mat"     
            log.info "                                |   ├-- output1Warp.nii.gz"  
            log.info "                                |   └-- metrics"
            log.info "                                |       └-- METRIC_NAME.nii.gz  [Optional]"
            log.info "                                └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"         
            log.info "                                    ├-- *bval"                  
            log.info "                                    ├-- *bvec"                  
            log.info "                                    ├-- *t1.nii.gz"           
            log.info "                                    ├-- *.trk"                   
            log.info "                                    ├-- *labels.nii.gz"         
            log.info "                                    ├-- *peaks.nii.gz"       
            log.info "                                    ├-- *fodf.nii.gz"         
            log.info "                                    ├-- OGenericAffine.mat"   
            log.info "                                    ├-- output1Warp.nii.gz"  
            log.info "                                    └-- metrics"
            log.info "                                        └-- METRIC_NAME.nii.gz  [Optional]"
            error "Please resubmit your command with the previous file structure."
        }

        input = file(params.input)

        // Loading all files.
        tracking_channel = Channel.fromFilePairs("$input/**/*.trk", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        labels_channel = Channel.fromFilePairs("$input/**/*labels.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        dwi_peaks_channel = Channel.fromFilePairs("$input/**/{*dwi.nii.gz,*.bval,*.bvec,*peaks.nii.gz}", size: 4, flat: true)
            { fetch_id(it.parent, input) }
        fodf_channel = Channel.fromFilePairs("$input/**/*fodf.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        metrics_channel = Channel.fromFilePairs("$input/**/metrics/*.nii.gz", size: -1, maxDepth: 2)
            { it.parent.parent.name }
        t1_channel = Channel.fromFilePairs("$input/**/*t1.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        transfos_channel = Channel.fromFilePairs("$input/**/{0GenericAffine.mat,output1Warp.nii.gz}", size: 2, flat: true)
            { fetch_id(it.parent, input) }

        // Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility.
        dwi_peaks_channel = dwi_peaks_channel.map{sid, bvals, bvecs, dwi, peaks -> tuple(sid, dwi, bvals, bvecs, peaks)}

        emit:
            trk = tracking_channel
            labels = labels_channel
            dwi_peaks = dwi_peaks_channel
            fodf = fodf_channel
            metrics = metrics_channel
            anat = t1_channel
            transfos = transfos_channel
}

// ** Fetching data for -profile connectomics,infant ** //
workflow get_data_connectomics_infant {
    main:
        if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "     --input=/path/to/[input_folder]             Input folder containing multiple subjects"
            log.info ""
            log.info "                                [Input]"
            log.info "                                ├-- S1"
            log.info "                                |   ├-- *dwi.nii.gz"            
            log.info "                                |   ├-- *dwi.bval"            
            log.info "                                |   ├-- *dwi.bvec"                
            log.info "                                |   ├-- *t2w.nii.gz"             
            log.info "                                |   ├-- *.trk"                  
            log.info "                                |   ├-- *labels.nii.gz"          
            log.info "                                |   ├-- *peaks.nii.gz"          
            log.info "                                |   ├-- *fodf.nii.gz"            
            log.info "                                |   ├-- OGenericAffine.mat"     
            log.info "                                |   ├-- output1Warp.nii.gz"  
            log.info "                                |   └-- metrics"
            log.info "                                |       └-- METRIC_NAME.nii.gz  [Optional]"
            log.info "                                └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"         
            log.info "                                    ├-- *bval"                  
            log.info "                                    ├-- *bvec"                  
            log.info "                                    ├-- *t2w.nii.gz"           
            log.info "                                    ├-- *.trk"                   
            log.info "                                    ├-- *labels.nii.gz"         
            log.info "                                    ├-- *peaks.nii.gz"       
            log.info "                                    ├-- *fodf.nii.gz"         
            log.info "                                    ├-- OGenericAffine.mat"   
            log.info "                                    ├-- output1Warp.nii.gz"  
            log.info "                                    └-- metrics"
            log.info "                                        └-- METRIC_NAME.nii.gz  [Optional]"
            error "Please resubmit your command with the previous file structure."
        }

        input = file(params.input)

        // Loading all files.
        tracking_channel = Channel.fromFilePairs("$input/**/*.trk", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        labels_channel = Channel.fromFilePairs("$input/**/*labels.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        dwi_peaks_channel = Channel.fromFilePairs("$input/**/{*dwi.nii.gz,*.bval,*.bvec,*peaks.nii.gz}", size: 4, flat: true)
            { fetch_id(it.parent, input) }
        fodf_channel = Channel.fromFilePairs("$input/**/*fodf.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        metrics_channel = Channel.fromFilePairs("$input/**/metrics/*.nii.gz", size: -1, maxDepth: 2)
            { it.parent.parent.name }
        t2w_channel = Channel.fromFilePairs("$input/**/*t2w.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        transfos_channel = Channel.fromFilePairs("$input/**/{0GenericAffine.mat,output1Warp.nii.gz}", size: 2, flat: true)
            { fetch_id(it.parent, input) }

        // Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility.
        dwi_peaks_channel = dwi_peaks_channel.map{sid, bvals, bvecs, dwi, peaks -> tuple(sid, dwi, bvals, bvecs, peaks)}

        emit:
            trk = tracking_channel
            labels = labels_channel
            dwi_peaks = dwi_peaks_channel
            fodf = fodf_channel
            metrics = metrics_channel
            anat = t2w_channel
            transfos = transfos_channel
}

workflow get_data_template {
    main:
    if ( !params.input ) {
            log.info "You must provide an input folder containing all images using:"
            log.info "        --input=/path/to/[input_folder]             Input folder containing multiple subjects for tracking"
            log.info ""
            log.info "                               [Input]"
            log.info "                               ├-- S1"
            log.info "                               |   ├-- *dwi.nii.gz"
            log.info "                               |   ├-- *dwi.bvec"
            log.info "                               |   ├-- *fa.nii.gz"
            log.info "                               |   ├-- *t2w.nii.gz"
            log.info "                               └-- S2"
            log.info "                                    ├-- *dwi.nii.gz"
            log.info "                                    ├-- *dwi.bvec"
            log.info "                                    ├-- *fa.nii.gz"
            log.info "                                    └-- *t2w.nii.gz"
            log.info "                               [References]"
            log.info "                               ├-- *fa_ref.nii.gz"
            log.info "                               └-- *t2w_ref.nii.gz"
            error "Please resubmit your command with the previous file structure."
        }
        
        input = file(params.input)
        references = file(params.references)

        // Loading all files.
        dwi_channel = Channel.fromFilePairs("$input/**/*dwi.{nii.gz,bvec}", size: 2, flat: true)
            { fetch_id(it.parent, input) }
        fa_channel = Channel.fromFilePairs("$input/**/*fa.nii.gz", size:1, flat: true)
            { fetch_id(it.parent, input) }
        anat_channel = Channel.fromFilePairs("$input/**/*t2w.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        anat_ref = Channel.fromPath("$references/*t2w_ref.nii.gz")
        fa_ref = Channel.fromPath("$references/*fa_ref.nii.gz")

        // Setting up dwi channel in this order : sid, dwi, bval, bvec for lisibility.
        dwi_channel = dwi_channel.map{sid, bvecs, dwi -> tuple(sid, dwi, bvecs)}
        
    emit:
        dwi = dwi_channel
        anat = anat_channel
        fa = fa_channel
        anat_ref = anat_ref
        fa_ref = fa_ref
}

def display_run_info () {
    log.info ""
    log.info "ChildBrainFlow pipeline"
    log.info "========================"
    log.info "ChildBrainFlow is an end-to-end pipeline that performs tractography, t1 reconstruction and connectomics. " 
    log.info "It is essentially a merged version of multiple individual pipeline to avoid the handling of inputs/outputs "
    log.info "between flows with some parameters tuned for pediatric brain scans. "
    log.info ""
    log.info "Start time: $workflow.start"
    log.info ""

    log.debug "[Command-line]"
    log.debug "$workflow.commandLine"
    log.debug ""

    log.info "[Git Info]"
    log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
    log.info ""

    log.info "[Inputs]"
    log.info "Input: $params.input"
    log.info "Output Directory: $params.output_dir"
    log.info ""

    if ( params.run_tracking ) {
        log.info "[Tracking Options]"
        log.info ""
        log.info "GLOBAL OPTIONS"
        log.info "Threshold for b0: $params.b0_thr"
        log.info "DWI Shell Tolerance: $params.dwi_shell_tolerance"
        log.info "Skip DWI preprocessing: $params.skip_dwi_preprocessing"
        log.info ""
        log.info "BET DWI OPTIONS"
        log.info "Initial fractional value for BET: $params.initial_bet_f"
        log.info "Finale fractional value for BET: $params.final_bet_f"
        log.info ""
        if ( params.infant_config ) {
            log.info "BET T2W OPTIONS"
            log.info "Run BET on T2W image: $params.run_bet_anat"
            log.info "Fractional value for T2W BET: $params.bet_anat_f"
        }
        else {
            log.info "BET T1W OPTIONS"
            log.info "T1 Tempalte: $params.template_t1"
        }
        log.info ""
        log.info "EDDY AND TOPUP OPTIONS"
        log.info "Configuration for topup: $params.topup_config"
        log.info "Encoding direction: $params.encoding_direction"
        log.info "Readout: $params.readout"
        log.info "Topup prefix: $params.topup_prefix"
        log.info "Topup BET fractional value: $params.topup_bet_f"
        log.info "Eddy command: $params.eddy_cmd"
        log.info "Run slice drop correction: $params.use_slice_drop_correction"
        log.info ""
        log.info "NORMALIZE OPTIONS"
        log.info "FA threshold for masking: $params.fa_mask_threshold"
        log.info ""
        log.info "RESAMPLE ANAT OPTIONS"
        log.info "Resampling resolution for Anatomical file: $params.anat_resolution"
        log.info "Interpolation method for Anatomical file: $params.anat_interpolation"
        log.info "Interpolation method for masks: $params.mask_interpolation"
        log.info ""
        log.info "RESAMPLE DWI OPTIONS"
        log.info "Resampling resolution for DWI: $params.dwi_resolution"
        log.info "Interpolation method for DWI: $params.dwi_interpolation"
        log.info "Interpolation method for DWI mask: $params.mask_dwi_interpolation"
        log.info ""
        log.info "EXTRACT DWI SHELLS OPTIONS"
        if ( params.dti_shells ) {
            log.info "DTI Shells: $params.dti_shells"
        } else {
            log.info "Maximum DTI shell value: $params.max_dti_shell_value"
        }
        log.info ""
        log.info "SH FITTING OPTIONS"
        log.info "Run SH fitting: $params.sh_fitting"
        log.info "SH fitting order: $params.sh_fitting_order"
        log.info "SH fitting basis: $params.sh_fitting_basis"
        log.info ""
        log.info "FODF OPTIONS"
        if ( params.fodf_shells ) {
            log.info "FODF Shells: $params.fodf_shells"
        } else {
            log.info "Minimum fODF shell value: $params.min_fodf_shell_value"
        }
        log.info "FODF Metrics A factor: $params.fodf_metrics_a_factor"
        log.info "Maximum FA value in ventricles: $params.max_fa_in_ventricle"
        log.info "Minimum MD value in ventricles: $params.min_md_in_ventricle"
        log.info "Relative threshold (RT): $params.relative_threshold"
        log.info "SH basis: $params.basis"
        log.info "SH order: $params.sh_order"
        log.info ""
        log.info "FRF OPTIONS"
        log.info "Run mean FRF: $params.mean_frf"
        log.info "FA threshold for single fiber voxel: $params.fa"
        log.info "Minimum FA for selecting voxel: $params.min_fa"
        log.info "Minimum number of voxels: $params.min_nvox"
        log.info "ROI radius: $params.roi_radius"
        log.info "Set FRF: $params.set_frf"
        log.info "Manual FRF: $params.manual_frf"
        log.info ""
        log.info "SEGMENT TISSUES OPTIONS"
        log.info "Number of tissues: $params.number_of_tissues"
        log.info ""
        log.info "SEEDING AND TRACKING OPTIONS"
        if ( params.run_pft_tracking ) {
            log.info "Algorithm for tracking: $params.pft_algo"
            log.info "Number of seeds per voxel: $params.pft_nbr_seeds"
            log.info "Seeding method: $params.pft_seeding"
            log.info "Step size: $params.pft_step_size"
            log.info "Theta threshold: $params.pft_theta"
            log.info "Minimum fiber length: $params.pft_min_len"
            log.info "Maximum fiber length: $params.pft_max_len"
            log.info "Compression: $params.pft_compress_streamlines"
        }
        else {
            log.info "Algorithm for tracking: $params.local_algo"
            log.info "Number of seeds per voxel: $params.local_nbr_seeds"
            log.info "Seeding method: $params.local_seeding"
            log.info "Step size: $params.local_step_size"
            log.info "Theta threshold: $params.local_theta"
            log.info "Minimum fiber length: $params.local_min_len"
            log.info "Maximum fiber length: $params.local_max_len"
            log.info "Compression: $params.local_compress_streamlines"
        }
        log.info ""
        log.info "PROCESSES PER TASKS"
        if ( !params.infant_config ) {
            log.info "Processes for denoising T1: $params.processes_denoise_t1"
            log.info "Processes for BET T1: $params.processes_bet_t1"
        }
        log.info "Processes for denoising DWI: $params.processes_denoise_dwi"
        log.info "Processes for EDDY: $params.processes_eddy"
        log.info "Processes for registration: $params.processes_registration"
        log.info "Processes for FODF: $params.processes_fodf"
        log.info ""
    }

    if ( params.run_freesurfer ) {
        log.info "[Freesurfer Options]"
        log.info ""
        log.info "Use Recon-all: $params.recon_all"
        log.info "Use FastSurfer + Recon-Surf: $params.recon_surf"
        log.info "Atlas utils folder: $params.atlas_utils_folder"
        log.info "Compute FS, BN, GL, SF: $params.compute_FS_BN_GL_SF"
        log.info "Compute BN Child: $params.compute_BN_child"
        log.info "Compute lobes: $params.compute_lobes"
        log.info "Compute lausanne multiscale: $params.compute_lausanne_multiscale"
        log.info "Number of threads: $params.nb_threads"
        log.info ""
        log.info "ATLAS SELECTION"
        log.info "Use Freesurfer atlas: $params.use_freesurfer_atlas"
        log.info "Use Brainnetome atlas: $params.use_brainnetome_atlas"
        log.info "Use Brainnetome Child atlas: $params.use_brainnetome_child_atlas"
        log.info "Use Glasser atlas: $params.use_glasser_atlas"
        log.info "Use Schaefer 100 atlas: $params.use_schaefer_100_atlas"
        log.info "Use Schaefer 200 atlas: $params.use_schaefer_200_atlas"
        log.info "Use Schaefer 400 atlas: $params.use_schaefer_400_atlas"
        log.info "Use Lausanne 1 atlas: $params.use_lausanne_1_atlas"
        log.info "Use Lausanne 2 atlas: $params.use_lausanne_2_atlas"
        log.info "Use Lausanne 3 atlas: $params.use_lausanne_3_atlas"
        log.info "Use Lausanne 4 atlas: $params.use_lausanne_4_atlas"
        log.info "Use Lausanne 5 atlas: $params.use_lausanne_5_atlas"
        log.info "Use dilated labels: $params.use_dilated_labels"
        log.info ""
    }

    if ( params.run_connectomics ) {
        log.info "[Connectomics Options]"
        log.info ""
        log.info "DECOMPOSE OPTIONS"
        log.info "No pruning: $params.no_pruning"
        log.info "No remove loops: $params.no_remove_loops"
        log.info "No remove outliers: $params.no_remove_outliers"
        log.info "Minimal length: $params.min_length"
        log.info "Maximal length: $params.max_length"
        log.info "Maximum looping angle: $params.loop_max_angle"
        log.info "Outlier treshold: $params.outlier_threshold"
        log.info ""
        log.info "COMMIT OPTIONS"
        log.info "Run COMMIT: $params.run_commit"
        log.info "Use COMMIT2: $params.use_commit2"
        log.info "Use both COMMIT: $params.use_both_commit"
        log.info "COMMIT on trk: $params.commit_on_trk"
        log.info "B-value threshold: $params.b_thr"
        log.info "Number of directions: $params.nbr_dir"
        log.info "Ball and stick: $params.ball_stick"
        log.info "Parallel diffusivity: $params.para_diff"
        log.info "Perpendicular diffusivity: $params.perp_diff"
        log.info "Isotropic diffusivity: $params.iso_diff"
        log.info ""
        log.info "PROCESSES OPTIONS"
        log.info "Number of processes for COMMIT: $params.processes_commit"
        log.info "Number of processes for AFD_FIXEL: $params.processes_afd_fixel"
        log.info "Number of processes for CONNECTIVITY: $params.processes_connectivity"
        log.info "" 
    }
}