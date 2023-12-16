#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FREESURFER {
    cpus params.nb_threads

    input:
        tuple val(sid), path(anat)
    output:
        tuple val(sid), path("$sid/"), emit: folders
        tuple val(sid), path("${sid}__final_t1.nii.gz"), emit: final_t1

    script:
    """
    export SUBJECTS_DIR=.
    recon-all -i $anat -s $sid -all -parallel -openmp $params.nb_threads
    mri_convert $sid/mri/antsdn.brain.mgz ${sid}__final_t1.nii.gz
    """
}

process RECON_SURF {
    cpus params.nb_threads

    input:
        tuple val(sid), path(anat)
    output:
        tuple val(sid), path("$sid/"), emit: folders
        tuple val(sid), path("${sid}__final_t1.nii.gz"), emit: final_t1

    script:
    // ** Adding a registration to .gca atlas to generate the talairach.m3z file (subcortical atlas segmentation ** //
    // ** wont work without it). A little time consuming but necessary. For FreeSurfer 7.3.2, RB_all_2020-01-02.gca ** //
    // ** is the default atlas. Update when bumping FreeSurfer version. ** //
    """
    mkdir output/
    bash /FastSurfer/run_fastsurfer.sh --sd \$(readlink -f ./) --sid $sid \\
        --t1 \$(readlink -f $anat) \
        --fs_license /freesurfer/license.txt \
        --parallel --device cpu --threads $params.nb_threads --allow_root
    mri_ca_register -align-after -nobigventricles -mask $sid/mri/brainmask.mgz \
        -T $sid/mri/transforms/talairach.lta -threads $params.nb_threads $sid/mri/norm.mgz \
        \${FREESURFER_HOME}/average/RB_all_2020-01-02.gca $sid/mri/transforms/talairach.m3z
    mri_convert $sid/mri/antsdn.brain.mgz ${sid}__final_t1.nii.gz
    """
}