#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process AVERAGE_VOLUMES_ANAT {
    label "AVERAGE_VOLUMES_ANAT"
    cpus 4
    publishDir = params.Pop_Avg_Publish_Dir

    input:
        path(volumes)
    output:
        tuple path("population_avg_anat.nii.gz"), 
        path("population_avg_anat_bet.nii.gz"), 
        path("population_avg_anat_bet_mask.nii.gz"), emit: popavg
    script:
    """
    scil_image_math.py mean $volumes population_avg_anat.nii.gz
    bet population_avg_anat.nii.gz temp.nii.gz -f 0.7 -R
    bet temp.nii.gz population_avg_anat_bet -m -R
    """
}

process AVERAGE_DWI {
    label "AVERAGE_DWI"
    cpus 4
    publishDir = params.Pop_Avg_Publish_Dir

    input:
        path(volumes)
    output:
        path("population_avg_dwi.nii.gz"), emit: dwipopavg

    script:
    """
    mrmath $volumes mean population_avg_dwi.nii.gz
    """
}