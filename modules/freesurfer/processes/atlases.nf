#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FS_BN_GL_SF {
    cpus params.nb_threads

    input:
        tuple val(sid), path(folder)

    output:
        tuple val(sid), path("*freesurfer_v5.nii.gz"), emit: freesurfer
        tuple val(sid), path("*freesurfer_v5_dilate.nii.gz"), emit: freesurfer_dilated
        tuple val(sid), path("*brainnetome_v5.nii.gz"), emit: brainnetome
        tuple val(sid), path("*brainnetome_v5_dilate.nii.gz"), emit: brainnetome_dilated
        tuple val(sid), path("*glasser_v5.nii.gz"), emit: glasser
        tuple val(sid), path("*glasser_v5_dilate.nii.gz"), emit: glasser_dilated
        tuple val(sid), path("*schaefer_100_v5.nii.gz"), emit: schaefer_100
        tuple val(sid), path("*schaefer_100_v5_dilate.nii.gz"), emit: schaefer_100_dilated
        tuple val(sid), path("*schaefer_200_v5.nii.gz"), emit: schaefer_200
        tuple val(sid), path("*schaefer_200_v5_dilate.nii.gz"), emit: schaefer_200_dilated
        tuple val(sid), path("*schaefer_400_v5.nii.gz"), emit: schaefer_400
        tuple val(sid), path("*schaefer_400_v5_dilate.nii.gz"), emit: schaefer_400_dilated
        path("*[brainnetome,freesurfer,glasser,schaefer]*.txt")
        path("*[brainnetome,freesurfer,glasser,schaefer]*.json")

    when:
        params.compute_FS_BN_GL_SF

    script:
    """
    ln -s $params.atlas_utils_folder/fsaverage \$(dirname ${folder})/
    bash $params.atlas_utils_folder/freesurfer_utils/generate_atlas_FS_BN_GL_SF_v5.sh \$(dirname ${folder}) \
        ${sid} ${params.nb_threads} FS_BN_GL_SF_Atlas/
    cp $sid/FS_BN_GL_SF_Atlas/* ./
    """
}

process BN_CHILD {
    cpus params.nb_threads

    input:
        tuple val(sid), path(folder)

    output:
        tuple val(sid), path("*brainnetome_child_v1.nii.gz"), emit: brainnetome_child
        tuple val(sid), path("*brainnetome_child_v1_dilate.nii.gz"), emit: brainnetome_child_dilated
        path("*[brainnetome_child]*.txt")
        path("*[brainnetome_child]*.json")

    when:
        params.compute_BN_child

    script:
    """
    ln -s $params.atlas_utils_folder/fsaverage \$(dirname ${folder})/
    bash $params.atlas_utils_folder/freesurfer_utils/generate_atlas_BN_child.sh \$(dirname ${folder}) \
        ${sid} ${params.nb_threads} Child_Atlas/
    cp $sid/Child_Atlas/* ./
    """
}

process LOBES {
    cpus params.nb_threads

    input:
        tuple val(sid), path(folder)

    output: 
        path("*lobes*.nii.gz"), emit: lobes
        path("*lobes*.txt")
        path("*lobes*.json")

    when:
        params.compute_lobes
    
    script:
    """
    mri_convert ${folder}/mri/rawavg.mgz rawavg.nii.gz

    mri_convert ${folder}/mri/wmparc.mgz wmparc.nii.gz
    scil_reshape_to_reference.py wmparc.nii.gz rawavg.nii.gz wmparc.nii.gz --interpolation nearest -f 
    scil_image_math.py convert wmparc.nii.gz wmparc.nii.gz --data_type uint16 -f
    
    mri_convert ${folder}/mri/brainmask.mgz brain_mask.nii.gz
    scil_image_math.py lower_threshold brain_mask.nii.gz 0.001 brain_mask.nii.gz --data_type uint8 -f
    scil_image_math.py dilation brain_mask.nii.gz 1 brain_mask.nii.gz -f
    scil_reshape_to_reference.py brain_mask.nii.gz rawavg.nii.gz brain_mask.nii.gz --interpolation nearest -f 
    scil_image_math.py convert brain_mask.nii.gz brain_mask.nii.gz --data_type uint8 -f

    scil_combine_labels.py atlas_lobes_v5.nii.gz --volume_ids wmparc.nii.gz 1003 1012 1014 1017 1018 1019 1020 \
        1024 1027 1028 1032 --volume_ids wmparc.nii.gz 1008 1022 1025 1029 1031 --volume_ids wmparc.nii.gz 1005 \
        1011 1013 1021 --volume_ids wmparc.nii.gz 1001 1006 1007 1009 1015 1015 1030 1033 --volume_ids wmparc.nii.gz \
        1002 1010 1023 1026 --volume_ids wmparc.nii.gz 8 --volume_ids wmparc.nii.gz 10 11 12 13 17 18 26 28 \
        --volume_ids wmparc.nii.gz 2003 2012 2014 2017 2018 2019 2020 2024 2027 2028 2032 --volume_ids wmparc.nii.gz \
        2008 2022 2025 2029 2031 --volume_ids wmparc.nii.gz 2005 2011 2013 2021 --volume_ids wmparc.nii.gz 2001 2006 \
        2007 2009 2015 2015 2030 2033 --volume_ids wmparc.nii.gz 2002 2010 2023 2026 --volume_ids wmparc.nii.gz 49 50 \
        51 52 53 54 58 60 --volume_ids wmparc.nii.gz 47 --volume_ids wmparc.nii.gz 16 --merge
    scil_dilate_labels.py atlas_lobes_v5.nii.gz atlas_lobes_v5_dilate.nii.gz --distance 2 \
        --labels_to_dilate 1 2 3 4 5 6 8 9 10 11 12 14 15 --mask brain_mask.nii.gz
    cp $params.atlas_utils_folder/freesurfer_utils/*lobes_v5* ./
    """
}

process LAUSANNE {
    cpus 1

    input:
        tuple val(sid), path(folder)
        each scale

    output:
        tuple val(sid), path("lausanne_2008_scale_1*.nii.gz"), emit: lausanne_1, optional: true
        tuple val(sid), path("lausanne_2008_scale_2*.nii.gz"), emit: lausanne_2, optional: true
        tuple val(sid), path("lausanne_2008_scale_3*.nii.gz"), emit: lausanne_3, optional: true
        tuple val(sid), path("lausanne_2008_scale_4*.nii.gz"), emit: lausanne_4, optional: true
        tuple val(sid), path("lausanne_2008_scale_5*.nii.gz"), emit: lausanne_5, optional: true
        path("*.txt")
        path("*.json")
    
    when:
        params.compute_lausanne_multiscale

    script:
    """
    ln -s $params.atlas_utils_folder/fsaverage \$(dirname ${folder})/
    freesurfer_home=\$(dirname \$(dirname \$(which mri_label2vol)))
    /usr/bin/python $params.atlas_utils_folder/lausanne_multi_scale_atlas/generate_multiscale_parcellation.py \
        \$(dirname ${folder}) ${sid} \$freesurfer_home --scale ${scale} --dilation_factor 0 --log_level DEBUG

    mri_convert ${folder}/mri/rawavg.mgz rawavg.nii.gz
    scil_image_math.py lower_threshold rawavg.nii.gz 0.001 mask.nii.gz --data_type uint8
    scil_reshape_to_reference.py ${folder}/mri/lausanne2008.scale${scale}+aseg.nii.gz mask.nii.gz \
        lausanne_2008_scale_${scale}.nii.gz --interpolation nearest
    scil_image_math.py convert lausanne_2008_scale_${scale}.nii.gz lausanne_2008_scale_${scale}.nii.gz \
        --data_type int16 -f
    scil_dilate_labels.py lausanne_2008_scale_${scale}.nii.gz lausanne_2008_scale_${scale}_dilate.nii.gz \
        --distance 2 --mask mask.nii.gz
    
    cp $params.atlas_utils_folder/lausanne_multi_scale_atlas/*.txt ./
    cp $params.atlas_utils_folder/lausanne_multi_scale_atlas/*.json ./
    """
}
