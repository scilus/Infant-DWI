#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process FODF_SHELL {
    cpus 3

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec)
    output:
        tuple val(sid), path("${sid}__dwi_fodf.nii.gz"), path("${sid}__bval_fodf"),
        path("${sid}__bvec_fodf"), emit: dwi_fodf
    script:
    if (params.fodf_shells)
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_extract_dwi_shell.py $dwi \
            $bval $bvec $params.fodf_shells ${sid}__dwi_fodf.nii.gz \
            ${sid}__bval_fodf ${sid}__bvec_fodf -t $params.dwi_shell_tolerance -f
    """
    else
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1

        shells=\$(awk -v min_fodf="$params.min_fodf_shell_value" -v b0_thr="$params.b0_thr" '{for (i = 1; i <= NF; i++) 
                {v = int(\$i);if (v >= min_fodf || v <= b0_thr) shells[v] = 1;}}
                END {
                    for (v in shells) print v;
                    }
                ' "$bval" | sort -n | tr '\n' ' ')  

        scil_extract_dwi_shell.py $dwi \
        $bval $bvec \$shells ${sid}__dwi_fodf.nii.gz \
        ${sid}__bval_fodf ${sid}__bvec_fodf -t $params.dwi_shell_tolerance -f
    """
}

process COMPUTE_FRF {
    cpus 3

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(b0_mask)
    output:
        tuple val(sid), path("${sid}__frf.txt"), emit: frf
    script:
    if (params.set_frf)
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_ssst_frf.py $dwi $bval $bvec frf.txt --mask $b0_mask\
        --fa $params.fa --min_fa $params.min_fa --min_nvox $params.min_nvox\
        --roi_radii $params.roi_radius --force_b0_threshold
        scil_set_response_function.py frf.txt $params.manual_frf ${sid}__frf.txt
        """
    else
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_ssst_frf.py $dwi $bval $bvec ${sid}__frf.txt --mask $b0_mask\
        --fa $params.fa --min_fa $params.min_fa --min_nvox $params.min_nvox\
        --roi_radii $params.roi_radius --force_b0_threshold
        """
}

process MEAN_FRF {
    cpus 1
    publishDir = "${params.output_dir}/MEAN_FRF"

    input:
        path(all_frf)
    output:
        path("mean_frf.txt"), emit: mean_frf
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_mean_frf.py $all_frf mean_frf.txt
    """
}

process FODF_METRICS {
    cpus params.processes_fodf

    input:
        tuple val(sid), path(dwi), path(bval), path(bvec), path(b0_mask), path(fa), path(md), path(frf)
    output:
        tuple val(sid), path("${sid}__fodf.nii.gz"), emit: fodf
        tuple val(sid), path("${sid}__peaks.nii.gz"), emit: peaks
        tuple val(sid), 
        path("${sid}__afd_total.nii.gz"),
        path("${sid}__nufo.nii.gz"), emit: afd_and_nufo
        tuple val(sid), 
        path("${sid}__peak_indices.nii.gz"),
        path("${sid}__afd_max.nii.gz"),
        path("${sid}__afd_sum.nii.gz"),
        path("${sid}__peak_values.nii.gz"),
        path("${sid}__rgb.nii.gz"), emit: fodf_metrics
    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_ssst_fodf.py $dwi $bval $bvec $frf ${sid}__fodf.nii.gz\
        --sh_order $params.sh_order --sh_basis $params.basis --force_b0_threshold\
        --mask $b0_mask --processes $task.cpus
    scil_compute_fodf_max_in_ventricles.py ${sid}__fodf.nii.gz $fa $md\
        --max_value_output ventricles_fodf_max_value.txt --sh_basis $params.basis\
        --fa_t $params.max_fa_in_ventricle --md_t $params.min_md_in_ventricle\
        -f
    
    a_threshold=\$( echo " $params.fodf_metrics_a_factor * `awk '{for(i=1;i<=NF;i++) if(\$i>maxval) maxval=\$i;}; END { print maxval;}' ventricles_fodf_max_value.txt`" | bc )

    scil_compute_fodf_metrics.py ${sid}__fodf.nii.gz --mask $b0_mask --sh_basis $params.basis\
        --peaks ${sid}__peaks.nii.gz --peak_indices ${sid}__peak_indices.nii.gz\
        --afd_max ${sid}__afd_max.nii.gz --afd_total ${sid}__afd_total.nii.gz\
        --afd_sum ${sid}__afd_sum.nii.gz --nufo ${sid}__nufo.nii.gz\
        --peak_values ${sid}__peak_values.nii.gz --rgb ${sid}__rgb.nii.gz\
        --rt $params.relative_threshold --at \$a_threshold
    """
}