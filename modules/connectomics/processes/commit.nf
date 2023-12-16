#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process COMMIT {
    cpus params.processes_commit
    memory params.commit_memory_limit

    input:
        tuple val(sid), path(h5), path(dwi), path(bval), path(bvec), path(peaks), path(para_diff), path(iso_diff), path(perp_diff)
    output:
        tuple val(sid), path("${sid}__decompose_commit.h5"), emit: h5_commit, optional: true
        tuple val(sid), path("${sid}__essential_tractogram.trk"), emit: trk_commit, optional: true
        tuple val(sid), path("${sid}__results_bzs/"), optional: true
        tuple val(sid), path("${sid}__results_bzs_1/"), optional: true
        tuple val(sid), path("${sid}__results_bzs_2/"), optional: true
    when:
        params.run_commit
    
    script:
    def para_diff_arg = para_diff ? "--para_diff \$(cat $para_diff)" : "--para_diff $params.para_diff"
    def iso_diff_arg = iso_diff ? "--iso_diff \$(cat $iso_diff)" : "--iso_diff $params.iso_diff"
    def perp_diff_arg = perp_diff ? "--perp_diff \$(cat $perp_diff)" : "--perp_diff $params.perp_diff"
    def ball_stick_arg = params.ball_stick ? "--ball_stick" : ""

    if ( params.use_commit2 && !params.use_both ) {
    """
    scil_run_commit.py $h5 $dwi $bval $bvec "${sid}__results_bzs/" --ball_stick --commit2 \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir\
        $para_diff_arg $iso_diff_arg
    mv "${sid}__results_bzs/commit_2/decompose_commit.h5" "./${sid}__decompose_commit.h5"
    """
    }
    else if ( params.use_both ) {
    """
    scil_run_commit.py $h5 $dwi $bval $bvec "${sid}__results_bzs_1/" --ball_stick --commit2 \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir\
        $para_diff_arg $iso_diff_arg
    scil_run_commit.py ${sid}__results_bzs_1/commit_2/essential_tractogram.trk $dwi $bval $bvec "${sid}__results_bzs_2/"\
        --in_peaks $peaks --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir\
        $para_diff_arg $iso_diff_arg $perp_diff_arg
    mv "${sid}__results_bzs_2/commit_1/essential_tractogram.trk" "./${sid}__essential_tractogram.trk"
    """
    }
    else {
    """
    scil_run_commit.py $h5 $dwi $bval $bvec "${sid}__results_bzs/" --in_peaks $peaks \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir $ball_stick_arg \
        $para_diff_arg $iso_diff_arg $perp_diff_arg
    mv "${sid}__results_bzs/commit_1/decompose_commit.h5" "./${sid}__decompose_commit.h5"
    """
    }
}

process COMMIT_ON_TRK {
    cpus params.processes_commit
    memory params.commit_memory_limit

    input:
        tuple val(sid), path(trk_h5), path(dwi), path(bval), path(bvec), path(peaks)
    output:
        tuple val(sid), path("${sid}__essential_tractogram.trk"), emit: trk_commit
        tuple val(sid), path("${sid}__results_bzs/")
    when:
        params.run_commit

    script:
    ball_stick_arg=""
    perp_diff_arg=""
    if ( params.ball_stick ) {
        ball_stick_arg="--ball_stick"
    }
    else {
        perp_diff_arg="--perp_diff $params.perp_diff"
    }
    """
    scil_run_commit.py $trk_h5 $dwi $bval $bvec "${sid}__results_bzs/" --in_peaks $peaks \
        --processes $params.processes_commit --b_thr $params.b_thr --nbr_dir $params.nbr_dir $ball_stick_arg \
        --para_diff $params.para_diff $perp_diff_arg --iso_diff $params.iso_diff
    mv "${sid}__results_bzs/commit_1/essential_tractogram.trk" "./${sid}__essential_tractogram.trk"
    """
}

process COMPUTE_PRIORS {
    cpus 1

    input:
        tuple val(sid), path(fa), path(md), path(ad), path(rd)
    output:
        tuple val("Priors"), path("${sid}__para_diff.txt"), emit: para_diff
        tuple val("Priors"), path("${sid}__iso_diff.txt"), emit: iso_diff
        tuple val("Priors"), path("${sid}__perp_diff.txt"), emit: perp_diff
        tuple val(sid), path("${sid}__mask_1fiber.nii.gz"), emit: mask_1fiber
        tuple val(sid), path("${sid}__mask_ventricles.nii.gz"), emit: mask_ventricles

    when:
        params.run_commit && params.compute_priors

    script:
    """
    scil_compute_NODDI_priors.py $fa $ad $rd $md \
        --out_txt_1fiber_para ${sid}__para_diff.txt \
        --out_txt_ventricles ${sid}__iso_diff.txt \
        --out_txt_1fiber_perp ${sid}__perp_diff.txt \
        --out_mask_1fiber ${sid}__mask_1fiber.nii.gz \
        --out_mask_ventricles ${sid}__mask_ventricles.nii.gz \
        --fa_min $params.fa_min_priors --fa_max $params.fa_max_priors \
        --md_min $params.md_min_priors --roi_radius $params.roi_radius_priors
    """
}

process AVERAGE_PRIORS {
    cpus 1

    input:
        tuple val(sid), path(para_diff), path(iso_diff), path(perp_diff)
    
    output:
        path("mean_para_diff.txt"), emit: mean_para_diff
        path("mean_iso_diff.txt"), emit: mean_iso_diff
        path("mean_perp_diff.txt"), emit: mean_perp_diff
    
    script:
    """
    cat $para_diff > all_para_diff.txt
    awk '{ total += \$1; count++ } END { print total/count }' all_para_diff.txt > mean_para_diff.txt
    cat $iso_diff > all_iso_diff.txt
    awk '{ total += \$1; count++ } END { print total/count }' all_iso_diff.txt > mean_iso_diff.txt
    cat $perp_diff > all_perp_diff.txt
    awk '{ total += \$1; count++ } END { print total/count }' all_perp_diff.txt > mean_perp_diff.txt
    """
}