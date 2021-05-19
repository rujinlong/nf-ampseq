process TRAIN_CLASSIFIER {
    publishDir "$baseDir", mode: 'copy', overwrite: false

    output:
    path("classifier.qza"), emit: classifier_ch

    when:
    params.mode == "all"

    """
    qiime tools import \
        --type 'FeatureData[Taxonomy]' \
        --input-format HeaderlessTSVTaxonomyFormat \
        --input-path $params.SILVA_taxonomy \
        --output-path ref_taxonomy.qza

    qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path $params.SILVA_sequence \
        --output-path ref_seqs.qza

    qiime feature-classifier extract-reads \
        --i-sequences ref_seqs.qza \
        --p-f-primer $params.primer_forward \
        --p-r-primer $params.primer_reverse \
        --p-min-length $params.min_reference_length \
        --p-max-length $params.max_reference_length \
        --o-reads ref_amplicon_seqs.qza \
        --p-n-jobs $task.cpus \
        --verbose \
        &> z01_SILVA_99_otus_16S_V4_seqs.log

    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-taxonomy ref_taxonomy.qza \
        --i-reference-reads ref_amplicon_seqs.qza \
        --o-classifier classifier.qza \
        --verbose \
        &> classifier.log
    """
}