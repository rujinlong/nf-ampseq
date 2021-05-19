process TAXONOMY {
    publishDir "$params.outdir/p03_taxonomy"
    publishDir "$params.report", pattern: "taxonomy_barplots"
    publishDir "$params.report", pattern: "taxonomy_rep_seqs"
    publishDir "$params.report", pattern: "taxa_abundance.tsv"

    input:
    path(classifier)
    path(repseq)
    path(table)
    path(metadata)
    path(feature_table)

    output:
    path("taxa_abundance.tsv")
    path("taxonomy_rep_seqs")
    path("taxonomy_barplots")
    path("taxonomy_rep_seqs.qza"), emit: taxa_repseq

    when:
    params.mode == "all"

    """
    qiime feature-classifier classify-sklearn \
        --i-classifier $classifier \
        --i-reads $repseq \
        --p-n-jobs $task.cpus \
        --o-classification taxonomy_rep_seqs.qza

    qiime metadata tabulate \
        --m-input-file taxonomy_rep_seqs.qza \
        --o-visualization taxonomy_rep_seqs.qzv

    qiime taxa barplot \
        --i-table $table \
        --i-taxonomy taxonomy_rep_seqs.qza \
        --m-metadata-file $metadata \
        --o-visualization taxonomy_barplots.qzv

    qiime tools export --input-path taxonomy_barplots.qzv --output-path taxonomy_barplots
    qiime tools export --input-path taxonomy_rep_seqs.qzv --output-path taxonomy_rep_seqs
    create_taxa_abundance_table.py -f $feature_table -t taxonomy_rep_seqs/metadata.tsv -o taxa_abundance.tsv
    """
}