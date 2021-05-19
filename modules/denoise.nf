process DENOISE {
    label "big"
    publishDir "$params.outdir/p02_denoise"
    publishDir "$params.report", pattern: "denoising_stats*"
    publishDir "$params.report", pattern: "representative_sequences*"
    publishDir "$params.report", pattern: "table*"
    publishDir "$params.report", pattern: "feature_table.tsv"

    input:
    path(clean_reads)

    output:
    path("denoising_stats")
    path("table")
    path("representative_sequences")
    path("repseqs.fasta"), emit: repseq
    path("representative_sequences.qza"), emit: repseq_ch
    path("table.qza"), emit: table_ch
    path("feature_table.tsv"), emit: table2taxa_abundance
    

    when:
    params.mode == "all"

    """
    qiime dada2 denoise-paired \
        --p-n-threads $task.cpus \
        --i-demultiplexed-seqs $clean_reads \
        --p-trunc-len-f $params.trunc_len_forward \
        --p-trunc-len-r $params.trunc_len_reverse \
        --p-max-ee-f 3 \
        --p-max-ee-r 3 \
        --p-min-fold-parent-over-abundance $params.mfpoa \
        --o-table table.qza \
        --o-representative-sequences representative_sequences.qza \
        --o-denoising-stats denoising_stats.qza \
        --verbose &> log_denoise.log

    # Denoising stats
    qiime metadata tabulate \
        --m-input-file denoising_stats.qza \
        --o-visualization denoising_stats.qzv

    # Representative sequences
    qiime feature-table tabulate-seqs \
        --i-data representative_sequences.qza \
        --o-visualization representative_sequences.qzv

    # Feature table
    qiime feature-table summarize \
        --i-table table.qza \
        --o-visualization table.qzv
    
    qiime tools export --input-path denoising_stats.qzv --output-path denoising_stats
    qiime tools export --input-path representative_sequences.qzv --output-path representative_sequences
    ln -s representative_sequences/sequences.fasta repseqs.fasta
    qiime tools export --input-path table.qza --output-path table
    qiime tools export --input-path table.qzv --output-path table
    biom convert -i table/feature-table.biom -o feature_table.tsv --to-tsv
    """
}