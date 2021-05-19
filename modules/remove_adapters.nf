process REMOVE_ADAPTERS {
    label "small"
    publishDir "$params.outdir/p01_clean_reads"

    input:
    path(manifest)

    output:
    path("reads_clean.qza"), emit: clean_reads_ch
    path("reads_clean.qzv")

    when:
    params.mode == "all" || params.mode == "remove_adapters"

    """
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path $manifest \
        --output-path reads_raw.qza \
        --input-format PairedEndFastqManifestPhred33V2

    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences reads_raw.qza \
        --p-cores $task.cpus \
        --p-front-f $params.primer_forward \
        --p-front-r $params.primer_reverse \
        --o-trimmed-sequences reads_clean.qza \
        --verbose &> log_trim.txt

    qiime demux summarize --i-data reads_clean.qza --o-visualization reads_clean.qzv
    """
}