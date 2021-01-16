#!/usr/bin/env nextflow
// Usage: nextflow run ampliseq.nf -profile hpc_slurm --mode "remove_adapters"
//        nextflow run ampliseq.nf -profile hpc_slurm --mode "all" -resume

if (params.classifier == false) {
    process train_classifier {
        publishDir "$baseDir", mode: 'copy', overwrite: false

        output:
        file("classifier.qza") into classifier_ch

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
} else {
    Channel
        .fromPath(params.classifier)
        .set { classifier_ch }
}


process prepare_input {
    publishDir "$baseDir", mode: 'copy', overwrite: true
    publishDir "$params.report", pattern: "metadata.tsv"

    output:
    file("manifest.tsv") into manifest_ch
    file("metadata.tsv") into metadata2taxonomy_ch
    file("metadata.tsv") into metadata2phyloseq_ch

    """
    echo "sampleid,forward,reverse" > manifest.csv

    for sid in \$(ls $params.raw_seq_dir | sed 's/_R.*//' | sort -u);do
        echo "\$sid,\${sid}_R1${params.raw_seq_ext},\${sid}_R2${params.raw_seq_ext}" >> manifest.csv
    done

    create_manifest.py -i manifest.csv -d $params.raw_seq_dir -o manifest.tsv
    create_metadata.py -i $params.metadata -o metadata.tsv
    """
}


process remove_adapters {
    label "small"
    publishDir "$params.outdir/p01_clean_reads"

    input:
    file(manifest) from manifest_ch

    output:
    file("reads_clean.qza") into clean_reads_ch
    file("reads_clean.qzv")

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


process denoise {
    label "big"
    publishDir "$params.outdir/p02_denoise"
    publishDir "$params.report", pattern: "denoising_stats*"
    publishDir "$params.report", pattern: "representative_sequences*"
    publishDir "$params.report", pattern: "table*"
    publishDir "$params.report", pattern: "feature_table.tsv"

    input:
    file(clean_reads) from clean_reads_ch

    output:
    file("denoising_stats")
    file("table")
    file("representative_sequences")
    file("repseqs.fasta") into repseq2phyloseq_ch
    file("representative_sequences.qza") into repseq_ch1
    file("representative_sequences.qza") into repseq_ch2
    file("table.qza") into table_ch
    file("table.qza") into table2phyloseq_ch
    file("feature_table.tsv") into table2taxa_abundance
    

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


process taxonomy {
    publishDir "$params.outdir/p03_taxonomy"
    publishDir "$params.report", pattern: "taxonomy_barplots"
    publishDir "$params.report", pattern: "taxonomy_rep_seqs"
    publishDir "$params.report", pattern: "taxa_abundance.tsv"

    input:
    file(classifier) from classifier_ch
    file(repseq) from repseq_ch1
    file(table) from table_ch
    file(metadata) from metadata2taxonomy_ch
    file(feature_table) from table2taxa_abundance

    output:
    file("taxa_abundance.tsv")
    file("taxonomy_rep_seqs")
    file("taxonomy_barplots")
    file("taxonomy_rep_seqs.qza") into taxseq2phyloseq_ch

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


process tree {
    publishDir "$params.outdir/p03_tree"
    publishDir "$params.report", pattern: "rooted_tree.qza"

    input:
    file(repseq) from repseq_ch2

    output:
    file("tree/*")
    file("rooted_tree.qza") into tree2phyloseq_ch

    when:
    params.mode == "all"

    """
    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences $repseq \
        --output-dir tree \
        --p-n-threads $task.cpus \
        --verbose \
        &> log_tree.log
    ln -s tree/rooted_tree.qza .
    """
}

process phyloseq {
    publishDir "$params.outdir"
    publishDir "$params.report"

    input:
    file(table) from table2phyloseq_ch
    file(tree) from tree2phyloseq_ch
    file(taxseq) from taxseq2phyloseq_ch
    file(metadata) from metadata2phyloseq_ch
    file(repseq) from repseq2phyloseq_ch

    output:
    file("physeq.Rdata")

    """
    qiime2phyloseq.R -f $table -t $tree -a $taxseq -m $metadata -s $repseq -o physeq.Rdata
    """
}


workflow.onComplete { 
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
