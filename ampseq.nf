#!/usr/bin/env nextflow
// Usage: nextflow run ampseq.nf -profile hpc_slurm -resume --mode "remove_adapters"
//        nextflow run ampseq.nf -profile hpc_slurm --mode "all" -resume

nextflow.enable.dsl=2

process train_classifier {
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

process prepare_input {
    publishDir "$baseDir", mode: 'copy', overwrite: true
    publishDir "$params.report", pattern: "metadata.tsv"

    output:
    path("manifest.tsv"), emit: manifest_ch
    path("metadata.tsv"), emit: metadata_ch

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

process denoise {
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

process taxonomy {
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

process tree {
    label "big"
    publishDir "$params.outdir/p03_tree"
    publishDir "$params.report", pattern: "rooted_tree.qza"

    input:
    path(repseq)

    output:
    path("tree/*")
    path("rooted_tree.qza"), emit: tree_ch

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
    path(table)
    path(tree)
    path(taxseq)
    path(metadata)
    path(repseq)

    output:
    path("physeq.Rdata")

    """
    qiime2phyloseq.R -f $table -t $tree -a $taxseq -m $metadata -s $repseq -o physeq.Rdata
    """
}

workflow {
    prepare_input()
    remove_adapters(prepare_input.out.manifest_ch)
    denoise(remove_adapters.out.clean_reads_ch)

    if( !file(params.classifier).exists() ) {
        train_classifier()
        classifier_ch = train_classifier.out.classifier_ch
    }
    else {
        classifier_ch = params.classifier
    }

    taxonomy(classifier_ch, denoise.out.repseq_ch, denoise.out.table_ch, prepare_input.out.metadata_ch, denoise.out.table2taxa_abundance)
    tree(denoise.out.repseq_ch)
    phyloseq(denoise.out.table_ch, tree.out.tree_ch, taxonomy.out.taxa_repseq, prepare_input.out.metadata_ch, denoise.out.repseq)
}
