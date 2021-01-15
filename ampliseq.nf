#!/usr/bin/env nextflow
// Usage: nextflow run ampliseq.nf -profile hpc_slurm --mode "remove_adapters"
//        nextflow run ampliseq.nf -profile hpc_slurm --mode "all" -resume

if (params.classifier == false) {
    process train_classifier {
        publishDir "$params.outdir"

        output:
        file("classifier.qza") into classifier_ch

        when:
        params.mode == "all"

        """
        qiime tools import \
            --type 'FeatureData[Taxonomy]' \
            --input-format HeaderlessTSVTaxonomyFormat \
            --input-path $params.SILVA_taxonomy \
            --output-path p01_SILVA_99_otus_16S_taxonomy.qza

        qiime tools import \
            --type 'FeatureData[Sequence]' \
            --input-path $params.SILVA_sequence \
            --output-path p01_SILVA_99_otus_16S_seqs.qza

        qiime feature-classifier extract-reads \
            --i-sequences p01_SILVA_99_otus_16S_seqs.qza \
            --p-f-primer $params.primer_forward \
            --p-r-primer $params.primer_reverse \
            --p-min-length $params.min_reference_length \
            --p-max-length $params.max_reference_length \
            --o-reads p02_SILVA_99_otus_16S_Vx_seqs.qza \
            --p-n-jobs $task.cpus \
            --verbose \
            &> z01_SILVA_99_otus_16S_V4_seqs.log

        qiime feature-classifier fit-classifier-naive-bayes \
            --i-reference-taxonomy p01_SILVA_99_otus_16S_taxonomy.qza \
            --i-reference-reads p02_SILVA_99_otus_16S_Vx_seqs.qza \
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
    publishDir "$params.outdir"
    publishDir "$params.report", pattern: "*.tsv"

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
    publishDir "$params.outdir/p01_no_adapters"

    input:
    file(manifest) from manifest_ch

    output:
    file("p02_primer_trimmed.qza") into clean_reads_ch
    file("p02_primer_trimmed.qzv")

    when:
    params.mode == "all" || params.mode == "remove_adapters"

    """
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path $manifest \
        --output-path p01_raw_reads.qza \
        --input-format PairedEndFastqManifestPhred33V2

    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences p01_raw_reads.qza \
        --p-cores $task.cpus \
        --p-front-f $params.primer_forward \
        --p-front-r $params.primer_reverse \
        --o-trimmed-sequences p02_primer_trimmed.qza \
        --verbose &> z01_trim.log

    qiime demux summarize --i-data p02_primer_trimmed.qza --o-visualization p02_primer_trimmed.qzv
    """
}


process denoise {
    label "big"
    publishDir "$params.outdir/p02_denoise"

    input:
    file(clean_reads) from clean_reads_ch

    output:
    file("table/*")
    file("representative_sequences/*")
    file("representative_sequences/dna-sequences.fasta") into repseq2phyloseq_ch
    file("dada2_output/representative_sequences.qza") into repseq_ch1
    file("dada2_output/representative_sequences.qza") into repseq_ch2
    file("dada2_output/table.qza") into table_ch
    file("dada2_output/table.qza") into table2phyloseq_ch
    file("denoising_stats/*")

    when:
    params.mode == "all"

    """
    qiime dada2 denoise-paired --p-n-threads $task.cpus \
        --i-demultiplexed-seqs $clean_reads \
        --p-trunc-len-f $params.trunc_len_forward \
        --p-trunc-len-r $params.trunc_len_reverse \
        --p-max-ee-f 3 \
        --p-max-ee-r 3 \
        --p-min-fold-parent-over-abundance 3 \
        --output-dir dada2_output \
        --verbose &> z02_denoise.log

    # Denoising stats
    qiime metadata tabulate \
        --m-input-file dada2_output/denoising_stats.qza \
        --o-visualization denoising_stats.qzv

    # Representative sequences
    qiime feature-table tabulate-seqs \
        --i-data dada2_output/representative_sequences.qza \
        --o-visualization representative_sequences.qzv

    # Feature table
    qiime feature-table summarize \
        --i-table dada2_output/table.qza \
        --o-visualization table.qzv
    
    qiime tools export --input-path denoising_stats.qzv --output-path denoising_stats
    qiime tools export --input-path dada2_output/table.qza --output-path table
    qiime tools export --input-path dada2_output/representative_sequences.qza --output-path representative_sequences
    """
}


process taxonomy {
    label "medium"
    publishDir "$params.outdir/p03_taxonomy"

    input:
    file(classifier) from classifier_ch
    file(repseq) from repseq_ch1
    file(table) from table_ch
    file(metadata) from metadata2taxonomy_ch

    output:
    file("taxonomy_rep_seqs*")
    file("taxonomy_barplots*")
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
    """
}


process tree {
    label "small"
    publishDir "$params.outdir/p04_tree"
    publishDir "$params.report", pattern: "*.stats"

    input:
    file(repseq) from repseq_ch2

    output:
    file("tree/*")
    file("tree/rooted_tree.qza") into tree2phyloseq_ch

    when:
    params.mode == "all"

    """
    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences $repseq \
        --output-dir tree \
        --p-n-threads $task.cpus \
        --verbose \
        &> z01_phylogenetic_tree_generation.log
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

