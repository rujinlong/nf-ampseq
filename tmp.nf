#!/usr/bin/env nextflow
// Usage: nextflow main.nf --reads_dir fq_seq_dir --prefix mysample --outdir outputdir --bin_path <bin_path> -resume
// params.bin_path path of the software, ex: ~/conda3/envs/binfo/bin
// params.mode ['QA', 'QC', 'complete']

// TOOLS
DIR_SCRIPT="${baseDir}/script"
DIR_BIN = params.bin_path
params.dir_diversity = "diversity_41000"
params.meta_col = "bee_type"

//Creates working dir
workingpath = params.outdir + "/" + params.prefix
workingdir = file(workingpath)
if( !workingdir.exists() ) {
    if( !workingdir.mkdirs() ) {
        exit 1, "Cannot create working directory: $workingpath"
    } 
}

// import all fastq reads into a qza
reads_dir = Channel.value(file(params.reads_dir, type: 'dir'))
process import_fq {
    input:
    file(dir_reads) from reads_dir

    output:
    file("reads.qza") into to_dada2_qa 
    file("reads.qza") into to_dada2_denoise_qza

    """
    ${DIR_BIN}/qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ${dir_reads} --output-path reads.qza
    """
}


process dada2_qa {
    input:
    file(reads) from to_dada2_qa

    output:
    file("qual_viz.qzv") into to_dada2_denoise_qzv

    """
    ${DIR_BIN}/qiime demux summarize --p-n 10000 --i-data ${reads} --o-visualization qual_viz
    """
}

// Time-consuming
process dada2_denoise {
    input:
    file("reads") from to_dada2_denoise_qza
    file("qual_viz") from to_dada2_denoise_qzv

    output:
    file("table.qza") into to_viz_sps
    file("table.qza") into to_filter
    file("table.qza") into to_rare
    file("representative_sequences.qza") into to_taxa
    file("representative_sequences.qza") into to_phylo


    """
    ${DIR_BIN}/qiime dada2 denoise-paired --i-demultiplexed-seqs ${reads} --o-table table --o-representative-sequences representative_sequences --p-trunc-len-f 151 --p-trunc-len-r 140 --p-trim-left-f 19 --p-trim-left-r 20 --p-n-threads 0 --verbose
    """
}


// This visualization shows us the sequences/sample spread
process viz_sps_spread {
    input:
    file(table) from to_viz_sps

    output:
    file("table_summary.qzv") into to_other1

    """
    ${DIR_BIN}/qiime feature-table summarize --i-table ${table} --o-visualization table_summary
    """
}


// Filter out sequences with few samples
process filter_samples {
    input:
    file(table) from to_filter

    output:
    file("filtered_table.qza") into to_taxa_table
    file("filtered_table.qza") into to_diversity_table
    file("filtered_table_summary.qzv") into to_other2

    """
    ${DIR_BIN}/qiime feature-table filter-samples --i-table table.qza --p-min-frequency 5000 --o-filtered-table filtered_table
    ${DIR_BIN}/qiime feature-table summarize --i-table filtered_table.qza --o-visualization filtered_table_summary
    """
}

// Download green_gene dataset from https://docs.qiime2.org/2019.1/data-resources
ch_green_gene = Channel.value(file(params.green_gene))
Channel.value(file(params.metadata)).into{ch_taxa_metadata; ch_diversity_metadata; ch_bg_metadata}
process taxa_classify {
    input:
    file("green_gene") from ch_green_gene
    file("filtered_table") from to_taxa_table
    file("rep_seqs") from to_taxa
    file("metadata") from ch_taxa_metadata

    output:
    file("taxa-bar-plots.qzv") into to_other3

    """
    ${DIR_BIN}/qiime feature-classifier classify-sklearn --i-classifier ${green_gene} --i-reads ${rep_seqs} --o-classification taxonomy
    ${DIR_BIN}/qiime taxa barplot --i-table ${filtered_table} --i-taxonomy taxonomy.qza --m-metadata-file ${metadata} --o-visualization taxa-bar-plots
    
    """
}


process phylo_tree {
    input:
    file("rep_seqs") from to_phylo

    output:
    file("rooted_tree.qza") into to_diversity

    """
    ${DIR_BIN}/qiime alignment mafft --i-sequences ${rep_seqs} --o-alignment aligned_representative_sequences
    ${DIR_BIN}/qiime alignment mask --i-alignment aligned_representative_sequences.qza --o-masked-alignment masked_aligned_representative_sequences
    ${DIR_BIN}/qiime phylogeny fasttree --i-alignment masked_aligned_representative_sequences.qza --o-tree unrooted_tree
    ${DIR_BIN}/qiime phylogeny midpoint-root --i-tree unrooted_tree.qza --o-rooted-tree rooted_tree
    """
}


process diversity_measures {
    input:
    file(root_tree) from to_diversity
    file(filtered_table) from to_diversity_table
    file(metadata) from ch_diversity_metadata

    output:
    file(params.dir_diversity) into to_diff_bg
    file(params.dir_diversity) into to_diff_rare

    """
    ${DIR_BIN}/qiime diversity core-metrics-phylogenetic --i-phylogeny ${root_tree} --i-table ${filtered_table} --p-sampling-depth 41000 --output-dir params.dir_diversity --m-metadata-file ${metadata}
    """
}


// Test for between-group differences
// Alpha rarefaction curves show taxon accumulation as a function of sequence depth
process diff_between_group {
    input:
    file(dir_diversity) from to_diff_bg
    file(metadata) from ch_bg_metadata
    file(table) from to_rare

    output:
    file("rooted_tree2.qza") into qiime_others

    """
    ${DIR_BIN}/qiime diversity alpha-group-significance --i-alpha-diversity ${dir_diversity}/faith_pd_vector.qza --m-metadata-file ${metadata} --o-visualization ${dir_diversity}/alpha_PD_significance
    ${DIR_BIN}/qiime diversity alpha-group-significance --i-alpha-diversity ${dir_diversity}/shannon_vector.qza --m-metadata-file ${metadata} --o-visualization ${dir_diversity}/alpha_shannon_significance
    ${DIR_BIN}/qiime diversity beta-group-significance --i-distance-matrix ${dir_diversity}/bray_curtis_distance_matrix.qza --m-metadata-file ${metadata} --m-metadata-column params.meta_col --o-visualization ${dir_diversity}/beta_bray_beetype_significance
    ${DIR_BIN}/qiime diversity alpha-rarefaction --i-table ${table} --p-max-depth 41000 --o-visualization ${dir_diversity}/alpha_rarefaction.qzv --m-metadata-file ${metadata} --i-phylogeny rooted_tree2.qza
    """
}


workflow.onComplete { 
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
