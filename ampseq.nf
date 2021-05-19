#!/usr/bin/env nextflow
// Usage: nextflow run ampseq.nf -profile hpc_slurm -resume --mode "remove_adapters"
//        nextflow run ampseq.nf -profile hpc_slurm --mode "all" -resume

nextflow.enable.dsl=2

include { TRAIN_CLASSIFIER } from './modules/train_classifier'
include { PREPARE_INPUT } from './modules/prepare_input'
include { REMOVE_ADAPTERS } from './modules/remove_adapters'
include { DENOISE } from './modules/denoise'
include { TAXONOMY } from './modules/taxonomy'
include { TREE } from './modules/phylogenetic_tree'
include { PHYLOSEQ } from './modules/phyloseq'

workflow {
    PREPARE_INPUT()
    REMOVE_ADAPTERS(PREPARE_INPUT.out.manifest_ch)
    DENOISE(REMOVE_ADAPTERS.out.clean_reads_ch)

    if( !file(params.classifier).exists() ) {
        TRAIN_CLASSIFIER()
        classifier_ch = TRAIN_CLASSIFIER.out.classifier_ch
    }
    else {
        classifier_ch = params.classifier
    }

    TAXONOMY(classifier_ch, DENOISE.out.repseq_ch, DENOISE.out.table_ch, PREPARE_INPUT.out.metadata_ch, DENOISE.out.table2taxa_abundance)
    TREE(DENOISE.out.repseq_ch)
    PHYLOSEQ(DENOISE.out.table_ch, TREE.out.tree_ch, TAXONOMY.out.taxa_repseq, PREPARE_INPUT.out.metadata_ch, DENOISE.out.repseq)
}
