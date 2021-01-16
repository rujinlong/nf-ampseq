#!/usr/bin/env Rscript

library(optparse)
library(qiime2R)
library(phyloseq)
library(tidyverse)

option_list = list(
    make_option(c("-f", "--features"), type="character", default="table.qza", 
                help="Abundance table", metavar="character"),
    make_option(c("-t", "--tree"), type="character", default="rooted_tree.qza", 
                help="Rooted tree", metavar="character"),
    make_option(c("-a", "--taxa"), type="character", default="taxonomy_rep_seqs.qza", 
                help="Taxonomy of rep seqs", metavar="character"),
    make_option(c("-m", "--meta"), type="character", default="metadata.tsv", 
                help="Metadata file", metavar="character"),
    make_option(c("-s", "--repseq"), type="character", default="dna-sequences.fasta", 
                help="Repseq FASTA", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="physeq.RData", 
                help="Output physeq object file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

physeq1 <- qza_to_phyloseq(
    features=opt$features,
    tree=opt$tree,
    taxonomy=opt$taxa,
    metadata=opt$meta
)

rep_seqs <- Biostrings::readDNAStringSet(opt$repseq, format = "fasta")
physeq2 <- phyloseq(rep_seqs)
physeq <- merge_phyloseq(physeq1, physeq2)
save(physeq, file = opt$out)
