process TREE {
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