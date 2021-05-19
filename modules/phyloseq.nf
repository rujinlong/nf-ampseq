process PHYLOSEQ {
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