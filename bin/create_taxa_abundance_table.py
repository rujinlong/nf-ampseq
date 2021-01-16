#!/usr/bin/env python

import click
import pandas as pd


@click.command()
@click.option("--feature_abundance", '-f', help="feature-table.tsv")
@click.option("--taxa_table", '-t', help="taxonomy_rep_seqs/metadata.tsv")
@click.option("--taxa_abundance", '-o', help='taxa_abundance.tsv')
def main(feature_abundance, taxa_table, taxa_abundance):
    df_feature = pd.read_csv(feature_abundance, sep='\t', skiprows=1)
    df_feature.rename(columns={"#OTU ID": "OTU_ID"}, inplace=True)
    
    df_taxa = pd.read_csv(taxa_table, sep='\t', skiprows=2, names=['OTU_ID', 'Taxon', 'Confidence'])
    
    tbl = df_taxa.merge(df_feature, on='OTU_ID', how='inner')
    tbl.to_csv(taxa_abundance, index=False, sep='\t')


if __name__ == '__main__':
    main()