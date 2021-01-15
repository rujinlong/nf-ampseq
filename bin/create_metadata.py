#!/usr/bin/env python

import pandas as pd
import click

def csv2tsv(df):
    first_col = df.columns[0]
    df.rename(columns={first_col:'sample-id'}, inplace=True)
    return df

@click.command()
@click.option("--fin", '-i', help="metadata.csv")
@click.option("--fout", '-o', help="metadata.tsv")
def main(fin, fout):
    df = pd.read_csv(fin)
    df = csv2tsv(df)
    df.to_csv(fout, sep='\t', index=False)


if __name__ == '__main__':
    main()