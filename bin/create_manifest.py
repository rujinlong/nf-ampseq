#!/usr/bin/env python

import pandas as pd
import click
import os

@click.command()
@click.option("--fin", '-i', help="manifest.csv")
@click.option("--data_dir", '-d', help="data folder")
@click.option("--fout", '-o', help="manifest.tsv")
def main(fin, data_dir, fout):
    if data_dir[-1] == '/':
        data_dir = data_dir[:-1]
    data_path = os.path.join(os.getcwd(), data_dir)
    df = pd.read_csv(fin)
    df.columns = ["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"]
    df["forward-absolute-filepath"] = df.apply(lambda x:os.path.join(data_path, x["forward-absolute-filepath"]), axis=1)
    df["reverse-absolute-filepath"] = df.apply(lambda x:os.path.join(data_path, x["reverse-absolute-filepath"]), axis=1)
    df.to_csv(fout, sep='\t', index=False)

if __name__ == '__main__':
    main()