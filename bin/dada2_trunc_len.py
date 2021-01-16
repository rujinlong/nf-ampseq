#!/usr/bin/env python

# Author: Jinlong Ru <jinlong.ru@gmail.com>
# Usage: dada2_trunc_len.py --help

import click

@click.command()
@click.option("--fwd_98", '-f', type=int, help="98% the forward length")
@click.option("--rvs_98", '-r', type=int, help="98% the reverse length")
@click.option("--len_insert", '-i', type=int, help="Average insert length of primer")
@click.option("--len_overlap", '-o', type=int, default=20, help="Overlap length")
@click.option("--fwd_cut_pct", '-p', type=float, default=0.25, help="Percent of cut base in forward")
def main(fwd_98, rvs_98, len_insert, len_overlap, fwd_cut_pct):
    len_cut = fwd_98 + rvs_98 - len_insert - len_overlap
    cut_fwd_base = int(len_cut * fwd_cut_pct)
    cut_rvs_base = len_cut - cut_fwd_base
    trunc_len_f = fwd_98 - cut_fwd_base
    trunc_len_r = rvs_98 - cut_rvs_base
    print("--p-trunc-len-f: {}\n--p-trunc-len-r: {}".format(trunc_len_f, trunc_len_r))
    

if __name__ == '__main__':
    main()