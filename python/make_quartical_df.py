# coding: utf-8
import numpy as np
import pandas as pd

from daskms.experimental.zarr import xds_from_zarr

def get_args():
    import argparse
    agp   = argparse.ArgumentParser('make_quartical_df', description="make_dataframe_from_quartical. We assume timechunk=0 when solving.")
    add   = agp.add_argument
    add ('path', help='Path to quartical solution')
    add ('-g', '--gain', help='Gain term name given in yaml (G or KCROSS)', dest='term', default='G')
    add ('-o', '--ofile', help='Output tag', dest='otag', required=True)
    return agp.parse_args()

if __name__ == "__main__":
    args     = get_args()
    gains    = xds_from_zarr(f"{args.path}::{args.term}")

    for i,gg in enumerate(gains):
        ## 
        ## time, freq, antenna, direction, correlation
        gf       = gg.gains[0,:,:,0,:].to_dataframe (dim_order=['antenna','correlation','gain_freq'])
        gflags   = gg.gain_flags[0,:,:,0].to_dataframe(dim_order=['antenna','gain_freq'])
        gf       = gf.join ( gflags['gain_flags'] )
        gf.to_csv(f"{args.otag}_{args.term}_{i:d}.csv")


