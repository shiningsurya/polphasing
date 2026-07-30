# coding: utf-8
import numpy as np
import pandas as pd

def get_args():
    import argparse
    agp = argparse.ArgumentParser("parse_bldata")
    add = agp.add_argument
    add ('bldata', help='Output of :write_baselines:')
    return agp.parse_args()


if __name__ == "__main__":
    args = get_args ()
    bl   = pd.read_csv (args.bldata, sep='\\s+')
    bl['complex'] = bl['complex'].apply(complex)
    gl   = bl.groupby(['ant1','ant2','correlation']).agg(list)
    #####
    gl.to_pickle(f"{args.bldata}_df.pkl")

