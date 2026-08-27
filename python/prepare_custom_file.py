# coding: utf-8
import numpy as np
import pandas as pd


df = pd.read_pickle("3C138_unphased_df.pkl")


"""
drop E06 it has nan
"""
mm = (df['ant1'] == 'E06') | (df['ant2'] == 'E06')

df = pd.DataFrame(df[~mm])

df = df.set_index(['ant1','ant2']).sort_index()

OFILE = "3C138_unphased.custom"

gg = np.array([])
ix = []
for idf in df.itertuples():
    ix.append ( f"{idf.Index[0]} {idf.Index[1]} rr\n" )
    gg = np.append(gg,idf.rr)
    ix.append ( f"{idf.Index[0]} {idf.Index[1]} rl\n" )
    gg = np.append(gg,idf.rl)
    ix.append ( f"{idf.Index[0]} {idf.Index[1]} lr\n" )
    gg = np.append(gg,idf.lr)
    ix.append ( f"{idf.Index[0]} {idf.Index[1]} ll\n" )
    gg = np.append(gg,idf.ll)

uu = np.complex64(gg)
uu.tofile(f"{OFILE}.raw")

with open (f"{OFILE}.htxt", 'w') as _f:
    _f.write("".join(ix))
