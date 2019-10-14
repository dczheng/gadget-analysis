#!/usr/bin/env python3

import sys
import pandas as pd

if len( sys.argv ) < 5:
    exit()

fn = sys.argv[1]
ds = sys.argv[2:-2]
N  = int( sys.argv[-2] )
w  = int( sys.argv[-1] )

print( fn )
print( '*' * len(fn) )

fmt1 = "%%%is "%w
fmt2 = "%%%ig "%w
print( "%6s "%('index'), end='' )
for i in ds:
    print( fmt1%i, end='' )
print()

d = pd.read_csv( fn )[ds].to_numpy()
m, n = d.shape
print( '-' * ( 6 + (w+1)*len(ds) ) )
for i in range(m):
    if N>0 and i>=N:
        break
    print( "%6i "%(i), end='' )
    for j in range(n):
        print( fmt2%d[i,j], end='' )
    print()

print( '-' * ( 6 + (w+1)*len(ds) ) )
print( "%6s "%('index'), end='' )
for i in ds:
    print( fmt1%i, end='' )
print()
