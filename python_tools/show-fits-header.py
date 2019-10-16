#!/usr/bin/env python3

from astropy.io import fits
import sys

h = fits.open( sys.argv[1] )[0].header
removed_keys = [ 'HISTORY', 'COMMENT' ]

ks = []
vs = []

for k in h.keys():
    if k in removed_keys:
        continue
    ks.append( k )
    vs.append( str(h[k]) )
lk = max([ len(k) for k in ks ])
lv = max([ len(v) for v in vs ])

#print( lk, lv )

fmt='%%%is   %%-%is'%( lk+1, lv+1 )
print( fmt%("key", 'value') )
for i in range(len(ks)):
    print( fmt%(ks[i], vs[i]) )
 

