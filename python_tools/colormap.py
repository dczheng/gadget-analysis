#!/usr/bin/env python3

from matplotlib import use
use( 'agg' )
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

plt.rc('text', usetex=False)

r = 10
N = 1000 
a=np.outer( np.ones(N//r), np.arange(0,1,1/N))

maps=[m for m in cm.datad if not m.endswith("_r")]
maps.sort()
l=len(maps)
w = 5
fig, axs = plt.subplots( l, 1, figsize=(w, w/r*l) ) 
for i, m in enumerate(maps):
    axs[i].axis( 'off' )
    axs[i].imshow(a,cmap=plt.get_cmap(m))
    axs[i].text( N, 0.4*N/r, m )
plt.savefig("colormaps.png",dpi=100 )
