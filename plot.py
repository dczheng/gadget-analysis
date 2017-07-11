#!/usr/bin/env python2
import numpy as np
import matplotlib.pyplot as plt
import sys

data = np.loadtxt( sys.argv[1] )
plt.imshow( data )
plt.colorbar()
plt.savefig( sys.argv[2] )
