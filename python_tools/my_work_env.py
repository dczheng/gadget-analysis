#!/usr/bin/env python3

import matplotlib
matplotlib.use( 'agg' )

import numpy as np
import matplotlib.pyplot as plt
import h5py
import tools_and_constants as tc
import os
import sys
from matplotlib import cm
import matplotlib.colors as mplc

plt.rc( 'text', usetex=True )
plt.rc( 'font', family='serif' )

output_dir = "/mnt/ddnfs/data_users/dczheng/simulation/plots/"

