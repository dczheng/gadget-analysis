#!/usr/bin/env python3

import sys
import pandas as pd

fn_csv =  sys.argv[1]
fn_tex =  sys.argv[2]

cols  = [ 'idx_cre', 'mtot_cre', 'mtot_r', 'mgas_cre', 'mgas_r', 'mstar_cre', \
        'mstar_r', 'mhot_cre', 'mhot_r', 'mwarmhot_cre', 'mwarmhot_r']
names = [\
    'ID',\
    '$M_{\\rm tot}$',\
    '$M_{\\rm tot}$-diff',\
    '$M_{\\rm gas}$',\
    '$M_{\\rm gas}$-diff',\
    '$M_{\\rm star}$',\
    '$M_{\\rm star}$-diff',\
    '$R_{\\rm hot}$',\
    '$R_{\\rm hot}$-diff',\
    '$R_{\\rm warm-hot}$',\
    '$R_{\\rm warm-hot}$-diff',\
            ]
caption = "\n\
    Impact of CRE process on the component mass of\
    31 massive cluster.\n\
    $M_{\\rm tot}$ is the total mass of cluster,\n\
    $M_{\\rm gas}$ is the total gas mass of cluster,\n\
    $M_{\\rm star}$ is the total star mass of cluster,\n\
    $R_{\\rm warm-hot} = M_{\\rm warm-hot}/M_{\\rm gas}$,\n\
        where $M_{\\rm warm-hot}$ is the warm-hot gas mass of cluster,\n\
    $R_{\\rm hot} = M_{\\rm hot}/M_{\\rm gas}$,\n\
        where $M_{\\rm hot}$ is the hot gas mass of cluster.\n \
    Note: the unit of mass is $10^{13} \\, h^{-1}M_\\odot$ and \n \
    all mass is obtained from simulation SIM-CRE.\n\
    See Fig.~\\ref{fig:g_profile} and \\ref{fig:g_cre_profile}\
    of appendix \\ref{sec:visual} for\
    the gas density profile of those clusters in our two simulations.\
    "
d = pd.read_csv( fn_csv )
dd = d[ cols ]
#print( dd )
m, n = dd.shape

f_tex = open( fn_tex, 'w' )
f_tex.write( "\\begin{table*}\n" )
f_tex.write( "\\begin{adjustbox}{width=0.9\\textwidth}\n" )
t = 'c'
for i in range(1, len(cols)):
    if i%2:
        t += '|c'
    else:
        t += '|r'
l = '\\begin{tabular}{%s}'%t
f_tex.write( l+'\n' )

f_tex.write( '\\hline\n' )

l = ' & '.join(names) + " \\\\"
f_tex.write( l+'\n' )

f_tex.write( '\\hline\n' )

for i in range(m):
    a = dd.loc[i].to_numpy()
    l = [ '%i'%(a[0]+1) ]
    for j in range(1,len(a)):
        if j%2:
            l.append( '%5.2f'%a[j] )
        else:
            l.append( '%5.2f\\%%'%a[j] )
    l = ' & '.join(l) + " \\\\"
    f_tex.write( l+'\n' )

f_tex.write( '\\hline\n' )
a = dd.mean( axis=0 ).to_numpy()

l = [ 'min' ]
for j in range(1, len(names)):
    if j%2:
        t = ' --- '
    else:
        t = dd[ cols[j] ]
        t =  '%5.2f\\%%'%t.min()
    t = t.replace( ' ', '\\ ' )
    l.append(t)
l = ' & '.join(l) + " \\\\"
f_tex.write( l+'\n' )

l = [ 'max' ]
for j in range(1, len(names)):
    if j%2:
        t = ' --- '
    else:
        t = dd[ cols[j] ]
        t =  '%5.2f\\%%'%t.max()
    t = t.replace( ' ', '\\ ' )
    l.append(t)
l = ' & '.join(l) + " \\\\"
f_tex.write( l+'\n' )

l = [ 'mean' ]
for j in range(1, len(names)):
    if j%2:
        t = ' --- '
    else:
        t =  '%5.2f\\%%'%a[j] 
    t = t.replace( ' ', '\\ ' )
    l.append(t)
l = ' & '.join(l) + " \\\\"
f_tex.write( l+'\n' )


l = [ '+' ]
for j in range(1, len(names)):
    if j%2:
        t = ' --- '
    else:
        t = dd[ cols[j] ]
        t =  ' %2i '%(len(t[t>0])) 
    t = t.replace( ' ', '\\ ' )
    l.append(t)
l = ' & '.join(l) + " \\\\"
f_tex.write( l+'\n' )

l = [ '-' ]
for j in range(1, len(names)):
    if j%2:
        t = ' --- '
    else:
        t = dd[ cols[j] ]
        t =  ' %2i '%(len(t[t<0])) 
    t = t.replace( ' ', '\\ ' )
    l.append(t)
l = ' & '.join(l) + " \\\\"

print( l )
f_tex.write( l+'\n' )

f_tex.write( '\\hline\n' )

f_tex.write( '\\end{tabular}\n' )
f_tex.write( "\\end{adjustbox}\n" )
f_tex.write( '\\caption{%s}\n'%caption )
f_tex.write( '\\label{tab:g}\n' )
f_tex.write( "\\end{table*}\n" )
f_tex.close()
