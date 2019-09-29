#!/usr/bin/env python3

from my_work_env import *

f = h5py.File( sys.argv[1], 'r' )
z = f['Header'].attrs[ "Redshift" ]
h = f['/PartType0/SmoothingLength'][()]
T = f['/PartType0/Temperature'][()]
mm = f['Header'].attrs[ "MassTable" ][0]
if mm==0:
    m = f[ 'PartType0/Masses' ][()]
else:
    m = mm
d = f[ 'PartType0/Density' ][()]
divb = f[ 'PartType0/DivergenceOfMagneticField' ][()]
b = f[ 'PartType0/MagneticField' ][()]
b = np.sqrt(b[:,0]**2 + b[:,1]**2 + b[:,2]**2)
f.close()

def put_sep():
    print( '-' * 30 + '%.2f'%z + '-'*30 )

put_sep()
print( sys.argv[1] )

idx = (b>0) * (divb>0) 
h = h[idx]
if len(m>1):
    m = m[idx]
b = b[idx]
d = d[idx]
T = T[idx]
divb = np.fabs(divb[idx])
bmu = b * 1e6   # to muG

v = m / d
divberr = divb * h / b 

def f( v, v1, v2=None ):
    if v2:
        idx = (v>v1) * (v<v2)
    else:
        idx = (v>v1)
    vv = v[idx]
    if v2:
        print( "[%8.2f~%-8.2f]: %7.3f%% %i "%( v1, v2, len(vv)/len(v) * 100 , len(vv) ) )
    else:
        print( "[%8.2f~%-8s]: %7.3f%% %i"%( v1, '?', len(vv)/len(v) * 100, len(vv) ) )

print( "B info:" )
for i in range(-2, 4):
    f( bmu, 10**(i), 10**(i+1) )
f( bmu, 10**4 )
print( "Bmax: %g"%bmu.max() )

print( "divBerr info:" )
f( divberr, 0.1, 1 )
f( divberr, 1, 10 )
f( divberr, 10 )

print( "min: %g, max: %g, mean: %g"%(\
                divberr.min(), divberr.max(), (divberr*v).sum() / v.sum()) )

#idx = divberr>100
#dd = d[idx] / mycc.rho_bar_in_gadget(z)
#TT = T[idx]
#divberr2 = divberr[idx]
#for i in range(len( dd )):
#    print( divberr2[i], dd[i], TT[i] )

put_sep()
