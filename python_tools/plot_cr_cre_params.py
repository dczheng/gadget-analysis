#!/usr/bin/env python3

from my_work_env import *
import tools_and_constants as tc

f = h5py.File( sys.argv[1] )

labels = [ \
    r'$log_{10}{C}$', \
    r'$log_{10}{e}$', \
    r'$log_{10}{n}$', \
    r'$log_{10}(q_{\rm min})$', \
    r'$log_{10}(q_{\rm max})$'
]

dn = [ 'C', 'e', 'n', 'qmin', 'qmax' ]

N = 5

print( 'load CRE data ...' )
cre = [ f[ '/PartType0/CRE_C' ].value, \
        f[ '/PartType0/CRE_e' ].value, \
        f[ '/PartType0/CRE_n' ].value, \
        f[ '/PartType0/CRE_qmin' ].value, \
        f[ '/PartType0/CRE_qmax' ].value ]

print( 'load CR data ...' )
cr = [ f[ '/PartType0/CR_C0' ].value, \
    f[ '/PartType0/CR_E0' ].value, \
    f[ '/PartType0/CR_n0' ].value, \
    f[ '/PartType0/CR_q0' ].value ]

cre_alpha = f[ '/PartType0/CRE_Alpha' ].value


dens = f[ '/PartType0/Density' ].value

print( 'proc CRE params ...' )
print( tc.m_e / tc.gadget_mass_in_g )
print( tc.gadget_length_in_cm )

cre[0] = cre[0] * dens / ( tc.m_e / tc.gadget_mass_in_g * tc.gadget_length_in_cm ** 3 )
cre[1] = cre[1] * dens * tc.gadget_energy_in_erg / ( tc.gadget_length_in_cm ** 3 )
cre[2] = cre[2] * dens / ( tc.m_e / tc.gadget_mass_in_g * tc.gadget_length_in_cm ** 3 )

print( 'proc CR params ...' )
cr[0] = cr[0] * np.power( dens, (2.75-1)*0.333333 ) / ( tc.m_p / tc.gadget_mass_in_g * tc.gadget_length_in_cm ** 3 )
cr[1] = cr[1] * dens * tc.gadget_energy_in_erg  / ( tc.gadget_length_in_cm ** 3 )
cr[2] = cr[2] * dens / ( tc.m_p / tc.gadget_mass_in_g * tc.gadget_length_in_cm ** 3 )


for i in range( N ):
    t = cre[i][np.isnan( cre[i] )]
    tN = len(t)
    if tN != 0:
        print( "WARNING: cre %s contain `nan` [%i]"%(dn[i], tN) )
        cre[i] = cre[i][ np.isfinite( cre[i] ) ]

    if i == N-1:
        break

    t = cr[i][np.isnan( cr[i] )]
    tN = len(t)
    if tN != 0:
        print( "WARNING: cr %s contain `nan` [%i]"%(dn[i], tN) )
        cr[i] = cr[i][ np.isfinite( cr[i] ) ]

for i in range( N-1 ):
    cre[i] = np.log10( cre[i][cre[i]>0] )
    cr[i] = np.log10( cr[i][cr[i]>0] )
cre[N-1] = np.log10( cre[N-1][cre[N-1]>0] )


fig, axs = plt.subplots( 2, 5, figsize=(3*N, 3) )

for i in range( N ):
    print( "plot cre [%s]"%(dn[i]) )
    print( "min: %g, max: %g"%( cre[i].min(), cre[i].max() ) )
    axs[0, i].hist( cre[i], bins=20 )
    axs[0, i].set_yscale( 'log' )
    axs[0, i].set_xlabel( labels[i] )
    if i == 0:
        axs[0, i].set_ylabel( r'$N$' )
        #axs[0, i].set_xlim( -20, 0 )
axs[1, N-1].hist( cre_alpha[cre_alpha>0], bins=20 )
axs[1, N-1].set_yscale( 'log' )
axs[1, N-1].set_xlabel( r'$\alpha$' )

for i in range( N-1 ):
    print( "plot cr [%s]"%(dn[i]) )
    print( "min: %g, max: %g"%( cr[i].min(), cr[i].max() ) )
    axs[1, i].hist( cr[i], bins=20 )
    axs[1, i].set_yscale( 'log' )
    axs[1, i].set_xlabel( labels[i] )
    if i == 0:
        axs[1, i].set_ylabel( r'$N$' )
        #axs[1, i].set_xlim( -20, 0 )
fig.tight_layout()

if ( len(sys.argv) > 2 ):
    fn = output_dir + sys.argv[2] + '_cr_cre_params.png'
else:
    fn = output_dir + sys.argv[1][:-5] + '_cr_cre_params.png'

print( "savefig to %s"%fn )
plt.savefig( fn )
