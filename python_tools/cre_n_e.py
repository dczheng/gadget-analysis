#!/usr/bin/env python3

from elec_cooling import *

c = float( sys.argv[1] )
a = float( sys.argv[2] )
q1 = float( sys.argv[3] )
q2 = float( sys.argv[4] )

print( tc.m_e / tc.gadget_mass_in_g )
print( 'n: %g, e: %g, n2: %g, e2: %g, T: %g' \
        %( \
            cre_n( c, a, q1 ), cre_e_gadget( c, a, q1 ),
            cre_n2( c, a, q1, q2 ), cre_e_gadget_2( c, a, q1, q2 ),
            cre_Tbar_gadget( a, q1 ) ) )
