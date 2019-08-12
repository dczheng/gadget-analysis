#!/usr/bin/env python3
import re

lines = ''.join(open( 'allvars.h' ).readlines())

i, j = re.search( 'GlobalParams.*GlobalParams;', lines, re.DOTALL ).span()
s = lines[i+1:j-1]
i, j = re.search( '\{.*\}', s, re.DOTALL ).span()
s = s[i+1:j-1]
s = ''.join( ''.join(s).split() )
s = s.split( ';' )[:-1]

ts = [ 'char', 'int', 'double' ]
tt = {'char':'STRING', 'int':'INT', 'double':'REAL'}
my_vars = {}
for i in range(len(s)):
    for j in range( len(ts) ):
        if ts[j] in s[i]:
            if ts[j] in my_vars.keys():
                my_vars[ ts[j] ] += s[i][len(ts[j]):].split( ',' )
            else:
                my_vars[ ts[j] ] = s[i][len(ts[j]):].split( ',' )

f = open( 'add_params.h', 'w' )
f.write( '#define ADD_PARAMS(){\\\n' )
for k in my_vars.keys():
    for v in my_vars[k]:
        t = re.search( '\[.*\]', v )
        if t:
            v = v[:t.span()[0]]
        f.write( "\tstrcpy( tag[nt], \"%s\" );\\\n"%v )
        f.write( "\taddr[nt] =&All.%s;\\\n"%v )
        f.write( "\tid[nt++] = %s;\\\n\\\n"%tt[k] )
f.write( '}\n' )
f.close()

