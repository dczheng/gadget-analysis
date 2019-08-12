#!/usr/bin/env python3
import re

lines = ''.join(open( 'allvars.h' ).readlines())

i, j = re.search( 'extern_start.*extern_end', lines, re.DOTALL ).span()
s = lines[i:j]
s = s.split( 'extern' )
f = open( 'allvars.c', 'w' )
f.write( '#include "allvars.h"\n\n' )
for l in s:
    if ';' in l:
        f.write( l.strip() + '\n' )
f.close()

