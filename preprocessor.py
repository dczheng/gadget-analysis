#!/usr/bin/env python3
import re
    
debug = 0
def myprint( s ):
    if debug:
        print(s)

def gen_add_params():

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
    
    #print( my_vars )
    n = 0
    for k in my_vars.keys():
        n += len( my_vars[k] )
    #print( n )

    f = open( 'add_params.h', 'w' )
    f.write( '#define MAXTAGS %i\n'%n )
    f.write( '#define ADD_PARAM_SINGLE(s, t){\\\n' )
    f.write( "\tstrcpy( tag[nt], #s );\\\n" )
    f.write( "\taddr[nt] =&All.s;\\\n" )
    f.write( "\tid[nt++] = t;\\\n\\\n" )
    f.write( '}\n\n' )
    f.write( '#define ADD_PARAMS(){\\\n' )
    for k in my_vars.keys():
        for v in my_vars[k]:
            t = re.search( '\[.*\]', v )
            if t:
                v = v[:t.span()[0]]
            f.write( "\tADD_PARAM_SINGLE( %s, %s );\\\n"%( v, tt[k] ) )
    f.write( '}\n' )
    f.close()
    
def gen_allvars():

    lines = ''.join(open( 'allvars.h' ).readlines())
    
    i, j = re.search( 'extern_start.*extern_end', lines, re.DOTALL ).span()
    s = lines[i:j]
    s = s.split( 'extern' )
    f = open( 'allvars.c', 'w' )
    f.write( '#include "allvars.h"\n\n' )
    for l in s:
        if ';' in l:
            f.write( l.strip().split( '//' )[0] + '\n' )
    f.close()

def make_protos_single( f ):

    ms = [ '#include', '#define', '#if', '#un', '#endif' ]
    lines = open( f ).readlines()

    ll = []
    for l in lines:
        if l.split() == []:
            continue
        ll.append( l )
    lines = ll

    n = len(lines)
    i = 0
    ss = []
    while i < n:
        '''
        remove macros and `\\` comments
        '''
        l = lines[i].strip()

        flag = 0
        for m in ms:
            if m in l:
                flag = 1

        i += 1
        if flag:
            if l[-1] != '\\':
                continue

            for x in range(i, n):
                t = lines[x].strip()
                if t[-1] == '\\':
                    i += 1
                else:
                    break
            i += 1
            continue

        t = re.search( '//', l )
        if t:
            ss.append( l[:t.span()[0]] )
        else:
            ss.append( l )
    s = ''.join( ss ) 
    myprint( s + '\n' + '-' * 80 )

    def remove_some_char( s, c0, c1, e=None, l=1 ):
        '''remove block'''
        ss = []
        n = len(s)
        i = 0
        while i < n:
            c = s[i:i+l]
            if c == c0 or c == c1:
                i += l 
            else:
                c = s[i]
                i += 1

            if c != c1:
                ss.append( c )
                continue

            cc = ss.pop()
            while cc != c0 and len(ss) != 0:
                cc = ss.pop()
            if e:
               ss.append( e )
        return ''.join(ss)

    s = remove_some_char( s, '{', '}', e=';' )
    myprint( s + '\n' + '-' * 80 )

    s = remove_some_char( s, r'/*', r'*/', l=2  )
    myprint( s + '\n' + '-' * 80 )

    s = s.split( ';' )
    ss = []
    for l in s:
        if '(' in l:
            myprint( l )
            ss.append( l.strip() + ';' )
    myprint( ss ) 
    myprint( '\n' + '-' * 80 )

    return ss
    

def make_protos():
    
    from os import listdir

    fs = listdir( '.' )

    #s = make_protos_single( "./check_flag.c" )
    #myprint( s )
    #exit()

    ff = open( 'protos.h', 'w' )
    for f in fs:
        if '.c' not in f:
            continue

        #print( f )
        info = '//%%------>>>>>>file : %s\n//%%\n'%f
        ff.write( info )
        s = make_protos_single( f )

        for l in s:
            ff.write( "%s\n"%l )
        ff.write( '\n\n' )
    ff.close()


make_protos()
gen_add_params()
gen_allvars()
