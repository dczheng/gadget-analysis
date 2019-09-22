#!/usr/bin/env python3
import re
import sys
import os
    
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

def deps1( fdh, fdc, As, B ):
    fdh.write( "#if defined(%s) "%As[0] )
    fdc.write( "#if defined(%s) "%As[0] )
    for a in As[1:]:
        fdh.write( "|| defined(%s) "%a )
        fdc.write( "|| defined(%s) "%a )
    fdh.write( "\n#ifndef %s\n"%B )
    fdh.write( "#define %s\n"%B )
    fdh.write( "#endif\n" )
    fdh.write( "#endif\n\n" )

    fdc.write( "\n\twritelog(\"%s\\n\");\n"%B )
    fdc.write( "#endif\n\n" )


def deps2( fdh, fdc, A, Bs ):
    fdh.write( "#ifdef %s\n"%A )
    fdc.write( "#ifdef %s\n"%A )
    for b in Bs:
        fdh.write( "#ifndef %s\n"%b )
        fdh.write( "#define %s\n"%b )
        fdh.write( "#endif\n" )
        fdc.write( "\twritelog(\"%s\\n\");\n"%b )
    fdh.write( "#endif\n\n" )
    fdc.write( "#endif\n\n" )


def gen_config( param_file, run_dir, new_param_file ):

    lines = open( param_file ).readlines()
    fn_h = run_dir + '/gadget-analysis-config.h'
    fn_c = run_dir + '/compile-time-info.c'
    fn_in = new_param_file
    fd_h = open( fn_h, "w" )
    fd_c = open( fn_c, "w" )
    fd_in = open( fn_in, "w" )

    fd_c.write( "#include \"allvars.h\"\n" )
    fd_c.write( "void compile_time_info (){\n" )
    fd_c.write( "\tput_header( \"compile_time_info:\" );\n" )
    
    for l in lines:
        if l[0] == '%' or len(l.strip()) == 0:
            continue
        t = [tt.strip() for tt in l.split()]
        if t[1] == 'on':
            fd_h.write( "#define %s\n"%(t[0].upper()) )
            fd_c.write( "\twritelog( \"%s\\n\" );\n"%(t[0].upper()) )
            continue
        if t[1] == 'off':
            continue
        fd_in.write( l )
    fd_in.close()

    fd_h.write( "\n\n"  )
    deps1( fd_h,fd_c, ( "GROUPTEMP", "GROUPU", "GROUPSFR", "GROUPB", "GROUPMACH",\
                "GROUPCRE", "GROUPRAD", "GROUPSPEC", "GROUPELECSPEC",\
                "GROUPTEMPPROFILE", "GROUPTEMPSTACK"), "GROUP" )
    deps1( fd_h, fd_c, ("GROUPSFR","BPDF", "PHASE"), "READSFR" )
    deps1( fd_h, fd_c, ("GROUPU","UTPFD", "GROUPCRE"), "READU" )
    deps1( fd_h, fd_c, ("GROUPB", "BPDF", "BDENSPDF", "DIVBERRPDF"), "READB" )
    deps1( fd_h, fd_c, ("DIVBERRPDF", "DIVBERRPDf", "DIVBERRDENSPDF"), "READDIVB" )
    deps1( fd_h, fd_c, ("GROUPMACH",), "READMACH" )
    deps1( fd_h, fd_c, ("GROUPCRE","GROUPELECSPEC", "CREPPDF", "CRENSLICE", "CREESLICE",\
                 "CRENTPDF" ), "READCRE" )
    deps1( fd_h, fd_c, ("GROUPTEMP","TEMPSLICE", "PDFTDIFFDENS", "PHASE","CRENTPDF",\
                    "TPDF", "GASRATIO", "HSMLTPDf", "UTPDF", "BPDF" ), "COMPUTETEMP" )
    deps1( fd_h, fd_c, ("GROUPRAD","RADSLICE", "TOTSPEC", "GROUPSPEC"), "RADSPEC" )
    deps1( fd_h, fd_c, ("HSMLTPDF","HSMLDENSPDF", "RADSLICE", "SMOOTH", "DIVBERRPDF", "DIVBERRDENSPDF"), "READHSML" )
    deps1( fd_h, fd_c, ("MF",), "FOF" )
    deps1( fd_h, fd_c, ("BSMOOTH",), "SMOOTH" )

    deps2( fd_h, fd_c, "RADSPEC", ("READB", "READCRE", "READHSML") )
    deps2( fd_h, fd_c, "GROUP", ("FOF", "TREE") )

    deps2( fd_h, fd_c, "FOF", ("TREE",) )

    deps1( fd_h, fd_c, ("FOF",), "READVEL" )

    fd_h.write( "#ifdef COMPUTETEMP\n" )
    fd_h.write( "#ifndef READTEMP\n" )
    fd_h.write( "#define READU\n" )
    fd_h.write( "#endif\n" )
    fd_h.write( "#endif\n" )
    fd_c.write( "}\n" )
    
    fd_h.close()
    fd_c.close()

def main():

    if len(sys.argv) < 2:
        print( "give parameter file" )
        exit()
    param_file = sys.argv[1]
    new_param_file = param_file + '.new'
    
    run_dir = os.path.dirname( os.path.realpath( sys.argv[0] ) )
    cur_dir = os.getcwd()
    #print( run_dir )
    #print( cur_dir )
    
    gen_config( param_file, run_dir, new_param_file )

    os.chdir( run_dir )
    gen_allvars()
    gen_add_params()
    make_protos()
    os.system( 'make delete' )
    os.system( 'make -j20' )
    os.chdir( cur_dir )

    if len(sys.argv) == 2:
        cmd = '%s/bin/gadget-analysis %s'%(run_dir, new_param_file )

    if len(sys.argv) > 2:
        a = int( sys.argv[2] )
        if len(sys.argv) == 3:
            b = 1
        else:
            b = int(sys.argv[3])
        cmd = 'mpirun -np %i %s/bin/gadget-analysis %s %i'\
            %(a, run_dir, new_param_file, b)

    print( "\nRUN: `%s`\n"%cmd )

    os.system( cmd )

main()
