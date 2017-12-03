#!/usr/bin/env python3
import sys
import os

pic_path = sys.argv[1]
gif_name = sys.argv[2]
pic_name = os.listdir( pic_path )
z = []
for n in pic_name:
    tmp = n.split( '_' )
    z.append( float(tmp[-1][:-4]) )
z_num = len( z )
for i in range( z_num-1 ):
    for j in  range( i, z_num ):
        if ( z[i] < z[j] ):
            tmp = z[i]
            z[i] = z[j]
            z[j] = tmp
            tmp = pic_name[i]
            pic_name[i] = pic_name[j]
            pic_name[j] = tmp

png_list = ' '

for i in range( z_num ):
    pic_name[i] = pic_path + '/' + pic_name[i]

png_list = png_list.join( pic_name )
convert_cmd = 'convert -delay 40  -loop 10 ' + png_list + ' ' + gif_name
print( convert_cmd )
os.system( convert_cmd )
