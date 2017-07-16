#!/bin/bash
make clean
make 
#rm -rf png
#mkdir  png
mpirun -np 200 ./process_data ./param.in
:<<BLOCK
for i in `seq  0 5 360`
do
    echo $i
    sed -i 's/3D_AZ.*/3D_AZ '"$i"'/' params.in
    ./process_data ./params.in
done
BLOCK
make clean
