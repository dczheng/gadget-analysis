#!/bin/bash
make clean
make 
#rm -rf png
#mkdir  png
mpirun -np 1 ./process_data ./param.in
:<<BLOCK
for i in `cat index.dat`
do
    index=`echo $i|awk -F ':' '{print$1}'`
    z=`echo $i|awk -F ':' '{print$2}'`
    echo $index $z
    sed -i '2s/^.*$/REDSHIFT '"$z"'/' param.in
    sed -i '1s/^.*$/FILE_PREFIX \/mnt\/ddnfs\/data_users\/xllian\/Simulate_Radio\/Gadget3_Simulation\/cool_feedbak_init_mag_bh_fof_1677\/outputs\/snapshot_'"$index"'/' param.in
    mpirun -np 32 ./process_data ./param.in
done
BLOCK
make clean

