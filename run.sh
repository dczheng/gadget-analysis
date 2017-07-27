#!/bin/bash
make clean
make 
#rm -rf png
#mkdir  png
#mpirun -np 73 ./gadget_analysis ./param.in
#:<<BLOCK
for i in `cat index.dat`
do
    index=`echo $i|awk -F ':' '{print$1}'`
    z=`echo $i|awk -F ':' '{print$2}'`
    echo $index $z
    sed -i '2s/^.*$/REDSHIFT '"$z"'/' param.in
    sed -i '1s/^.*$/FILE_PREFIX \/mnt\/ddnfs\/data_users\/xllian\/Simulate_Radio\/Gadget3_Simulation\/cool_feedbak_init_mag_bh_fof_1677\/outputs\/snapshot_'"$index"'/' param.in
    nohup ./gadget_analysis ./param.in  >> ${index}_log  2>&1 &
    sleep 30
done
#BLOCK
make clean

