#!/bin/bash
for line in `cat ./index_z_1.txt`
do
index=`echo $line | awk -F ':' '{print$1}' `
z=`echo $line | awk -F ':' '{print$2}'`
echo "FILE_PREFIX /mnt/ddnfs/data_users/xllian/Simulate_Radio/Gadget3_Simulation/209_209_20/outputs/209_209_20_${index}" > param.in
echo "NUM_FILES 1" >> param.in
echo "OUT_FILE ./ic.hdf5" >> param.in
echo "OUT_PICTURE_PREFIX ${out_prefix}" >> param.in
echo "SLICE_NUM 30" >> param.in
echo "SLICE_INDEX_NUM 5" >> param.in
echo "SLICE_INDEX 1 4 8 16 25" >> param.in
echo "REDSHIFT ${z}">> param.in
echo "PIC_XSIZE 256" >> param.in
echo "PIC_YSIZE 256" >> param.in
echo "3D_BOX 512 512 512" >> param.in
echo "SCALAR_UNIT 1e-8" >> param.in
echo "3D_AL 30.0 1.0 30.0" >> param.in
echo "3D_AZ 0.0 5.0 360.0" >> param.in
echo "3D_CORNER1 0.0 0.0 0.0" >> param.in
echo "3D_CORNER2 20000.0 20000.0 20000.0" >> param.in
./process_data ./param.in 
done
