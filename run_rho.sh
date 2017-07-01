#!/bin/bash
cp process_data process_data_rho
for line in `cat ./index_z_1.txt`
do
index=`echo $line | awk -F ':' '{print$1}' `
z=`echo $line | awk -F ':' '{print$2}'`
#echo $index
#echo $z
rho_dir="./rho/"
if [ ! -d ${rho_dir} ]
then
    mkdir ${rho_dir}
fi
out_prefix="${rho_dir}${z}"
echo $out_prefix
if [ ! -d ${out_prefix} ]
then
    mkdir $out_prefix
fi
echo "FILE_PREFIX /mnt/ddnfs/data_users/xllian/Simulate_Radio/Gadget3_Simulation/209_209_20/outputs/209_209_20_${index}" > param_rho.in
echo "NUM_FILES 1" >> param_rho.in
echo "OUT_FILE ./ic.hdf5" >> param_rho.in
echo "OUT_PICTURE_PREFIX ${out_prefix}" >> param_rho.in
echo "SLICE_NUM 30" >> param_rho.in
echo "SLICE_INDEX_NUM 5" >> param_rho.in
echo "SLICE_INDEX 1 4 8 16 25" >> param_rho.in
echo "REDSHIFT ${z}">> param_rho.in
echo "PIC_XSIZE 256" >> param_rho.in
echo "PIC_YSIZE 256" >> param_rho.in
echo "3D_BOX 512 512 512" >> param_rho.in
echo "SCALAR_UNIT 1e-9" >> param_rho.in
echo "3D_AL 30.0 1.0 30.0" >> param_rho.in
echo "3D_AZ 0.0 5.0 360.0" >> param_rho.in
echo "3D_CORNER1 0.0 0.0 0.0" >> param_rho.in
echo "3D_CORNER2 20000.0 20000.0 20000.0" >> param_rho.in
nohup ./process_data_rho ./param_rho.in > "${rho_dir}/${z}.log" 2>&1 &
done
