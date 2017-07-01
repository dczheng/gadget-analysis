#!/bin/bash
cp process_data process_data_mag
for line in `cat ./index_z_1.txt`
do
index=`echo $line | awk -F ':' '{print$1}' `
z=`echo $line | awk -F ':' '{print$2}'`
#echo $index
#echo $z
mag_dir="./mag/"
if [ ! -d ${mag_dir} ]
then
    mkdir ${mag_dir}
fi
out_prefix="${mag_dir}${z}"
echo $out_prefix
if [ ! -d ${out_prefix} ]
then
    mkdir $out_prefix
fi
echo "FILE_PREFIX /mnt/ddnfs/data_users/xllian/Simulate_Radio/Gadget3_Simulation/209_209_20/outputs/209_209_20_${index}" > param_mag.in
echo "NUM_FILES 1" >> param_mag.in
echo "OUT_FILE ./ic.hdf5" >> param_mag.in
echo "OUT_PICTURE_PREFIX ${out_prefix}" >> param_mag.in
echo "SLICE_NUM 30" >> param_mag.in
echo "SLICE_INDEX_NUM 5" >> param_mag.in
echo "SLICE_INDEX 1 4 8 16 25" >> param_mag.in
echo "REDSHIFT ${z}">> param_mag.in
echo "PIC_XSIZE 256" >> param_mag.in
echo "PIC_YSIZE 256" >> param_mag.in
echo "3D_BOX 512 512 512" >> param_mag.in
echo "SCALAR_UNIT 1e-8" >> param_mag.in
echo "3D_AL 30.0 1.0 30.0" >> param_mag.in
echo "3D_AZ 30.0 5.0 30.0" >> param_mag.in
echo "3D_CORNER1 0.0 0.0 0.0" >> param_mag.in
echo "3D_CORNER2 20000.0 20000.0 20000.0" >> param_mag.in
nohup ./process_data_mag ./param_mag.in > "${mag_dir}/${z}.log" 2>&1 &
done
