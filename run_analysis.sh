#!/bin/bash
out_prefix="test_snapshot"
sed -i '1s/.*/FilePrefix  '"$out_prefix"'/' gadget-analysis.in
num=`ls ./$out_prefix* | wc -l`
echo "snapshot number: " $num
rm gadget-analysis.log -rf
mpirun -np $num ./gadget-analysis ./gadget-analysis.in

datas=" mag divB dBdt gas_rho" #rho_0 rho_1" # hge1"
for i in $datas
do
    rm $i.zip
    zip $i.zip $i -r
    rm $i -rf
done
