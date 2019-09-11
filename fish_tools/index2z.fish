#!/usr/bin/env fish

set N ( ls snap* | wc -l )
set N (math $N-1)
echo $N
printf "" > index2z.dat
for i in ( seq 0 $N )
    set fn (printf "snapshot_%03i.hdf5" $i)
    set z (gadget-header $fn 3 1 | grep "redshift" | awk -F ':' '{print $2}')
    printf "index: %03i, z: %f\n" $i $z
    printf "%03i %f\n" $i $z >> index2z.dat
end
