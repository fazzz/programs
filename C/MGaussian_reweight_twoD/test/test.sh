#!/bin/sh

dir=/home/yamamori/calspa/refcalc/UmbSam/AD/s_UmbSam_vac_2012-11-12_ff99SB/UmbAtPhsi_0_30_60_90_120_150_180_210_240/anl

r3phi_WHAM=( dummy 0   0.52360  1.0472 1.5708 2.0944  2.6180 3.1416 3.6652 4.1888 4.7124 5.2360 5.7596 )
r3psi_WHAM=( dummy 0   0.52360  1.0472 1.5708 2.0944  2.6180 3.1416 3.6652 4.1888 4.7124 5.2360 5.7596 )

rk3phi_WHAM=( dummy 20 20 20 20 20 20 20 20 20 20 20 20 )
rk3psi_WHAM=( dummy 20 20 20 20 20 20 20 20 20 20 20 20 )

MGRin=input/MGR.in
UMBin=input/umb.in

rm ${MGRin}
rm ${UMBin}

for i in `seq 1 12`; do
    for j in `seq 1 12`; do
	dtrj=${dir}/ADv_${i}_${j}_10ns.dtrj
	dtrjCV=${dir}/ADv_${i}_${j}_10ns.CV
	gawk '{if(NR%100==1){print $1 " " $2 " " $3}}' < ${dtrj} > temp.txt
	sed -e "1d" temp.txt > ${dtrjCV}
	rm temp.txt

	cat <<EOF >> ${MGRin}
${dtrjCV}
EOF

	cat <<EOF >> ${UMBin}
${r3phi_WHAM[$i]} ${rk3phi_WHAM[$i]} ${r3psi_WHAM[$j]} ${rk3psi_WHAM[$j]}
EOF

    done
done

#../bin/MGaussian_reweight_twoD 144 10 ${MGRin} ${UMBin} output/MGaussian12x12.txt output/pmf12x12.txt

#set env LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/yamamori/work/programs/ABAMD_2014_05_01/lib
#cd /home/yamamori/work/programs/MGaussian_reweight_twoD/test
#run 144 10 input/MGR.in input/umb.in output/MGaussian12x12.txt output/pmf12x12.txt
