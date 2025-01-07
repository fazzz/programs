#!~/bin/sh

cat <<EOF > input/MGaussian_reweight.in
/home/yamamori/calspa/refcalc/UmbSam/AD_oneD/s_UmbSam_vac_2014-07-19_ff99SB/Umb12_alPhi/anl/ADv_1_10ns.CV
/home/yamamori/calspa/refcalc/UmbSam/AD_oneD/s_UmbSam_vac_2014-07-19_ff99SB/Umb12_alPhi/anl/ADv_2_10ns.CV
/home/yamamori/calspa/refcalc/UmbSam/AD_oneD/s_UmbSam_vac_2014-07-19_ff99SB/Umb12_alPhi/anl/ADv_3_10ns.CV
/home/yamamori/calspa/refcalc/UmbSam/AD_oneD/s_UmbSam_vac_2014-07-19_ff99SB/Umb12_alPhi/anl/ADv_4_10ns.CV
/home/yamamori/calspa/refcalc/UmbSam/AD_oneD/s_UmbSam_vac_2014-07-19_ff99SB/Umb12_alPhi/anl/ADv_5_10ns.CV
/home/yamamori/calspa/refcalc/UmbSam/AD_oneD/s_UmbSam_vac_2014-07-19_ff99SB/Umb12_alPhi/anl/ADv_6_10ns.CV
/home/yamamori/calspa/refcalc/UmbSam/AD_oneD/s_UmbSam_vac_2014-07-19_ff99SB/Umb12_alPhi/anl/ADv_7_10ns.CV
/home/yamamori/calspa/refcalc/UmbSam/AD_oneD/s_UmbSam_vac_2014-07-19_ff99SB/Umb12_alPhi/anl/ADv_8_10ns.CV
/home/yamamori/calspa/refcalc/UmbSam/AD_oneD/s_UmbSam_vac_2014-07-19_ff99SB/Umb12_alPhi/anl/ADv_9_10ns.CV
EOF

cat <<EOF > input/umb.in
    0       10     
    0.52360 10
    1.0472  10
    1.5708  10
    2.0944  10
    2.6180  10
    3.1416  10
    3.6652  10
    4.1888  10
EOF

../bin/MGaussian_reweight 9 10 input/MGaussian_reweight.in input/umb.in output/MGaussian.txt output/pmf.txt
