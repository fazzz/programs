
libs=( dummy ABAb_backpass ABAb_calcattfrc ABAb_calcKineE  ABAb_gtree  ABAb  ABAb_integ ABAb_Inverse_backpass ABAb_Inverse  ABAb_Inverse_mainpass ABAb_mainpass ABAb_Nose-Hoover_chain ABAb_Nose-Hoover ABAb_Nose-Hoover_new ABAb_Nose-Hoover_new_mvV ABAb_pick_data  ABAb_prepass  ABAb_set_frc ABAb_set_imatb ABAb_set_imat ABAb_set_lref ABAb_set_rst ABAb_set_tmat ABAb_set_trans ABAb_update )

nlibs=`expr ${#libs[*]} - 1`

cd ~/work/programs/ABAb
for i in `seq 1 ${nlibs}`; do
    sh  /home/yamamori/work/programs/compile/makeLIB.sh ${libs[$i]} ${libs[$i]}
done   

