#!/bin/sh

########################################
# if [ -f output/all.scoverAll ]; then #
#     rm output/all.scoverAll	       #
# fi				       #
# 				       #
# touch output/all.scoverAll	       #
########################################

for f in `ls input/*.qscore.txt`; do
    name=$( basename ${f} .qscore.txt )											      
    # n=$( wc -l ${f} | awk '{print $1}' )											        #
    # 																        #
    # grep -E "^${name}\." scoverAll/all.scoverAll > output/${name}.scoverAll							        #
    # 																        #
    # for d in 0.80 0.85 0.90; do												        #
    # 	./clust_simple ${n} ${d} ${f} > output/${name}.qscore.clust.th${d}							        #
    # 																        #
    # 	paste output/${name}.scoverAll output/${name}.qscore.clust.th${d} > temp						        #
    # 	mv temp output/${name}.scoverAll											        #
    # done															        #
    # 																        #
    # printf "${name} " >> output/all.scoverAll											        #
    # sort -nrk 2,2 output/${name}.scoverAll | head -1 | awk '{printf("%8.3f %8.3f %8.3f ",$2,$3,$4)}' >> output/all.scoverAll	        #
    # sort -nrk 9,9 output/${name}.scoverAll | head -1 | awk '{printf("%8.3f %8.3f %8.3f ",$2,$3,$4)}' >> output/all.scoverAll	        #
    # sort -nrk 10,10 output/${name}.scoverAll | head -1 | awk '{printf("%8.3f %8.3f %8.3f ",$2,$3,$4)}' >> output/all.scoverAll        #
    # sort -nrk 11,11 output/${name}.scoverAll | head -1 | awk '{printf("%8.3f %8.3f %8.3f\n",$2,$3,$4)}' >> output/all.scoverAll       #
    
    
    ./plot_score_vs_score.sh output/${name}.scoverAll
done

