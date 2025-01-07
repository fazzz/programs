#!/bin/sh

if [ -f output_2/all.scoverAll ]; then
    rm output_2/all.scoverAll
fi

touch output_2/all.scoverAll

for f in `ls input_2/*.qscore.txt`; do
    name=$( basename ${f} .qscore.txt )
    n=$( wc -l ${f} | awk '{print $1}' )

    grep -E "^${name}\." scoverAll/all.scoverAll > output_2/${name}.scoverAll
    
    for d in 0.80 0.85 0.90; do
	./clust_simple ${n} ${d} ${f} > output_2/${name}.qscore.clust.th${d}

	paste output_2/${name}.scoverAll output_2/${name}.qscore.clust.th${d} > temp
	mv temp output_2/${name}.scoverAll
    done

    printf "${name} " >> output_2/all.scoverAll
    sort -nrk 2,2 output_2/${name}.scoverAll | head -1 | awk '{printf("%8.3f %8.3f %8.3f ",$2,$3,$4)}' >> output_2/all.scoverAll
    sort -nrk 9,9 output_2/${name}.scoverAll | head -1 | awk '{printf("%8.3f %8.3f %8.3f ",$2,$3,$4)}' >> output_2/all.scoverAll
    sort -nrk 10,10 output_2/${name}.scoverAll | head -1 | awk '{printf("%8.3f %8.3f %8.3f ",$2,$3,$4)}' >> output_2/all.scoverAll
    sort -nrk 11,11 output_2/${name}.scoverAll | head -1 | awk '{printf("%8.3f %8.3f %8.3f\n",$2,$3,$4)}' >> output_2/all.scoverAll
    
    ./plot_score_vs_score.sh output_2/${name}.scoverAll
done

