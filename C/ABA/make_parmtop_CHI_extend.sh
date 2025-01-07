 #!/bin/sh

pn=CHI
seq=( NGLY TYR ASP PRO GLU THR GLY THR TRP CGLY  )

dir=$(pwd)
								
leap=/home/appl/amber10/exe/tleap 			
pn=${pn}v

cat <<EOF > ${pn}_extend.cmd
pn=sequence{ ${seq[*]} }
saveAmberparm pn ${pn}_extend.top  ${pn}_extend.crd
quit
EOF

cat <<EOF > ${pn}_extend.seq
${seq[*]}
EOF

export AMBERHOME=/home/appl/amber10
${leap} -s -f //home/appl/amber10/dat/leap/cmd/leaprc.ff99SB -f ${dir}/${pn}_extend.cmd

cp ${pn}_extend.top  ${pn}_extend.crd ${pn}_extend.cmd ${INPUT}
