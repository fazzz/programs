 #!/bin/sh

pn=GLYD
seq=( ACE GLY NME )

dir=$(pwd)
								
leap=/home/appl/amber10/exe/tleap 			
pn=${pn}v

cat <<EOF > ${pn}.cmd
pn=sequence{ ${seq[*]} }
impose pn { { 1 3 } } { { \$N \$CA \$C \$N -40.0 } { \$C \$N \$CA \$C -60.0  } }
saveAmberparm pn ${pn}_helix.top ${pn}_helix.crd
quit
EOF

cat <<EOF > ${pn}.seq
${seq[*]}
EOF

cat << eof2 > qsub_leap.sh
#!/bin/csh
#$ -S /bin/sh 
#$ -cwd 
#$ -V
#$ -N leap
#$ -q all.q
#$ -l hostname=ajisai01
#$ -l hostname=ajisai03
#$ -e ${dir}/leap.e
#$ -o ${dir}/leap.o
 
cd ${dir}
export AMBERHOME=/home/appl/amber10
${leap} -s -f //home/appl/amber10/dat/leap/cmd/leaprc.ff99SB -f ${dir}/${pn}.cmd
eof2

cat qsub_leap.sh
chmod +x qsub_leap.sh
qsub qsub_leap.sh

