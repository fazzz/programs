#!/bin/python

import mdtraj as md
trr = '/data/hp150268_a/k02372/calspa/T4-Lysozyme_AlkylUrea/flexMD2/T4L_TMUW20K3M/sam_450K_00-30ns/T4L_TMUW20K3M_sam.trr'
gro = '/data/hp150268_a/k02372/calspa/T4-Lysozyme_AlkylUrea/flexMD2/T4L_TMUW20K3M/min/T4L_TMUW20K3M_min.gro'

trja = md.load(trr,top=gro)
trjb = md.load(gro)

rmsd_total = md.rmsd(trja, trjb, 0)
rmsd_bb    = md.rmsd(trja, trjb, 0)
rmsd_sc    = md.rmsd(trja, trjb, 0)
rmsd_ha    = md.rmsd(trja, trjb, 0)

for i in range(len(rmsd_total)):
    print '%f ' %(rmsd_total[i])
