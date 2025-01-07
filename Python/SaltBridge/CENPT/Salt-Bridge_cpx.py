#!/Users/yamamoriyuu/.pyenv/versions/anaconda3-4.1.0/bin/python

import argparse
import numpy as np
import mdtraj as md

parser = argparse.ArgumentParser(description='This script calculates distance of Salt-Bridge of Protein')
parser.add_argument('xtc', metavar='XTC', type=str, help='file name of XTC')
parser.add_argument('gro', metavar='GRO', type=str, help='file name of GRO')
args = parser.parse_args()

t = md.load(args.xtc, top=args.gro)

indices = []

atoms_01 = t.top.select("(resname == ASP or resname == GLU or resname == ASN or resname == GLN or resname == T1P or resname == T2P or resname == S1P or resname == S2P) and (name =~ '^O.{1,}')" )
atoms_02 = t.top.select("(resname == ARG or resname == LYS or resname == HIS) and (name =~ '^N.{1,}')")

k = 0
for i in t.top.residues:
    for j in t.top.residues:
        indx_res_i = t.top.atom(i).residue.index
        indx_res_j = t.top.atom(j).residue.index
        indices.append([i,j])

dist = md.compute_distances(t,indices)
        
for i in range(len(dist)):
    sbdist = min(dist[i])
    sbindx = np.where(dist[i] == sbdist)
    sbindx_atom01 = indices[(sbindx[0])[0]][0]
    sbindx_atom02 = indices[(sbindx[0])[0]][1]
    print('%8.3f %s--%s' % (sbdist*10.0, t.top.atom(sbindx_atom01), t.top.atom(sbindx_atom02)))

