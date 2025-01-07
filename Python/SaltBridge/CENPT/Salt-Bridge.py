#!/Users/yamamoriyuu/.pyenv/versions/anaconda3-4.1.0/bin/python

import argparse
import numpy as np
import mdtraj as md

#parser = argparse.ArgumentParser(description='calc. Distances of Salt-Bridge')
parser = argparse.ArgumentParser(description='This script calculates distance of Salt-Bridge between Resdiue i and j')
parser.add_argument('xtc', metavar='XTC', type=str, help='file name of XTC')
parser.add_argument('gro', metavar='GRO', type=str, help='file name of GRO')
parser.add_argument('i_ares', metavar='ARES', type=int, help='Acid Residue Num.')
parser.add_argument('j_bres', metavar='BRES', type=str, help='Basic Residue Num.')
args = parser.parse_args()

xtc = args.xtc
gro = args.gro
i_ares = args.i_ares
j_bres = args.j_bres

#print(i_ares)
#print(j_bres)

t = md.load(xtc, top=gro)

atoms_01 = t.top.select("(resid == i_ares) and (name =~ '^O') and sidechain")
atoms_02 = t.top.select("(resid == j_bres) and (name =~ '^N') and sidechain")

print(atoms_01)
print(atoms_02)

#indices = []

#for i in atoms_01:
#    for j in atoms_02:
#        indices.append([i,j])

#dist = md.compute_distances(t,indices)

#for i in range(len(dist)):
#    sbdist = min(dist[i])
#    sbindx = np.where(dist[i] == sbdist)
#    sbindx_atom01 = indices[(sbindx[0])[0]][0]
#    sbindx_atom02 = indices[(sbindx[0])[0]][1]
#    print('%8.3f %s--%s' % (sbdist*10.0, t.top.atom(sbindx_atom01), t.top.atom(sbindx_atom02)))
