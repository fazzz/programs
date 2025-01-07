#!/Users/yamamoriyuu/.pyenv/versions/anaconda3-4.1.0/bin/python

import argparse
import numpy as np
import mdtraj as md

parser = argparse.ArgumentParser(description='calc. Distances of Salt-Bridge')
parser.add_argument('xtc', metavar='XTC', type=str, help='file name of XTC')
parser.add_argument('gro', metavar='GRO', type=str, help='file name of GRO')
args = parser.parse_args()

xtc = args.xtc
gro = args.gro

t = md.load(xtc, top=gro)

atoms_01 = t.top.select("(resid == 72) and (name =~ 'OD[12]')")
atoms_02 = t.top.select("(resid == 74) and ((name =~ 'NE') or (name =~ 'NH[12]'))")

indices = []

for i in atoms_01:
    for j in atoms_02:
        indices.append([i,j])
        print('%s-' %(t.top.atom(i)), end="")
        print('%s ' %(t.top.atom(j)), end="")
        
dist = md.compute_distances(t,indices)

for i in range(len(dist)):
    for j in range(len(indices)):
        print('%8.3f ' % (dist[i][j]) ,end="")
    print(' ')

