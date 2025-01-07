#!/Users/yamamoriyuu/.pyenv/versions/anaconda3-4.1.0/bin/python

import argparse
import numpy as np
import mdtraj as md

parser = argparse.ArgumentParser(description='calc. dihedral angles of gro')
parser.add_argument('gro', metavar='GRO', type=str, help='file name of GRO')
args = parser.parse_args()

gro = args.gro

trj = md.load(gro)

atom_pairs=((12,59),(12,74),(12,88),(12,95),(12,109),(33,74),(33,88),(33,95),(33,109),(53,88),(53,95),(53,109),(59,95),(59,109),(74,109),(88,109))

atom_pairs = np.array(atom_pairs)

dist = md.compute_distances(trj,atom_pairs)

for i in range(len(dist[0])) :
    print("%8.3f " % (dist[0][i]),end="")
print(' ')
