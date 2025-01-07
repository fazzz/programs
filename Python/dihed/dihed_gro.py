#!/Users/yamamoriyuu/.pyenv/versions/anaconda3-4.1.0/bin/python

import argparse
import mdtraj as md

parser = argparse.ArgumentParser(description='calc. dihedral angles of gro')
parser.add_argument('gro', metavar='GRO', type=str, help='file name of GRO')
args = parser.parse_args()

gro = args.gro

trj = md.load(gro)

phi = md.compute_phi(trj)
psi = md.compute_psi(trj)

print(phi[0])

print(psi[0])

for i in range(len(phi[1])):
    for j in range(len(psi[0])):
        print('%8.3f %8.3f ' % (phi[1][i][j],psi[1][i][j]))
