#!/Users/yamamoriyuu/.pyenv/versions/anaconda3-4.1.0/bin/python

import numpy as np
import argparse
import mdtraj as md

parser = argparse.ArgumentParser(description='calc. dihedral angles of traj.')
parser.add_argument('xtc', metavar='XTC', type=str, help='file name of XTC')
parser.add_argument('gro', metavar='GRO', type=str, help='file name of GRO')
args = parser.parse_args()

trr = args.xtc
gro = args.gro

trj = md.load(trr,top=gro)

phi = md.compute_phi(trj)
psi = md.compute_psi(trj)

k = 0
for i in range(len(phi[1])):
    for j in range(len(phi[0])):
        sinphi = np.sin(phi[1][i][j])
        cosphi = np.cos(phi[1][i][j])
        sinpsi = np.sin(psi[1][i][j])
        cospsi = np.cos(psi[1][i][j])
        k = k + 1
        print('%8.3f ' % (sinphi),end="")
        if (k%5 == 0):
            print(' ')
        k = k + 1
        print('%8.3f ' % (cosphi),end="")
        if (k%5 == 0):
            print(' ')
        k = k + 1
        print('%8.3f ' % (sinpsi),end="")
        if (k%5 == 0):
            print(' ')
        k = k + 1
        print('%8.3f ' % (cospsi),end="")
        if (k%5 == 0):
            print(' ')
print(' ')
