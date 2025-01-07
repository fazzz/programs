#!/Users/yamamoriyuu/.pyenv/versions/anaconda3-4.1.0/bin/python

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

for i in range(len(phi[1])):
    for j in range(len(psi[0])):
        print('%8.3f %8.3f ' % (phi[1][i][j],psi[1][i][j]))
