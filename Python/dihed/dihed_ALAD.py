#!/Users/yamamoriyuu/.pyenv/versions/anaconda3-4.1.0/bin/python

import argparse
import numpy as np
import mdtraj as md

parser = argparse.ArgumentParser(description='calc. dihedral angles of traj.')
parser.add_argument('xtc', metavar='XTC', type=str, help='file name of XTC')
parser.add_argument('gro', metavar='GRO', type=str, help='file name of GRO')
args = parser.parse_args()

trr = args.xtc
gro = args.gro

trj = md.load(trr,top=gro)

indices = [[1,2,5,7],[2,5,7,9],[5,7,9,15],[7,9,11,12],[7,9,15,17],[9,15,17,19],[15,17,19,20]]
dih = md.compute_dihedrals(trj,indices)

##################################################
# k = 0                                          #
# for i in range(len(dih)):                      #
#      for j in range(len(indices)):             #
#           k = k + 1                            #
#           print('%8.3f ' % (dih[i][j]),end="") #
#           if (k%5 == 0):                       #
#                print(' ')                      #
# print(' ')                                     #
##################################################

for i in range(len(dih)):
     for j in range(len(indices)):
          print('%8.3f ' % (dih[i][j]),end="")
     print(' ')  
     
