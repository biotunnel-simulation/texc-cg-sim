import numpy as np
import os

cmd.load('/Users/6WD4/pymol/complex.py')

# If Tunnel1 is the correct tunnel after comparison
xyz=cmd.get_coords('Tunnel1',1)

r=[]
cmd.iterate_state(1,'Tunnel1','r.append(vdw)',space=locals(),atomic=0)

np.savetxt('/Users/tunnel_geo/tunnel_coordinates_6WD4.txt',xyz,fmt='\%.2f')
np.savetxt('/Users/tunnel_geo/tunnel_radius_6WD4.txt',r,fmt='\%.2f')
