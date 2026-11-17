import numpy as np
import os
file_path1 = '/Volumes/LaCie/project2/khanh_coord/pos_coordinates_1NKW.txt'
file_path2 = '/Volumes/LaCie/project2/khanh_coord/neg_coordinates_1NKW.txt'
file_path3 = '/Volumes/LaCie/project2/khanh_coord/p_coordinates_1NKW.txt'
if os.path.exists(file_path1):
    os.remove(file_path1)
if os.path.exists(file_path2):
    os.remove(file_path2)
if os.path.exists(file_path3):
    os.remove(file_path3)
cmd.load('/Volumes/LaCie/project2/mole_khanh/1NKW/pymol/complex.py')
cmd.remove('chain B')
cmd.select("arg_lys_tun_sel","resn ARG+LYS within 50 of Tunnel5")
cmd.select('carbons1', 'name CA in arg_lys_tun_sel')
xyz_pos=cmd.get_coords('carbons1',1)
cmd.select("asp_glu_tun_sel","resn ASP+GLU within 50 of Tunnel5")
cmd.select('carbons2', 'name CA in asp_glu_tun_sel')
xyz_neg=cmd.get_coords('carbons2',1)
cmd.select('p_tun_sel', 'name p within 50 of Tunnel5')
xyz_p=cmd.get_coords('p_tun_sel',1)
np.savetxt('/Volumes/LaCie/project2/khanh_coord/pos_coordinates_1NKW.txt',xyz_pos,fmt='%.2f')
np.savetxt('/Volumes/LaCie/project2/khanh_coord/neg_coordinates_1NKW.txt',xyz_neg,fmt='%.2f')
np.savetxt('/Volumes/LaCie/project2/khanh_coord/p_coordinates_1NKW.txt',xyz_p,fmt='%.2f')
cmd.create('pos', 'carbons1')
cmd.color('blue', 'pos')
cmd.create('neg', 'carbons2')
cmd.color('red', 'neg')
cmd.create('phos', 'p_tun_sel')
cmd.color('orange', 'phos')
