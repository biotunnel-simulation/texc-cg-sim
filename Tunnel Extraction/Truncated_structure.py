import numpy as np
import os

cmd.load("/Users/6WD4.cif")

cmd.select('ul4', 'chain d')
cmd.select('ul22', 'chain s')
cmd.select('ul23', 'chain t')

cmd.create('ul4', 'ul4')
cmd.show('spheres', 'ul4')
cmd.color('red', 'ul4')

cmd.create('ul22', 'ul22')
cmd.show('spheres', 'ul22')
cmd.color('yellow', 'ul22')

cmd.create('ul23', 'ul23')
cmd.show('spheres', 'ul23')
cmd.color('blue', 'ul23')

cmd.select('nucleotide_61', 'resi 61 and chain t')
cmd.create('U61', 'nucleotide_61')
coords = cmd.get_coords('nucleotide_61')
file_path = '/Users/exit_coords_6WD4.txt'
if os.path.exists(file_path):
   os.remove(file_path)
np.savetxt('/Users/exit_coords_6WD4.txt', coords[0])

obj1_name = "ul4"
obj2_name = "ul22"
cutoff = 15.0
cmd.select("close_atoms", f"{obj1_name} within {cutoff} of {obj2_name}")
cmd.create('closest', 'close_atoms')

atom_list = cmd.get_model("closest").atom
atoms = cmd.get_model("closest").atom
last_atom = "closest and id " + str(atom_list[-1].id)
cmd.select("last_atom", last_atom)
cmd.create('last_a', 'last_atom')

cmd.select("edited", "br. 6WD4 within 80 of last_a")
cmd.create("6WD4_edited", "edited")
cmd.show('spheres', '6WD4_edited')
cmd.save("/Users/6WD4_edited.pdb", "6WD4_edited")
