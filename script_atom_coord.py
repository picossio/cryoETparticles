import numpy as np
from numpy import random 
from Bio import PDB

from atom_coord import center_atomic_coord
from atom_coord import quaternion_rotation
from atom_coord import randomquat

from numpy.random import default_rng
rng = default_rng()

parser = PDB.PDBParser()    #for reading pdb files
io = PDB.PDBIO()    #for writing pdb files
struct = parser.get_structure('BSA','/Users/dkleebatt/Documents/Structures/4f5s.pdb') #retrieving structure from specified pdb file

n_pixels = 128     
pixel_size = 1.25   
target = 4

x_atom=[]   #initializing array to hold x,y,z position of all atoms
y_atom=[]
z_atom=[]

for chains in struct:
    for chain in chains:
        for residue in chain:                           
            for atom in residue:    #for loop to locate each atom sequentially in the structure 
                x_atom.append(atom.get_vector()[0])   #append position array with atom position data in pdb file 
                y_atom.append(atom.get_vector()[1])
                z_atom.append(atom.get_vector()[2])
                
x_atom=np.array(x_atom) #x_atom is a numpy array with the x-coordinate of each atom
y_atom=np.array(y_atom)
z_atom=np.array(z_atom)

#Center the coordinates for initial rotations
x_atom, y_atom, z_atom = center_atomic_coord (x_atom, y_atom, z_atom) #puts BSA at the origin

x_center = []
y_center = []

atom_pos = []
x_pos = []
y_pos = []
z_pos = []

num_protein = 0 #initializing counter

while num_protein <= target - 1:
    use=1
    q=randomquat()          #q is the output of the uniformly generated quaternion generator
    
    #Rotate them
    x_atomrot, y_atomrot, z_atomrot = quaternion_rotation(q, x_atom, y_atom, z_atom)

    xc=random.randint(int(n_pixels * pixel_size)) 
    yc=random.randint(int(n_pixels * pixel_size)) 

    if num_protein>0:
          for k in range(0,len(x_center)):   #for all k between 0 and the number of objects in x_center
                dis=np.sqrt((xc-x_center[k])**2+(yc-y_center[k])**2) ## distance calculation between proposed position and all other positions
                if 20 < dis < 125:
                    rand_val = rng.standard_normal(1)
                    if rand_val * dis < 20:
                        use=0
                if dis<20:  # volume exclusion... if this distance is less than 125... don't place the protein here # 125 x 40 x 40 angstrom
                    use=0   # assign use the value of 0

    if use==1:      # distance condition is met, append the x,y_center arrays with the generated random integer 
          x_center.append(xc)
          y_center.append(yc)

    #Displace them along x,y plane
    if use==1:
        #print("BSA", j)
        num_protein += 1
        for i in range(0,len(x_atom)): #for each atom 
            x_pos.append(x_atomrot[i]+xc)       # x_pos contains atomic coordinates after rotation and displacement according to center position
            y_pos.append(y_atomrot[i]+yc)
            z_pos.append(z_atomrot[i])
x_pos = np.array(x_pos)
y_pos = np.array(y_pos)
z_pos = np.array(z_pos)

atom_pos = np.vstack([x_pos,y_pos,z_pos])
print(np.shape(atom_pos))
print(str(num_protein) + ' proteins')
