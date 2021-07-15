import numpy as np
from numpy import random 
from Bio import PDB

from atom_coord import center_atomic_coord
from atom_coord import quaternion_rotation
from atom_coord import randomquat

#MAIN
#Reading BSA pdb
parser = PDB.PDBParser()    #for reading pdb files
io = PDB.PDBIO()    #for writing pdb files
struct = parser.get_structure('BSA','/Users/dkleebatt/Documents/Structures/4f5s.pdb') #retrieving structure from specified pdb file

num_protein = 0 #initializing counter

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

#print(x_atom,y_atom,z_atom)


#Center the coordinates for initial rotations
x_atom, y_atom, z_atom = center_atomic_coord (x_atom, y_atom, z_atom) #puts BSA at the origin

x_center = []
y_center = []

atom_pos = []
x_pos = []
y_pos = []
z_pos = []

for j in range(0,1):        # in range of 0 to 400
    use=1
    q=randomquat()          #q is the output of the uniformly generated quaternion generator
    #print(q)

    #Rotate them
    x_atomrot, y_atomrot, z_atomrot = quaternion_rotation(q, x_atom, y_atom, z_atom)

    xc=random.randint(int(1200*1.5)) # int makes an interger from a digit... randint generates uniform integers [0,1800]
    yc=random.randint(int(1200*1.5))      # I think this is the x and y coordinate where the center of the protein is being placed

    if j>0:
          for k in range(0,len(x_center)):   #for all k between 0 and the number of objects in x_center
                dis=np.sqrt((xc-x_center[k])**2+(yc-y_center[k])**2) ## distance calculation between proposed position and all other positions
                if dis<100:  # volume exclusion... if this distance is less than 125... don't place the protein here 
                    use=0   # assign use the value of 0

    if use==1:      # when j=0, append the x,y_center arrays with the generated random integer 
          x_center.append(xc)
          y_center.append(yc)

    #Displace them along x,y plane
    if use==1:
        #print("BSA", j)
        num_protein = num_protein + 1
        for i in range(0,len(x_atom)): #for each atom 
            #print(x_atomrot[i]+xc,y_atomrot[i]+yc,z_atomrot[i])
            x_pos.append(x_atomrot[i]+xc)
            y_pos.append(y_atomrot[i]+yc)
            z_pos.append(z_atomrot[i])
x_pos = np.array(x_pos)
y_pos = np.array(y_pos)
z_pos = np.array(z_pos)

atom_pos = np.vstack([x_pos,y_pos,z_pos])
print(np.shape(atom_pos))
print(str(num_protein) + ' proteins')
print(len(x_pos))
print(atom_pos)
