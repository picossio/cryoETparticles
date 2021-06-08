from Bio import PDB
import numpy as np
from numpy import random
from scipy.spatial.transform import Rotation 
import os
import json
import math


#First let's create some functions

def center_atomic_coord (x,y,z):
    
    x, y, z = x-np.mean(x), y-np.mean(y), z-np.mean(z)
    return(x, y, z) 

def quaternion_rotation(q, x, y, z):
    
    '''
    Performs a rotation using quaternions.
    
    If it is based from a quaternion q1 = w + xi + yj + zk, then the 
    quaternion array q should be q = [x, y, z, w]
    '''
    
    Q = Rotation.from_quat(q).as_dcm()
    
    coord = np.array([x,y,z])
    rot_coord = np.dot(Q, coord)
    
    x_rot = rot_coord[0]
    y_rot = rot_coord[1]
    z_rot = rot_coord[2]
    
    return(x_rot, y_rot, z_rot)

def randomquat( ):

    #Generates uniform quaternions in SO2
    a1=random.rand()
    a2=random.rand()
    a3=random.rand()

    th=math.pi*a1
    ps=2*math.pi*a2
    ph=2*math.pi*a3

    q1=math.cos(th/2)*math.cos(ps/2)
    q2=math.cos(th/2)*math.sin(ps/2)
    q3=math.sin(th/2)*math.cos(ph+(ps/2))
    q4=math.sin(th/2)*math.sin(ph+(ps/2))

    q=[q1, q2, q3, q4]
     
    return(q)


#MAIN
#Reading BSA pdb
parser = PDB.PDBParser()
io = PDB.PDBIO()
struct = parser.get_structure('ca-4f5s.pdb','ca-4f5s.pdb')

x_atom=[]
y_atom=[]
z_atom=[]

for chains in struct:
    for chain in chains:
        for residue in chain:                             
            for atom in residue:
                x_atom.append(atom.get_vector()[0])
                y_atom.append(atom.get_vector()[1])
                z_atom.append(atom.get_vector()[2])
                
x_atom=np.array(x_atom)
y_atom=np.array(y_atom)
z_atom=np.array(z_atom)

#print(x_atom,y_atom,z_atom)


#Center the coordinates for initial rotations
x_atom, y_atom, z_atom = center_atomic_coord (x_atom, y_atom, z_atom)

x_center=[]
y_center=[]

for j in range(0,400):
      use=1
      q=randomquat()
#      print(q)

      #Rotate them
      x_atomrot, y_atomrot, z_atomrot = quaternion_rotation(q, x_atom, y_atom, z_atom)

      xc=random.randint(int(1200*1.5))
      yc=random.randint(int(1200*1.5))      

      if j>0:
           for k in range(0,len(x_center)):
                 dis=np.sqrt((xc-x_center[k])**2+(yc-y_center[k])**2) 
                 if dis<125:
                      use=0

      if use==1:
           x_center.append(xc)
           y_center.append(yc)


      #Displace them along x,y plane
      if use==1:
          print("BSA", j)
          for i in range(0,len(x_atom)):
                print(x_atomrot[i]+xc,y_atomrot[i]+yc,z_atomrot[i])

