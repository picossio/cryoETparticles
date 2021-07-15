import numpy as np
from numpy import random 
from scipy.spatial.transform import Rotation 
import math 

def center_atomic_coord (x,y,z):
    
    x, y, z = x-np.mean(x), y-np.mean(y), z-np.mean(z) #Subtracts the average of all x, y and z atomic positions from the array with all x,y,z positions
    return(x, y, z)                                     # effectively, average position becomes 0,0,0... all coordinates are centered around the average being 0

def quaternion_rotation(q, x, y, z):
    
    '''
    Performs a rotation using quaternions.
    
    If it is based from a quaternion q1 = w + xi + yj + zk, then the 
    quaternion array q should be q = [x, y, z, w]
    '''
    
    Q = Rotation.from_quat(q).as_matrix() #initializes quanternions
    
    coord = np.array([x,y,z]) 
    rot_coord = np.dot(Q, coord)
    
    x_rot = rot_coord[0]
    y_rot = rot_coord[1]
    z_rot = rot_coord[2]
    
    return(x_rot, y_rot, z_rot)

def randomquat( ):

    #Generates uniform quaternions in SO3
    a1=random.rand() #uniformly distributed number between 0 and 1
    a2=random.rand()
    a3=random.rand()

    th=math.pi*a1
    ps=2*math.pi*a2 #constant pi * a random number
    ph=2*math.pi*a3

    q1=math.cos(th/2)*math.cos(ps/2)    # some constants generated from random number but filtered thru sin and cosine for some reason
    q2=math.cos(th/2)*math.sin(ps/2)
    q3=math.sin(th/2)*math.cos(ph+(ps/2))
    q4=math.sin(th/2)*math.sin(ph+(ps/2))

    q=[q1, q2, q3, q4]      # this is the output of a quaternion... is q1 or q4 the integer????
     
    return(q)

