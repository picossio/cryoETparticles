import numpy as np
import matplotlib.pyplot as plt
from script_atom_coord import atom_pos
from script_atom_coord import x_pos

class img_generator:
    def __init__(self, atom_pos, x_pos):

        self.atom_pos = atom_pos
        self.x_pos = x_pos

        self.align_system

    def align_system(self):
        center_pos = np.sum(atom_pos, axis = 1) / len(x_pos)
        aligned_pos = np.zeros((3,len(x_pos)))

        for k in range(len(x_pos)):
            aligned_pos[0,k] = atom_pos[0,k] - center_pos[0]
            aligned_pos[1,k] = atom_pos[1,k] - center_pos[1]
            aligned_pos[2,k] = atom_pos[2,k] - center_pos[2]
        
        self.aligned_pos = aligned_pos

    def generate_image(self, n_pixels, pixel_size, sigma, cutoff):
        
        self.system_atoms = self.aligned_pos
        grid_min = -pixel_size * (n_pixels + 1)*0.5
        grid_max = pixel_size * (n_pixels - 3)*0.5 + pixel_size

        x_grid = np.arange(grid_min, grid_max, pixel_size)
        y_grid = np.arange(grid_min, grid_max, pixel_size) # this is the pixel grid 

        Ixy = np.zeros((n_pixels, n_pixels))

        for atom in range(len(x_pos)): # for each atom in the array 
        
            #Values of the gaussians
            gauss_x = np.zeros((n_pixels)) #initializing array for gaussian
            gauss_y = np.zeros((n_pixels))
            
            #Selected neighbors
            x_neigh = np.where(np.absolute(x_grid - self.system_atoms[0,:][atom]) <= cutoff*sigma) # select neighbor pixels of a given atom when their distance <= cutoff*sigma
            y_neigh = np.where(np.absolute(y_grid - self.system_atoms[1,:][atom]) <= cutoff*sigma) #        
                
            gauss_x[x_neigh] = np.exp( -0.5 * ( ((x_grid[x_neigh] - self.system_atoms[0,:][atom]) /sigma)**2) ) # generates gaussian @all neighboring atoms depending on their distance from the given coordinate
            gauss_y[y_neigh] = np.exp( -0.5 * ( ((y_grid[y_neigh] - self.system_atoms[1,:][atom]) /sigma)**2) )
            
            Ixy = np.add(Ixy, gauss_x[:, None]*gauss_y)  # recursively adds all calculated gaussians 

        self.image = np.sqrt(2*np.pi) * sigma * Ixy
        

    def plot_image(self):

        fig, ax = plt.subplots()
        ax.imshow(self.image)
        plt.show()

        print(self.image.shape)
        print(self.image.dtype)
        print('printed image array')
        print(self.image)
        print('printed positions')
        print(self.system_atoms)
        print(len(x_pos))
        print(np.shape(x_pos))
        