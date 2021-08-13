import numpy as np
import matplotlib.pyplot as plt

class img_generator:
    def __init__(self, atom_pos, x_pos, num_protein):

        self.atom_pos = atom_pos
        self.x_pos = x_pos
        self.num_protein = num_protein
        self.align_system

    def align_system(self):

        center_pos = np.sum(self.atom_pos, axis = 1) / len(self.x_pos)
        self.center_pos = np.transpose(np.reshape(center_pos,(1,3)))
        aligned_pos = self.atom_pos[:,:] - self.center_pos[:,:]
        self.aligned_pos = aligned_pos
        
    def generate_image(self, n_pixels, pixel_size, sigma, cutoff):
        
        self.system_atoms = self.aligned_pos
        grid_min = -pixel_size * (n_pixels + 1)*0.5
        grid_max = pixel_size * (n_pixels - 3)*0.5 + pixel_size

        x_grid = np.arange(grid_min, grid_max, pixel_size)
        y_grid = np.arange(grid_min, grid_max, pixel_size) 

        Ixy = np.zeros((n_pixels, n_pixels))

        gauss_x = np.zeros(shape=(len(self.x_pos), n_pixels))
        gauss_y = np.zeros(shape=(len(self.x_pos), n_pixels))
     

        for atom in range(len(self.x_pos)): # for each atom in the array 
        
            gauss_x = np.zeros((n_pixels)) #initializing array for gaussian
            gauss_y = np.zeros((n_pixels))
            
            #Selected neighbors
            x_neigh = np.where(np.absolute(x_grid - self.system_atoms[0,:][atom]) <= cutoff*sigma) # select neighbor pixels of a given atom when their distance <= cutoff*sigma
            y_neigh = np.where(np.absolute(y_grid - self.system_atoms[1,:][atom]) <= cutoff*sigma) #        
                
            gauss_x[x_neigh] = np.exp( -0.5 * ( ((x_grid[x_neigh] - self.system_atoms[0,:][atom]) /sigma)**2) ) # generates gaussian @all neighboring atoms depending on their distance from the given coordinate
            gauss_y[y_neigh] = np.exp( -0.5 * ( ((y_grid[y_neigh] - self.system_atoms[1,:][atom]) /sigma)**2) )
            
            Ixy = np.add(Ixy, gauss_x[:, None]*gauss_y)  # recursively adds all calculated gaussians 

        self.image = np.sqrt(2*np.pi) * sigma * Ixy

    def gen_image(self, n_pixels, pixel_size, sigma, cutoff):

        self.system_atoms = self.aligned_pos
        grid_min = -pixel_size * (n_pixels + 1)*0.5
        grid_max = pixel_size * (n_pixels - 3)*0.5 + pixel_size

        p = np.arange(len(self.x_pos))
        q = np.arange(grid_min, grid_max, pixel_size)
        x_grid, y_grid = np.meshgrid(q,p, indexing='ij')

        print(x_grid)
        print(np.shape(x_grid))
        Ixy = np.zeros((n_pixels, n_pixels))

        atom_gauss_x = np.zeros(shape=(len(self.x_pos), n_pixels))
        atom_gauss_y = np.zeros(shape=(len(self.x_pos), n_pixels))
     
        x_neigh = np.where(np.absolute(x_grid - self.system_atoms[0,:]) <= cutoff*sigma) # select neighbor pixels of a given atom when their distance <= cutoff*sigma
        y_neigh = np.where(np.absolute(y_grid - self.system_atoms[1,:]) <= cutoff*sigma) #

        atom_gauss_x[x_neigh] = np.exp( -0.5 * ( ((x_grid[x_neigh] - self.system_atoms[0,:]) /sigma)**2) ) # generates gaussian @all neighboring atoms depending on their distance from the given coordinate
        atom_gauss_y[y_neigh] = np.exp( -0.5 * ( ((y_grid[y_neigh] - self.system_atoms[1,:]) /sigma)**2) )

        Ixy =  atom_gauss_x[:, None]*atom_gauss_y  # recursively adds all calculated gaussians 

        self.image = np.sqrt(2*np.pi) * sigma * Ixy


    def plot_image(self):

        fig, ax = plt.subplots()
        ax.imshow(self.image)
        ax.set_title(str(self.num_protein) + ' ' + 'Proteins')
        plt.show()
        plt.savefig('init.png')
