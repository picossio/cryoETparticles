from wpa import img_generator
from script_atom_coord import atom_pos
from script_atom_coord import x_pos

n_pixels = 256     # more pixel means the image will be able to include more information
pixel_size = 1.25    # increase pixel size to cover a larger area, more grainy image
sigma = 1
cutoff = 30

#Create the generator
gen = img_generator(atom_pos, x_pos)

#align coordinates
gen.align_system()

#Generate the image
gen.generate_image(n_pixels, pixel_size, sigma, cutoff)

#Plot the image
gen.plot_image()


