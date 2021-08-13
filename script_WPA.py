from WPA import img_generator
from script_atom_coord import atom_pos, x_pos, num_protein, x_center, y_center, n_pixels, pixel_size

x_center = x_center
y_center = y_center
 
sigma = 1
cutoff = 15

#Create the generator
gen = img_generator(atom_pos, x_pos, num_protein)

#align coordinates
gen.align_system()

#Generate the image
gen.generate_image(n_pixels, pixel_size, sigma, cutoff)

# gen.gen_image(n_pixels, pixel_size, sigma, cutoff)

#Plot the image
gen.plot_image()


