from script_WPA import n_pixels, pixel_size, num_protein, x_center, y_center

keV = 100   # accelerating voltage of microscope in keV
C_s = 2 * 10 ** 7 # spherical abberation constant in angstrom 
z = 10000 # defocus in angstrom
d = n_pixels * pixel_size # box size in angstrom  ... 3 x molecualar radius
sigma = 2 # standard deviation of gaussian noise
sigma2 = 5
name = str('doog')

from script_WPA import gen
    
phi = gen.image 

from CTF import transfer

fun = transfer(keV, C_s, z, phi, n_pixels, pixel_size, d, sigma, sigma2, name, num_protein, x_center, y_center)

fun.Lens_effects()

fun.CTF_1d()

fun.plot_imageCTF()

fun.plot_CTF_1d()

fun.intensity()

fun.gaussian_noise_machine()

fun.gaussian_noise_machine2()

fun.poisson_noise_machine()

fun.plot_image_out()

fun.plot_image_g_noise()

fun.plot_image_p_noise()

fun.output_data()  


# plots a circle to check how the ctf is working 
# import numpy as np
# import matplotlib.pyplot as plt

# # phi = np.zeros(shape = (n_pixels, n_pixels))

# # for x in range(n_pixels):

# #     for y in range(n_pixels):

# #         if (x - n_pixels/2) ** 2 + (y - n_pixels/2) ** 2 < 50 ** 2:
            
# #             phi[x,y] = 1 

# # fig, ax = plt.subplots()
# # ax.imshow(phi)
# # plt.show()
