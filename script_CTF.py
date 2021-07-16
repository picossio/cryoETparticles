from script_wpa import n_pixels
from script_wpa import pixel_size

keV = 200   # accelerating voltage of microscope in keV
C_s = 2 * 10 ** 7 # spherical abberation constant in angstrom 
z = 10000 # defocus in angstrom
d = n_pixels * pixel_size # box size in angstrom  ... 3 x molecualar radius
sigma = .75 # standard deviation of noise 

from script_wpa import gen

# phi = gen.image 

import numpy as np 
import matplotlib.pyplot as plt

phi = np.zeros(shape = (n_pixels, n_pixels))

for x in range(n_pixels):

    for y in range(n_pixels):

        if (x - n_pixels/2) ** 2 + (y - n_pixels/2) ** 2 < 50 ** 2:
            
            phi[x,y] = 1 

fig, ax = plt.subplots()
ax.imshow(phi)
plt.show()

from CTF import transfer

fun = transfer(keV, C_s, z, phi, n_pixels, pixel_size, d, sigma)

fun.Lens_effects()

fun.plot_imageCTF()

fun.noise_machine()

fun.intensity()

fun.intensity_noise()

fun.plot_image_out()

fun.plot_image_noise()