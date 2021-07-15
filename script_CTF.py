from script_wpa import n_pixels
from script_wpa import pixel_size

keV = 100   # accelerating voltage of microscope in keV
C_s = 2 * 10 ** 7 # spherical abberation constant in angstrom 
z = 10000 # defocus in angstrom
d = 220 # box size in angstrom  ... 3 x molecualar radius

from script_wpa import gen

phi = gen.image 

from CTF import transfer

fun = transfer(keV, C_s, z, phi, n_pixels, pixel_size, d)

fun.Lens_effects()

fun.plot_imageCTF()

fun.noise_machine()

fun.intensity()

fun.intensity_noise()

fun.plot_image_out()

fun.plot_image_noise()