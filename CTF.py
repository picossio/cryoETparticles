# plan
# take phi
# Fourier(phi)*CTF
# FT(phi)*CTF + colored noise 
# iFT(FT(phi)*CTF + colored noise)
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc 

class transfer:
    def __init__(self, keV, C_s, z, phi, n_pixels, pixel_size, d):
        
        self.C_s = C_s  
        self.z = z  
        self.phi = phi
        self.n_pixels = n_pixels
        self.pixel_size = pixel_size
        self.D = d

        self.eV = keV * 1000
        self.lam = sc.h / np.sqrt(2 * sc.electron_mass * self.eV * sc.electron_volt) / np.sqrt(1 + self.eV * sc.electron_volt / (2 * sc.electron_mass * sc.speed_of_light ** 2)) * 10 ** 10   #wavelength in angstrom
        self.k = 1 / self.lam  # wavenumber in 1 / angstrom
        self.theta_0 = 0.16 * self.lam  # unitless number descibiing angular acceptance of aperature

        print('lambda = ' + str(self.lam))
        print('k = ' + str(self.k))
        print('theta_0= ' + str(self.theta_0))

    def Lens_effects(self):
        
        FT_phi = np.fft.fft2(self.phi)   # outputs N x N 2d fourier space with each dimension being a spatial frequency
                                            # spatial frequency is discretized in 1/pixel_size units

        xi = np.zeros(shape = (self.n_pixels, self.n_pixels))

        for x in range(self.n_pixels):
            for y in range(self.n_pixels):

                xi[x,y] = np.sqrt(((x - self.n_pixels / 2) / self.pixel_size) ** 2 + (((y - self.n_pixels / 2) / self.pixel_size)) ** 2)    # assigned absolute value of spatial frequency for each pixel in fourier space

        theta = xi / (self.k * self.D)

        chi = 1/2 * self.k * self.z * theta ** 2 + self.k * self.C_s * theta ** 4 / 8 
        
        w_2 = 0.10  # constant descibing percentage of amplitude contrast in the image
        w_1 = np.sqrt(1 - w_2 ** 2)
        B = 1 #np.exp(- theta ** 2 / self.theta_0 ** 2)     # envelope function which descibes the decay of CTF amplitude

        self.CTF = B * (w_1 * np.sin(chi) - w_2 * np.cos(chi))
        self.phi_CTF = FT_phi * self.CTF
        
    def noise_machine(self):

        rng = np.random.default_rng()
        w_noise = rng.normal(0, 200, size = (self.n_pixels, self.n_pixels)) 
        w_noise = np.reshape(w_noise, (self.n_pixels, self.n_pixels))
        FT_w_noise = np.fft.fft2(w_noise)   

        rad_freq = np.zeros(shape = (self.n_pixels, self.n_pixels))
        f = np.zeros(shape = (self.n_pixels, self.n_pixels))
    
        for x in range(self.n_pixels):
            for y in range(self.n_pixels):

                rad_freq [x,y] = np.sqrt(((x - self.n_pixels / 2) / self.pixel_size) ** 2 + (((y - self.n_pixels / 2) / self.pixel_size)) ** 2)
                f[x,y] = 1 / np.sqrt((1 + rad_freq[x,y] ** 2))
        
        c_noise = f * FT_w_noise
        self.phi_CTF_noise = self.phi_CTF + c_noise

        print(np.shape(self.phi_CTF_noise))
        print(self.phi_CTF_noise)

    def plot_imageCTF(self):

            fig, ax = plt.subplots()
            ax.imshow(self.CTF)
            plt.show()

    def intensity(self):

        incident_wave = np.fft.ifft2(self.phi_CTF)
        self.wave_intensity = abs(incident_wave) ** 2 

    def intensity_noise(self):

        incident_wave_noise = np.fft.ifft2(self.phi_CTF_noise)
        self.wave_intensity_noise = abs(incident_wave_noise) ** 2 

    def plot_image_out(self):

        fig, ax = plt.subplots()
        ax.imshow(self.wave_intensity)
        plt.show()

    def plot_image_noise(self):

        fig, ax = plt.subplots()
        ax.imshow(self.wave_intensity_noise)
        plt.show()