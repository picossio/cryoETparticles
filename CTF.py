import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc

class transfer:
    def __init__(self, keV, C_s, z, phi, n_pixels, pixel_size, d, sigma, sigma2, name, num_protein, x_center, y_center):
        
        self.C_s = C_s  
        self.z = z  
        self.phi = phi
        self.n_pixels = n_pixels
        self.pixel_size = pixel_size
        self.D = d
        self.sigma = sigma
        self.sigma2 = sigma2 
        self.name = name
        self.x_center = x_center
        self.y_center = y_center
        self.num_protein = num_protein

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

        # for x in range(self.n_pixels):
        #     for y in range(self.n_pixels):

                # xi[x,y] = np.sqrt(((x - self.n_pixels / 2) / self.pixel_size) ** 2 + (((y - self.n_pixels / 2) / self.pixel_size)) ** 2)    # assigned absolute value of spatial frequency for each pixel in fourier space
        for x in range(int(self.n_pixels/2)):
            for y in range(int(self.n_pixels/2)):

                xi[x,y] = np.sqrt(((x) / self.pixel_size) ** 2 + (((y) / self.pixel_size)) ** 2)
                xi[self.n_pixels -1 - x, self.n_pixels - 1 - y] = np.sqrt(((x) / self.pixel_size) ** 2 + (((y) / self.pixel_size)) ** 2)
                xi[self.n_pixels -1 - x, y] = np.sqrt(((x) / self.pixel_size) ** 2 + (((y) / self.pixel_size)) ** 2)
                xi[x, self.n_pixels - 1 - y] = np.sqrt(((x) / self.pixel_size) ** 2 + (((y) / self.pixel_size)) ** 2)
        
        theta = xi / (self.k * self.D)

        chi = 1/2 * self.k * self.z * theta ** 2 + 1/8 * self.k * self.C_s * theta ** 4
        
        w_2 = 0.10  # constant descibing percentage of amplitude contrast in the image
        w_1 = np.sqrt(1 - w_2 ** 2)
        B = 1 #np.exp(- theta ** 2 / self.theta_0 ** 2)     # envelope function which descibes the decay of CTF amplitude

        self.CTF = B * (w_1 * np.sin(chi) - w_2 * np.cos(chi))
        self.phi_CTF = self.CTF * FT_phi

    def plot_imageCTF(self):

            fig, ax = plt.subplots()
            ax.imshow(self.CTF, cmap='gray')
            # plt.savefig(self.name + '_CTF.png')
            plt.show()
            
    def intensity(self):

        incident_wave = np.fft.ifft2(self.phi_CTF)
        self.wave_intensity = np.real(incident_wave)
        self.wave_intensity_range = np.clip(self.wave_intensity, -100,  0.01)

    def gaussian_noise_machine(self):

        rng = np.random.default_rng()
        g_noise = rng.normal(0,self.sigma, size = (self.n_pixels, self.n_pixels)) 
        g_noise = np.reshape(g_noise, (self.n_pixels, self.n_pixels))
        self.image_g_noise = self.wave_intensity + g_noise

    def gaussian_noise_machine2(self):

        rng = np.random.default_rng()
        g_noise = rng.normal(0,self.sigma2, size = (self.n_pixels, self.n_pixels)) 
        g_noise = np.reshape(g_noise, (self.n_pixels, self.n_pixels))
        self.image_g_noise2 = self.wave_intensity + g_noise

    def poisson_noise_machine(self):
        rng = np.random.default_rng()
        self.image_p_noise = np.zeros(shape=(self.n_pixels,self.n_pixels))
        minimum = np.min(self.wave_intensity)
    
        for x in range(self.n_pixels):
            for y in range(self.n_pixels):
                        
                p_noise = rng.poisson(self.wave_intensity[x,y] - minimum)
                self.image_p_noise[x,y] = self.wave_intensity[x,y] + p_noise

    def output_data(self):

        file_path = '/Users/dkleebatt/VS_Code_projects/Micrograph_Generator/txt_data_files/' + self.name + '.txt'

        with open(file_path, mode='a') as f:
            
            f.write('#' + str(self.num_protein) + ' Proteins' + ' @ (x,y)' + '\n' + '\n')

            for i in range(int(self.num_protein)):
                f.write(str(self.x_center[i]) + ' ' + str(self.y_center[i]))
                f.write('\n')

            f.write('\n' + '# ' + 'x, ' + 'y, ' + 'incident wave, ' + 'CTF, ' + 'image, ' + 'image + gaussian, ' + 'image + gaussian2, ' + 'image + poisson' + '\n' + '\n')

        for x in range(int(self.n_pixels)):
                for y in range(int(self.n_pixels)):

                    with open(file_path, mode='a') as f:

                        f.write(str(x) + ' ' + str(y) + ' ' + str(self.phi[x,y]) + ' ' + str(self.CTF[x,y]) + ' ' + str(self.wave_intensity[x,y]) + ' ' + str(self.image_g_noise[x,y]) +  ' ' + str(self.image_g_noise2[x,y]) + ' ' + str(self.image_p_noise[x,y]))
                        f.write('\n')
                       
    def plot_image_out(self):

        fig, ax = plt.subplots()
        ax.imshow(self.wave_intensity, cmap='gray')
        ax.set_title(self.name + ' CTF Convoluted Image')
        # plt.savefig(self.name + '_no_noise.png')
        plt.show()

    def plot_image_g_noise(self):

        fig,ax = plt.subplots()
        ax.imshow(self.image_g_noise, cmap='gray')
        ax.set_title(self.name + ' Gaussian Noise' + 'sigma = ' + str(self.sigma))
        plt.savefig(self.name + '_g_noise.png')
        plt.show()

    def plot_image_g_noise2(self):

        fig,ax = plt.subplots()
        ax.imshow(self.image_g_noise2, cmap='gray')
        ax.set_title(self.name + ' Gaussian Noise' + 'sigma = ' + str(self.sigma2))
        # plt.savefig(self.name + '_g_noise2.png')
        plt.show()

    def plot_image_p_noise(self):

        fig,ax = plt.subplots()
        ax.imshow(self.image_p_noise, cmap='gray')
        ax.set_title(self.name + ' Poisson Noise')
        # plt.savefig(self.name + '_p_noise.png')
        plt.show()

    def CTF_1d(self):
        
        self.frequency = np.zeros(self.n_pixels)

        for i in range (int(self.n_pixels)):
        
            self.frequency[i] = i / self.pixel_size

        theta = self.frequency / (self.k * self.D)
        chi = 1/2 * self.k * self.z * theta ** 2 + 1/8 * self.k * self.C_s * theta ** 4
        w_2 = 0.10  
        w_1 = np.sqrt(1 - w_2 ** 2)
        B = np.exp(- theta ** 2 / self.theta_0 ** 2)    

        self.CTF_1d = B * (w_1 * np.sin(chi) - w_2 * np.cos(chi))

        print(np.shape(self.plot_CTF_1d))

    def plot_CTF_1d(self):

        fig,ax = plt.subplots()
        ax.plot(self.frequency, self.CTF_1d)
        ax.set_title(self.name + ' 1D CTF' + ' z = ' + str(self.z ))
        ax.set_xlabel('Frequency')
        ax.set_ylabel('CTF')
        # plt.savefig(self.name + '_1d_CTF.png')
        plt.show()

    # definitions below this comment are not used in the current implementation

    def noise_machine(self):

        rng = np.random.default_rng()
        w_noise = rng.normal(0, self.sigma, size = (self.n_pixels, self.n_pixels)) 
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
