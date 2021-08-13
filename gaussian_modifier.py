import numpy as np
import matplotlib.pyplot as plt

num_protein = 50
name = 'a_' + str(num_protein) + '_proteins'
n_pixels = 1000
sigma = 100

wave_intensity = np.loadtxt('/Users/dkleebatt/VS_Code_projects/Micrograph_Generator/txt_data_files/' + name + '.txt', skiprows=num_protein+5, usecols=(0,1,4))

rng = np.random.default_rng()
g_noise = rng.normal(0, sigma,  size=n_pixels ** 2) 

g_noise_value = wave_intensity[:,2] + g_noise

image_g_noise = np.stack((wave_intensity[:,0], wave_intensity[:,1], g_noise_value), axis=1)

print(image_g_noise)

file_path = '/Users/dkleebatt/VS_Code_projects/Micrograph_Generator/txt_data_files/' + 'g'+ str(sigma)+ '_' + name + '.txt'

for i in range(n_pixels ** 2):
    with open(file_path, mode='a') as f:

        f.write(str(int(image_g_noise[i,0])) + ' ' + str(int(image_g_noise[i,1])) + ' ' + str(image_g_noise[i,2]))
        f.write('\n')

std_noise = sigma
std_image = np.std(wave_intensity[:,2])
SNR = std_image / std_noise
print(SNR)

fig,ax = plt.subplots()

grid = g_noise_value.reshape(1000,1000)

ax.imshow(grid, cmap='gray')

ax.set_title(name + '_' + str(sigma) + '_Gaussian Noise')
plt.savefig(name + '_' + str(sigma) + '_g_noise.png')
plt.show()