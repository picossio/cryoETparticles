import numpy as np
import matplotlib.pyplot as plt

cc = np.loadtxt('/Users/dkleebatt/VS_Code_projects/Micrograph_Generator/CC_a_50/CROSS_CORRELATION_50_sphere', skiprows=5, usecols=(2,3,4))

shape = np.shape(cc)
print(shape)

cc_map = np.zeros(shape=(100,100))

for i in range(int(shape[0])):
    x = int(cc[i,0]/10)
    y = int(cc[i,1]/10)
    cc_map[x,y] = cc[i,2]  

print(np.shape(cc_map))

# for x in range(100):
#     for y in range(100):
#         print(cc_map[x,y])

fig, ax = plt.subplots()
ax.set_title('Cross Correlation')
ax.imshow(cc_map)
plt.show()
