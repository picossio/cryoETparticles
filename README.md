# Cryo-EM Micrograph Simulator  

## About the project
This package of python code allows for the simulation of cryogenic electron microscope images. It is intended to create noisy, cell-like enviornments which are filled with proteins. The output of this code is formatted to be used in conjuction with BioEM software. BioEM will perform a cross-correlation calculation to itendify protein locations in the generated micrograph.  
This project was made by Devon Kleeblatt in the summer of 2021 while he was working as an intern for the Flatiron Institute.

---
## How to run the scripts
1. Download the following files
   - atom_coord.py
   - script_atom_coord.py
   - WPA<l>.py
   - script_WPA.py
   - CTF<l>.py
   - script_CTF.py
   - Protein databank (pdb) file for the desired protein
2. Set up a python3 enviornment with the following imports
   - numpy
   - scipy
   - math
   - matplotlib
   - biopython
3. Edit the python scripts as descibed in the next section
#4. Run script_CTF.py
5. There are several outputs of this program, including:  
   - 2D projection of the proteins... visualization of the weak phase approximation
   - Images of the 1D, 2D contrast transfer function
   - Images of each micrograph
   - .txt file which contains:  
        1. (x, y) position of each protein in the micrographs  
        2. image intensity at each pixel in for each generated micrograph

## Notes on each file to create the micrograph
---
### **atom_coord.py**  
_No modification is necessary_

---
### **script_atom_coord.py**  
  
_line 14:_  
`struct = parser.get_structure('x','y.pdb')`  
x = protein name  
y = path to pdb file with protein template     

_lines 16-18:_ 
```
n_pixels = A  
pixel_size = B  
target = C
```
         
A = number of pixels in 1 dimension ... **A** x **A** micrograph  
B = pixel size in angstrom  
C = number of desired proteins in mcirograph

_lines 59-63:_  
```
if min < dis < max:
                    rand_val = rng.standard_normal(1)
                    if rand_val * dis < min:
                        use=0
                if dis<min:
```
**min** = minumum allowed distance between two proteins  
**max** = maximum distance value which has the potential to cause protein overlap
     
Protein distances between **min** and **max** will be multiplied by a random number too see if they clear the minimum distance.  

---
### **WPA<l>.py**
_No modification is necessary_

---
### **script_WPA.py**  
_lines 7-8:_
```
sigma = A 
cutoff = B
```
A = standard deviation of the gaussian spheres which represent atoms in the proteins  
B = limits the distance of how far gaussians will be added

---
### **CTF<l>.py**  
_line 99:_
```
        file_path = '/A/' + self.name + '.txt'

```
A = path to the file where the code's text output will save

___
### **script_CTF.py**  
_lines 3-9:_
```
keV = A 
C_s = B 
z = C 
d = n_pixels * pixel_size # box size in angstrom  ... 3 x molecualar radius
sigma = E 
sigma2 = F
name = str('G')
```
A = accelerating voltage of microscope in keV ~ 100-200  
B = spherical abberation constant ~ 2 * 10^7  
C = defocus in angstrom _... this can varry widely_ ~ 10,000  
d .... do not change  
E = standard deviation of gaussian noise  
F = standard deviation of gaussian noise model #2  

---
