# Beam Truss Model

A MATLAB code to generate the mass and stiffnes matrices for a beam truss structure and conduct a modal analysis. The beam model is based on the Euler-Bernoulli assumptions. 

**Features**
- generation of mass and stiffness matrices 
- **modal analysis**
- geometry plot
- **dynamic visualization of the modes in a video**
- additional inertia of connection elements is accounted for

<video src="https://user-images.githubusercontent.com/121583390/225946434-ba5fbfa8-a669-41bf-bdda-a836b786e69b.mp4" controls="controls" autoplay loop>
</video>


### How to use

**1. Setup**
<p align="center">
<img src="https://user-images.githubusercontent.com/121583390/225952612-79bbb595-7614-460f-b84e-7864a2766d91.jpg" width="400">
</p>

Specify the beam truss using the matrices `P` and `Beam`.

The rows of P are the connection points between the beams. The first column numbers the connecton points, the second until fourth column describe the x- y- and z-coordinate. 

for example

```
P= [1         0   -0.6351    0.2935, 
    2   -0.2750   -0.1588    0.2935, 
    ...
    10        0         0         0];
```
The rows of `Beam` each stand for a beam of the truss. The beams are numbered in the first column. The second and third column specify the number of the start and end connection point of the respective beam. The numbers of the connection points had been previously defined in `P`.

for example

```
Beam = [1,  8,  1
        2,  6,  8
        ...
        24, 7,  9];
```
Generate an instance of the `beam_truss` class. Specifiy the number `n` of elements per beam.   

`my_truss = beam_truss(P,Beam,n);`

The geometry is plotted subsequently.

**2. Modal Analsysis**

Conduct a modal analysis using the command 

`modal_analysis(my_truss)`

and visualize a selected mode `i_mode` 

`disp_modal_analysis(my truss, i_mode)`

Please note that by default only the first 35 modes are calculated, but that can be changed easily in the code so you can also look into higher-order modes. Since the structure is not fixed, the first six eigenmodes will be rigid body modes. 

### Troubleshooting

The code should be fairly well documented. If you find a bug, please open an issue and I will look into it. 

### Literature

[1] Géradin, Michel, and Daniel J. Rixen. Mechanical vibrations: theory and application to structural dynamics. John Wiley & Sons, 2014. _**(highly recommended)**_

[2] Shabana, Ahmed A. Dynamics of multibody systems. Cambridge university press, 2020.

[![View Beam Truss on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://de.mathworks.com/matlabcentral/fileexchange/126420-beam-truss)
