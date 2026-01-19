% Implementation of the Noise Clinic algorithm.

# ABOUT

* Authors    : Marc Lebrun <marc.lebrun.ik@gmail.com>, Miguel Colom <colom@cmla.ens-cachan.fr>
* Copyright : (C) 2015 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

# OVERVIEW

This source code provides an implementation of the Noise Clinic blind denoising algorithm.

# UNIX/LINUX/MAC USER GUIDE

The code is compilable on Unix/Linux and Mac OS. 

- Compilation. 
Automated compilation requires the make program.

- Library. 
This code requires the libpng and the libfftw libraries.

- Image format. 
Only the PNG format is supported. 
 
-------------------------------------------------------------------------
Usage:
1. Download the code package and extract it. Go to that directory. 

2. Compile the source code (on Unix/Linux/Mac OS). 
There are two ways to compile the code. 
(1) RECOMMENDED, with Open Multi-Processing multithread parallelization 
(http://openmp.org/). Roughly speaking, it accelerates the program using the 
multiple processors in the computer. Run
make OMP=1

OR
(2) If the compiler does not support OpenMp, run
make

3. Run NoiseClinic image denoising.
./NoiseClinic
You can decide on the number of scales you want to use with the nbScales parameter. Use 2 by default (value required).
You can decide to multiply the estimated noise curves at each scale by a coefficient factor. Use 1.0 by default (value required).
You can decide to print on screen some information. Use verbose = 1 in this case (value required).
You need to specify an input noisy image (ImNoisy in the following example), an output denoised image (ImDenoised) and
an output image (ImDiff) containing the differences between ImNoisy and ImDenoised.
ALL IMAGES MUST HAVE A SPECIFIED .png EXTENSION.
 
Example, run
./NoiseClinic ImNoisy.png ImDenoised.png ImDiff.png 2 1.0 1

# ABOUT THIS FILE

Copyright 2015 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
