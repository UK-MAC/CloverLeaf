This directory contain a CUDA Fortran port of the serial Cloverleaf benchmark for a single 
GPU.

Managed Memory
--------------

This code utilizes the Managed Memory feature introduced in the CUDA 6.0 Toolkit
and in the PGI 14.7 compilers.  Managed Memory is a pool of memory accessible 
to both CPU and GPU using a single variable, which eliminates the need for separate host 
and device versions of data along with coding explicit data movement between host and 
device.  In CUDA Fortran, allocation into this pool of memory is achieved using the 
"managed" variable attribute in declarations (see definitions.cuf).  Managed Memory 
requires use of the 6.0 CUDA libraries and is available on devices with compute capability 
of 3.0 or higher.

One caveat in using Managed Memory is that on multi-GPU systems where any two GPUs are NOT 
peer-to-peer capable, the system falls back to using zero-copy memory which will typically
result in a large performance degradation.  Since this code utilizes a single GPU, the 
zero-copy fallback can be avoided by setting the environment variable CUDA_VISIBLE_DEVICES
to the appropriate GPU, or by setting the environment variable CUDA_MANAGED_FORCE_DEVICE_ALLOC
to a non-zero value.  (To determine the whether all pairs of GPUs on your system are 
peer-to-peer capable, compile and run the example code p2pAccess found in: 

/opt/pgi/*/2014/examples/CUDA-Fortran/CUDA-Fortran-Book/chapter4/P2P/p2pAccess/)


Texture Cache
-------------

On devices of compute capability 3.5 and higher, kernel arguments declared as 
"intent(in)" will be routed through the read-only texture cache via the LDG 
instruction.  As a result of this feature, there are neither use of explicit 
textures nor shared memory in any kernels.  On devices of compute capability 3.0, 
the code will run slower as the texture cache is not used.


Test Cases
----------

Four different input files are included with the code are those recommended to test for 
correctness.  Cloverleaf expects the input file to be named clover.in, so copy one of 
these files to clover.in before running a test.

The kinetic energy, as displayed in the output file clover.out every 10 timesteps and 
at the end of the run, is usually the most sensitive state variable and can be used as 
a correctness test.  The reference values for the kinetic energy for the last time step 
of these runs, along with approximate wall-clock times for a system with a K20 are listed 
below:

Input	       	     	  Kinetic Energy       Wall-clock time (s)
clover_bm_short.in	  0.1193E+01	       1.26 
clover_bm.in		  0.2590E+01	       43
clover_bm16_short.in	  0.3075E+00	       19
clover_bm16.in		  0.4854E+01	       650	  


