# CloverLeaf Compilation and Usage ![image](Clover3D_alpha_small.png "CloverLeaf")

## CloverLeaf Directory Structure

Below the top level Leaf directory there is a directory called CloverLeaf. The 
sub directories in this directory contain the various implementation of the 
code.

* Serial - contains a serial version with no MPI or OpenMP
* OpenMP - contains an OpenMP version only with no MPI
* MPI - contains an MPI only implementation
* OpenACC - contains an OpenACC/MPI implementation that works under the Cray compiler
* HMPP- contains another OpenACC/MPI implementation that works with the CAPS and Cray compiler
* Offload - contains an Intel Offload/MPI implementation
* CUDA - contains the CUDA/MPI implementation
* Ref - contains a hybrid OpenMP/MPI implemention. The Serial, OpenMP and MPI 
versions are extracted from this version so should not diverge from it apart 
from the removal of the relevant software models.

## CloverLeaf Build Procedure

Dependencies in CloverLeaf have been kept to a minimum to ease the build 
process across multiple platforms and environments. There is a single makefile 
for all compilers and adding a new compiler should be straight forward.

In many case just typing `make` in the required software directory will work. 
This is the case if the mpif90 and mpicc wrappers are available on the system. 
This is true even for the Serial and OpenMP versions.

If the MPI compilers have different names then the build process needs to 
notified of this by defining two environment variables, `MPI_COMPILER` and 
`C_MPI_COMPILER`. 

For example on some Intel systems:

`make MPI_COMPILER=mpiifort C_MPI_COMPILER=mpiicc`

Or on Cray systems:

`make MPI_COMPILER=ftn C_MPI_COMPILER=cc`

### OpenMP Build

All compilers use different arguments to invoke OpenMP compilation. A simple 
call to make will invoke the compiler with -O3. This does not usually include 
OpenMP by default. To build for OpenMP for a specific compiler a further 
variable must be defined, `COMPILER` that will then select the correct option 
for OpenMP compilation. 

For example with the Intel compiler:

`make COMPILER=INTEL`

Which then append the -openmp to the build flags.

Other supported compiler that will be recognised are:-

* CRAY
* SUN
* GNU
* XL
* PATHSCALE
* PGI

The default flags for each of these is show below:-

* INTEL: -O3 -ipo
* SUN: -fast
* GNU: -ipo
* XL: -O5
* PATHSCLE: -O3
* PGI: -O3 -Minline
* CRAY: -em  _Note: that by default the Cray compiler with pick the optimum 
options for performance._

### Other Flags

The default compilation with the COMPILER flag set chooses the optimal 
performing set of flags for the specified compiler, but with no hardware 
specific options or IEEE compatability.

To produce a version that has IEEE compatiblity a further flag has to be set on 
the compiler line.

`make COMPILER=INTEL IEEE=1`

This flag has no effect if the compiler flag is not set because IEEE options 
are always compiler specific.

For each compiler the flags associated with IEEE are shown below:-

* INTEL: -fp-model strict –fp-model source –prec-div –prec-sqrt
* CRAY: -hpflex_mp=intolerant
* SUN: -fsimple=0 –fns=no
* GNU: -ffloat-store
* PGI: -Kieee
* PATHSCALE: -mieee-fp
* XL: -qstrict –qfloat=nomaf

Note that the MPI communications have been written to ensure bitwise identical 
answers independent of core count. However under some compilers this is not 
true unless the IEEE flags is set to be true. This is certainly true of the 
Intel and Cray compiler. Even with the IEEE options set, this is not guarantee 
that different compilers or platforms will produce the same answers. Indeed a 
Fortran run can give different answers from a C run with the same compiler, 
same options and same hardware.

Extra options can be added without modifying the makefile by adding two further 
flags, `OPTIONS` and `C_OPTIONS`, one for the Fortran and one for the C options.

`make COMPILER=INTEL OPTIONS=-xavx C_OPTIONS=-xavx`

A build for a Xeon Phi would just need the -xavx option above replaced by -mmic.

Finally, a `DEBUG` flag can be set to use debug options for a specific compiler.

`make COMPILER=PGI DEBUG=1`

These flags are also compiler specific, and so will depend on the `COMPILER` 
environment variable.

So on a system without the standard MPI wrappers, for a build that requires 
OpenMP, IEEE and AVX this would look like so:-

```
make COMPILER=INTEL MPI_COMPILER=mpiifort C_MPI_COMPILER=mpiicc IEEE=1 \
OPTIONS="-xavx" C_OPTIONS="-xavx"
```

### Vectorisation

Fortran tends to vectorise without the use the specific pragmas due to the 
higher level definition of data compared to C. C almost always needs pragmas to 
ensure that the compiler knows loops are safe to vectorise. Unfortunately there 
is no common standard for vector pragmas across compiler vendors, though 
\#pragma ivdep works on many. This means that on some systems (e.g. IBM) the 
vector pragmas may need to be modified to attain peak performance of the C code. 
Care needs to be taken with forcing vectorisation because even loops without 
obviously data dependencies can calculate the wrong answers.

### OpenACC Build

The makefile for this build does not differ from the Base version. In this case 
it is important to have the correct environment loaded on the system of use. 
On the Cray systems this will usually involve loading some NVIDIA and 
accelerator modules. Without these the code will still compile and run but the 
OpenACC pragmas wil be ignored and the calculation will take place on the CPU.

### HMPP Build

The makefile for this build does not differ from the Base version. 
It is also important to have the correct environment loaded on the system of use. 
This will usually involve loading some NVIDIA and accelerator modules. Without 
these the code will still compile and run but the OpenACC pragmas wil be ignored
and the calculation will take place on the CPU. To use the CAPS HMPP compiler for
the OpenACC pragmas it is neccessary to modify the make line so that the MPI_COMPILER
variable is pree-fixed with "hmpp" plus the addition of any system specific flags, say
for example nvcc flags is the code is targetting a GPU.


### Co-Array Build

A Co-Array Fortran build of CloverLeaf will be included in a future release, 
and is currently under production and testing.

### Shmem Build

A Shmem build of CloverLeaf will be included in a future release, and is 
currently under production and testing.

### OpenCL Build

An OpenCL build of CloverLeaf will be included in a future release, and is 
currently under production and testing.

### CUDA Build

A CUDA build of CloverLeaf is now included in the repository, and has 
been tested on large GPU systems successfully. One extra command is required
to compile to take account on the NVIDIA GPU being used. Set NV_ARCH=FERMI for
M2090 cards and NV_ARCH=KEPLER for K20 cards.

## Running the Code

CloverLeaf takes no command line arguments. It expects to find a file called 
`clover.in` in the directory it is running in.

There are a number of input files that come with the code. To use any of these 
they simply need to be copied to `clover.in` in the run directory and 
CloverLeaf invoked. The invocation is system dependent. 

For example for a hybrid run:

```
export OMP_NUM_THREADS=4

mpirun -np 8 clover_leaf
```

### Weak and Strong Scaling

Note that with strong scaling, as the task count increases for the same size 
global problem, the memory use of each task decreases. Eventually, the mesh 
data starts to fit into the various levels of cache. So even though the 
communications overhead is increasing, super-scalar leaps in performance can be 
seen as task count increases. Eventually all cache benefits are gained and the 
communications dominate. 

For weak scaling, memory use stays close to constant and these super-scalar 
increases aren't seen but the communications overhead stays constant relative 
to the computational overhead, and scaling remains good.

### Other Issues to Consider

System libraries and settings can also have a significant effect on performance. 

The use of the `HugePage` library can make memory access more efficient. The 
implementation of this is very system specific and the details will not be 
expanded here. 

Variation in clock speed, such as `SpeedStep` or `Turbo Boost`, is also 
available on some hardware and care needs to be taken that the settings are 
known. 

Many systems also allow some level of hyperthreading at a core level. These 
usually share floating point units and for a code like CloverLeaf, which is 
floating point intensive and light on integer operations, are unlikely to 
produce a benefit and more likely to reduce performance.

CloverLeaf is considered a memory bound code. Most data does not stay in cache 
very long before it is replaced. For this reason, memory speed can have a 
significant effect on performance. For the same reason, the same is true of 
hardware caches.

## Testing the Results

Even though bitwise answers cannot be expected across systems, answers should be 
very close. A summary print of state variables is printed out by default every 
ten hydrodynamic steps and then at the end of the run. This print gives average 
value of the volume, mass, density, pressure, kinetic energy and internal 
energy. 

As all boundaries are reflective the volume, mass, density and total energy 
should remain constant, though the scheme is not exactly conservative in energy 
due to the nature of the staggered grid. Kinetic energy is usually the most 
sensitive state variable and this is most likely to show compiler/system 
differences. If mass and volume do not stay constant through a run, then 
something is seriously wrong.

There is a very simple, small self test include in the CloverLeaf. If the code is
invoked no clover.in input file present, this test will be run and the answer tested
against a "known" solution.

There are four standard input files that are recommended for testing. 
Initially it is suggested than `clover_bm_short.in` is run. This is not a very 
sensitive test and the kinetic energy at the end of this run should be 0.1193E+01.
It is quick to run, even on a single core, and should stop after 87 steps.

The second test to try is `clover_bm.in`. This runs for 2955 timesteps and is 
more sensitive than the first test. Through this simulation the whole 
computational mesh in traversed by a shock and so it is a good test of the 
parallel implementation because all internal boundaries will be crossed during 
the course of the simulation. The final kinetic energy should be 0.2590E+01.

The third test to try is `clover_bm16_short.in`. This is the "socket" test and 
has a much larger mesh size and therefore, memory footprint. The final kinetic 
energy should be 0.3075E+00.

The last test to run for validation purposes is `clover_bm16.in`. This is a 
fairly long, large mesh run and the kinetic energy at the final time should be 
0.4854E+01.
