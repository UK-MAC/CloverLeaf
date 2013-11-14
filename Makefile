PACKAGE=CloverLeaf

help:
	@echo "  CloverLeaf, a Lagrangian-Eulerian hydrodynamics mini-application   "
	@echo "                                                                     "
	@echo "  CloverLeaf comes in the following principal flavours:              "
	@echo "                                                                     "
	@echo "      - Ref                                                          "
	@echo "           This is the reference MPI/OpenMP implementation.          "
	@echo "      - Serial                                                       "
	@echo "           This is the serial implementation.                        "
	@echo "      - MPI                                                          "
	@echo "           This is the parallel MPI implementation.                  "
	@echo "      - OpenMP                                                       "
	@echo "           This version is threaded using MPI directives.            "
	@echo "      - OpenACC_KERNELS                                              "
	@echo "           This version uses OpenACC directives to utilise GPU       "
	@echo "           hardware.                                                 "
	@echo "      - OpenACC_LOOPS                                                "
	@echo "           This version uses OpenACC directives to utilise GPU       "
	@echo "           hardware.                                                 "
	@echo "      - OpenACC_CRAY                                                 "
	@echo "           This version uses OpenACC directives to utilise GPU       "
	@echo "           hardware.                                                 "
	@echo "      - Offload                                                      "
	@echo "      - CUDA                                                         "
	@echo "      - OpenMP 4                                                     "
	@echo "                                                                     "
	@echo "  Build a particular version by typing:                              "
	@echo "      \`make <flavour>\`                                             "
	@echo "  where <flavour> is the name of the version you want to make, all   "
	@echo "  in lowercase.                                                      "
	@echo "                                                                     "
	@echo "  Please use the COMPILER environment variable to specify the        "
	@echo "  compiler you wish to use. Supported compilers are:                 "
	@echo "      - GNU                                                          "
	@echo "      - CRAY                                                         "
	@echo "      - INTEL                                                        "
	@echo "      - PATHSCALE                                                    "
	@echo "      - PGI                                                    "
	@echo "      - SUN                                                    "
	@echo "      - XLF                                                   "

all: ref openmp mpi openacc_cray openacc_kernels openacc_loops cuda openmp4 serial offload

ref:
	cd $(PACKAGE)_ref; make

serial:
	cd $(PACKAGE)_Serial; make

openmp:
	cd $(PACKAGE)_OpenMP; make

mpi:
	cd $(PACKAGE)_MPI; make

openmp4:
	cd $(PACKAGE)_OpenMP4; make

openacc_cray:
	cd $(PACKAGE)_OpenACC_CRAY; make

openacc_kernels:
	cd $(PACKAGE)_OpenACC_KERNELS; make

openacc_loops:
	cd $(PACKAGE)_OpenACC_LOOPS; make

cuda:
	cd $(PACKAGE)_CUDA; make

offload:
	cd $(PACKAGE)_Offload; make

clean:
	cd $(PACKAGE)_ref; make clean
	cd $(PACKAGE)_Serial; make clean
	cd $(PACKAGE)_OpenMP; make clean
	cd $(PACKAGE)_MPI; make clean
	cd $(PACKAGE)_OpenMP4; make clean
	cd $(PACKAGE)_OpenACC_CRAY; make clean
	cd $(PACKAGE)_OpenACC_KERNELS; make clean
	cd $(PACKAGE)_OpenACC_LOOPS; make clean
	cd $(PACKAGE)_CUDA; make clean
	cd $(PACKAGE)_Offload; make clean
