#Crown Copyright 2012 AWE.
#
# This file is part of CloverLeaf.
#
# CloverLeaf is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or (at your option) 
# any later version.
#
# CloverLeaf is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details.
#
# You should have received a copy of the GNU General Public License along with 
# CloverLeaf. If not, see http://www.gnu.org/licenses/.

#  @brief Makefile for CloverLeaf
#  @author Wayne Gaudin, Andy Herdman
#  @details Agnostic, platform independent makefile for the Clover Leaf benchmark code.

# It is not meant to be clever in anyway, just a simple build out of the box script.
# Just make sure mpif90 is in your path. It uses mpif90 even for all builds because this abstracts the base
#  name of the compiler. If you are on a system that doesn't use mpif90, just replace mpif90 with the compiler name
#  of choice. The only mpi dependencies in this non-MPI version are mpi_wtime in timer.f90.

# There is no single way of turning OpenMP compilation on with all compilers.
# The known compilers have been added as a variable. By default the make
#  will use no options, which will work on Cray for example, but not on other
#  compilers.


FLAGS       = -O3 -Mcuda=cc35,6.0,maxregcount:64
ifdef DEBUG
  FLAGS       = -O0 -g -C -Mchkstk -Ktrap=fp -Mchkfpstk -Mchkptr
endif

MPI_COMPILER=mpif90 -I/opt/pgi/osx86-64/2014/mpi/mpich/include

clover_leaf: *.f90 *.cuf Makefile
	$(MPI_COMPILER) $(FLAGS)	\
	data.f90			\
	definitions.cuf			\
	pack_kernel.f90			\
	clover.f90			\
	report.f90			\
	timer.f90			\
	parse.f90			\
	read_input.f90			\
	initialise_chunk_kernel.f90	\
	initialise_chunk.f90		\
	build_field.f90			\
	update_halo_kernel.cuf		\
	update_halo.f90			\
	ideal_gas_kernel.cuf		\
	ideal_gas.cuf			\
	start.f90			\
	generate_chunk_kernel.f90	\
	generate_chunk.f90		\
	initialise.f90			\
	field_summary_kernel.cuf	\
	field_summary.f90		\
	viscosity_kernel.cuf		\
	viscosity.cuf			\
	calc_dt_kernel.cuf		\
	calc_dt.cuf			\
	timestep.f90			\
	accelerate_kernel.cuf		\
	accelerate.cuf			\
	revert_kernel.cuf		\
	revert.f90			\
	PdV_kernel.cuf			\
	PdV.cuf				\
	flux_calc_kernel.cuf		\
	flux_calc.f90			\
	advec_cell_kernel.cuf		\
	advec_cell_driver.cuf		\
	advec_mom_kernel.cuf		\
	advec_mom_driver.cuf		\
	advection.f90			\
	reset_field_kernel.cuf		\
	reset_field.f90			\
	hydro.cuf			\
	visit.f90			\
	clover_leaf.f90			\
	timer_c.c                       \
	-o clover_leaf; echo $(MESSAGE)

clean:
	rm -f a.out *.o *~ *.mod *genmod* *hmd* *.oo *.hmf *.lst *.cub *.ptx clover_leaf
