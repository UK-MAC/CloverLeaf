!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief CloverLeaf top level program: Invokes the main cycle
!>  @author Wayne Gaudin
!>  @details CloverLeaf in a proxy-app that solves the compressible Euler
!>  Equations using an explicit finite volume method on a Cartesian grid.
!>  The grid is staggered with internal energy, density and pressure at cell
!>  centres and velocities on cell vertices.
!>
!>  A second order predictor-corrector method is used to advance the solution
!>  in time during the Lagrangian phase. A second order advective remap is then
!>  carried out to return the mesh to an orthogonal state.
!>
!>
!>  It can use OpenMP, OpenACC on a compute device.
!>
!>  NOTE: that the proxy-app uses uniformly spaced mesh. The actual method will
!>  work on a mesh with varying spacing to keep it relevant to it's parent code.
!>  For this reason, optimisations should only be carried out on the software
!>  that do not change the underlying numerical method. For example, the
!>  volume, though constant for all cells, should remain array and not be
!>  converted to a scalar.
PROGRAM clover_leaf

  USE clover_module

  IMPLICIT NONE


  CALL clover_init_comms()

  IF(parallel%boss)THEN
      WRITE(*,*)
      WRITE(*,'(a15,f8.3)') 'Clover Version ',g_version
      WRITE(*,*)
      WRITE(0,*)
      WRITE(0,'(a15,f8.3)') 'Clover Version ',g_version
      WRITE(0,*)
  ENDIF

  CALL initialise

  CALL hydro
  
  ! Deallocate everything
  
END PROGRAM clover_leaf

