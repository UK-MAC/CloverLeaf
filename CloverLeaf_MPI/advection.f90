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

!>  @brief Top level advection driver
!>  @author Wayne Gaudin
!>  @details Controls the advection step and invokes required communications.

MODULE advection_module

CONTAINS

SUBROUTINE advection()

  USE clover_module
  USE advec_cell_driver_module
  USE advec_mom_driver_module
  USE update_halo_module

  IMPLICIT NONE

  INTEGER :: sweep_number,direction,c

  INTEGER :: xvel,yvel

  INTEGER :: fields(NUM_FIELDS)

  sweep_number=1
  IF(advect_x)      direction=g_xdir
  IF(.not.advect_x) direction=g_ydir
  xvel=g_xdir
  yvel=g_ydir

  fields=0
  fields(FIELD_ENERGY1)=1
  fields(FIELD_DENSITY1)=1
  fields(FIELD_VOL_FLUX_X)=1
  fields(FIELD_VOL_FLUX_Y)=1
  CALL update_halo(fields,2)

  DO c=1,number_of_chunks
    CALL advec_cell_driver(c,sweep_number,direction)
  ENDDO

  fields=0
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1
  fields(FIELD_MASS_FLUX_X)=1
  fields(FIELD_MASS_FLUX_y)=1
  CALL update_halo(fields,2)

  DO c=1,number_of_chunks
    CALL advec_mom_driver(c,xvel,direction,sweep_number) 
  ENDDO
  DO c=1,number_of_chunks
    CALL advec_mom_driver(c,yvel,direction,sweep_number) 
  ENDDO

  sweep_number=2
  IF(advect_x)      direction=g_ydir
  IF(.not.advect_x) direction=g_xdir

  DO c=1,number_of_chunks
    CALL advec_cell_driver(c,sweep_number,direction)
  ENDDO

  fields=0
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1
  fields(FIELD_MASS_FLUX_X)=1
  fields(FIELD_MASS_FLUX_y)=1
  CALL update_halo(fields,2)

  DO c=1,number_of_chunks
    CALL advec_mom_driver(c,xvel,direction,sweep_number) 
  ENDDO
  DO c=1,number_of_chunks
    CALL advec_mom_driver(c,yvel,direction,sweep_number) 
  ENDDO

END SUBROUTINE advection

END MODULE advection_module
