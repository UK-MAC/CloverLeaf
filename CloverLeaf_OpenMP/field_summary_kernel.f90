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

!>  @brief Fortran field summary kernel
!>  @author Wayne Gaudin
!>  @details The total mass, internal energy, kinetic energy and volume weighted
!>  pressure for the chunk is calculated.

MODULE field_summary_kernel_module

CONTAINS

SUBROUTINE field_summary_kernel(x_min,x_max,y_min,y_max, &
                                volume,                  &
                                density0,                &
                                energy0,                 &
                                pressure,                &
                                xvel0,                   &
                                yvel0,                   &
                                vol,mass,ie,ke,press     )

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: pressure
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
  REAL(KIND=8) :: vol,mass,ie,ke,press

  INTEGER      :: j,k,jv,kv
  REAL(KIND=8) :: vsqrd,cell_vol,cell_mass

  vol=0.0
  mass=0.0
  ie=0.0
  ke=0.0
  press=0.0

!$OMP PARALLEL
!$OMP DO PRIVATE(vsqrd,cell_vol,cell_mass) REDUCTION(+ : vol,mass,press,ie,ke)
  DO k=y_min,y_max
    DO j=x_min,x_max
      vsqrd=0.0
      DO kv=k,k+1
        DO jv=j,j+1
          vsqrd=vsqrd+0.25*(xvel0(jv,kv)**2+yvel0(jv,kv)**2)
        ENDDO
      ENDDO
      cell_vol=volume(j,k)
      cell_mass=cell_vol*density0(j,k)
      vol=vol+cell_vol
      mass=mass+cell_mass
      ie=ie+cell_mass*energy0(j,k)
      ke=ke+cell_mass*0.5*vsqrd
      press=press+cell_vol*pressure(j,k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE field_summary_kernel

END MODULE field_summary_kernel_module
