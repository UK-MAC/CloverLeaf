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

!>  @brief Fortran acceleration kernel
!>  @author Wayne Gaudin
!>  @details The pressure and viscosity gradients are used to update the 
!>  velocity field.

MODULE accelerate_kernel_module

CONTAINS

SUBROUTINE accelerate_kernel(x_min,x_max,y_min,y_max,dt,     &
                             xarea,yarea,                    &
                             volume,                         &
                             density0,                       &
                             pressure,                       &
                             viscosity,                      &
                             xvel0,                          &
                             yvel0,                          &
                             xvel1,                          &
                             yvel1,                          &
                             stepbymass                      )

  IMPLICIT NONE

  INTEGER               :: x_min,x_max,y_min,y_max
  REAL(KIND=8)          :: dt

  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: density0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+2) :: xarea
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+3) :: yarea
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: pressure
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: viscosity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: xvel0,yvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: xvel1,yvel1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: stepbymass

  INTEGER               :: j,k
  REAL(KIND=8)          :: nodal_mass

!$OMP PARALLEL

!$OMP DO PRIVATE(nodal_mass)
  DO k=y_min,y_max+1
    DO j=x_min,x_max+1

      nodal_mass=(density0(j-1,k-1)*volume(j-1,k-1)  &
                 +density0(j  ,k-1)*volume(j  ,k-1)  &
                 +density0(j  ,k  )*volume(j  ,k  )  &
                 +density0(j-1,k  )*volume(j-1,k  )) &
                 *0.25_8

      stepbymass(j,k)=0.5_8*dt/nodal_mass

    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min,y_max+1
    DO j=x_min,x_max+1

      xvel1(j,k)=xvel0(j,k)-stepbymass(j,k)*(xarea(j  ,k  )*(pressure(j  ,k  )-pressure(j-1,k  ))    &
                                            +xarea(j  ,k-1)*(pressure(j  ,k-1)-pressure(j-1,k-1)))
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min,y_max+1
    DO j=x_min,x_max+1

      yvel1(j,k)=yvel0(j,k)-stepbymass(j,k)*(yarea(j  ,k  )*(pressure(j  ,k  )-pressure(j  ,k-1))    &
                                            +yarea(j-1,k  )*(pressure(j-1,k  )-pressure(j-1,k-1)))

    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min,y_max+1
    DO j=x_min,x_max+1

      xvel1(j,k)=xvel1(j,k)-stepbymass(j,k)*(xarea(j  ,k  )*(viscosity(j  ,k  )-viscosity(j-1,k  )) &
                                            +xarea(j  ,k-1)*(viscosity(j  ,k-1)-viscosity(j-1,k-1)))

    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min,y_max+1
    DO j=x_min,x_max+1

      yvel1(j,k)=yvel1(j,k)-stepbymass(j,k)*(yarea(j  ,k  )*(viscosity(j  ,k  )-viscosity(j  ,k-1)) &
                                            +yarea(j-1,k  )*(viscosity(j-1,k  )-viscosity(j-1,k-1)))

    ENDDO
  ENDDO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE accelerate_kernel

END MODULE accelerate_kernel_module
