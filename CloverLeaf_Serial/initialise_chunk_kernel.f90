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

!>  @brief Fortran chunk initialisation kernel.
!>  @author Wayne Gaudin
!>  @details Calculates mesh geometry for the mesh chunk based on the mesh size.

MODULE initialise_chunk_kernel_module

CONTAINS

SUBROUTINE initialise_chunk_kernel(x_min,x_max,y_min,y_max,       &
                                   xmin,ymin,dx,dy,               &
                                   vertexx,                       &
                                   vertexdx,                      &
                                   vertexy,                       &
                                   vertexdy,                      &
                                   cellx,                         &
                                   celldx,                        &
                                   celly,                         &
                                   celldy,                        &
                                   volume,                        &
                                   xarea,                         &
                                   yarea                          )

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  REAL(KIND=8) :: xmin,ymin,dx,dy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexdx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexy
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexdy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: cellx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celly
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celldy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+2) :: xarea
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+3) :: yarea

  INTEGER      :: j,k

  DO j=x_min-2,x_max+3
     vertexx(j)=xmin+dx*float(j-x_min)
  ENDDO

  DO j=x_min-2,x_max+3
    vertexdx(j)=dx
  ENDDO

  DO k=y_min-2,y_max+3
     vertexy(k)=ymin+dy*float(k-y_min)
  ENDDO

  DO k=y_min-2,y_max+3
    vertexdy(k)=dy
  ENDDO

  DO j=x_min-2,x_max+2
     cellx(j)=0.5*(vertexx(j)+vertexx(j+1))
  ENDDO

  DO j=x_min-2,x_max+2
     celldx(j)=dx
  ENDDO

  DO k=y_min-2,y_max+2
     celly(k)=0.5*(vertexy(k)+vertexy(k+1))
  ENDDO

  DO k=y_min-2,y_max+2
     celldy(k)=dy
  ENDDO

  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
        volume(j,k)=dx*dy
     ENDDO
  ENDDO

  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
        xarea(j,k)=celldy(k)
     ENDDO
  ENDDO

  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
        yarea(j,k)=celldx(j)
     ENDDO
  ENDDO

END SUBROUTINE initialise_chunk_kernel

END MODULE initialise_chunk_kernel_module
