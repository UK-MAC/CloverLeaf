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

!>  @brief Fortran kernel to update the external halo cells in a chunk.
!>  @author Wayne Gaudin
!>  @details Updates halo cells for the required fields at the required depth
!>  for any halo cells that lie on an external boundary. The location and type
!>  of data governs how this is carried out. External boundaries are always
!>  reflective.

MODULE update_halo_kernel_module

CONTAINS

  SUBROUTINE update_halo_kernel(x_min,x_max,y_min,y_max,                            &
                        left,bottom,right,top,                                      &
                        left_boundary,bottom_boundary,right_boundary,top_boundary,  &
                        chunk_neighbours,                                           &
                        density0,                                                   &
                        energy0,                                                    &
                        pressure,                                                   &
                        viscosity,                                                  &
                        soundspeed,                                                 &
                        density1,                                                   &
                        energy1,                                                    &
                        xvel0,                                                      &
                        yvel0,                                                      &
                        xvel1,                                                      &
                        yvel1,                                                      &
                        vol_flux_x,                                                 &
                        vol_flux_y,                                                 &
                        mass_flux_x,                                                &
                        mass_flux_y,                                                &
                        fields,                                                     &
                        depth                                                       )
  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  INTEGER :: left,bottom,right,top
  INTEGER :: left_boundary,bottom_boundary,right_boundary,top_boundary
  INTEGER, DIMENSION(4) :: chunk_neighbours
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: pressure,viscosity,soundspeed
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1,energy1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1,yvel1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x,mass_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y,mass_flux_y
  INTEGER :: fields(:),depth

  ! These need to be kept consistent with the data module to avoid use statement
  INTEGER,      PARAMETER :: CHUNK_LEFT   =1    &
                            ,CHUNK_RIGHT  =2    &
                            ,CHUNK_BOTTOM =3    &
                            ,CHUNK_TOP    =4    &
                            ,EXTERNAL_FACE=-1

  INTEGER,      PARAMETER :: FIELD_DENSITY0   = 1         &
                            ,FIELD_DENSITY1   = 2         &
                            ,FIELD_ENERGY0    = 3         &
                            ,FIELD_ENERGY1    = 4         &
                            ,FIELD_PRESSURE   = 5         &
                            ,FIELD_VISCOSITY  = 6         &
                            ,FIELD_SOUNDSPEED = 7         &
                            ,FIELD_XVEL0      = 8         &
                            ,FIELD_XVEL1      = 9         &
                            ,FIELD_YVEL0      =10         &
                            ,FIELD_YVEL1      =11         &
                            ,FIELD_VOL_FLUX_X =12         &
                            ,FIELD_VOL_FLUX_Y =13         &
                            ,FIELD_MASS_FLUX_X=14         &
                            ,FIELD_MASS_FLUX_Y=15         &
                            ,NUM_FIELDS       =15

  INTEGER :: j,k

!$OMP PARALLEL

  ! Update values in external halo cells based on depth and fields requested
  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          density0(j,1-k)=density0(j,0+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          density0(j,y_max+k)=density0(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          density0(1-j,k)=density0(0+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          density0(x_max+j,k)=density0(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          density1(j,1-k)=density1(j,0+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          density1(j,y_max+k)=density1(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          density1(1-j,k)=density1(0+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          density1(x_max+j,k)=density1(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          energy0(j,1-k)=energy0(j,0+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          energy0(j,y_max+k)=energy0(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          energy0(1-j,k)=energy0(0+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          energy0(x_max+j,k)=energy0(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          energy1(j,1-k)=energy1(j,0+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          energy1(j,y_max+k)=energy1(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          energy1(1-j,k)=energy1(0+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          energy1(x_max+j,k)=energy1(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          pressure(j,1-k)=pressure(j,0+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          pressure(j,y_max+k)=pressure(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          pressure(1-j,k)=pressure(0+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          pressure(x_max+j,k)=pressure(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          viscosity(j,1-k)=viscosity(j,0+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          viscosity(j,y_max+k)=viscosity(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          viscosity(1-j,k)=viscosity(0+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          viscosity(x_max+j,k)=viscosity(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_XVEL0).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          xvel0(j,1-k)=xvel0(j,1+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          xvel0(j,y_max+1+k)=xvel0(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          xvel0(1-j,k)=-xvel0(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          xvel0(x_max+1+j,k)=-xvel0(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_XVEL1).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          xvel1(j,1-k)=xvel1(j,1+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          xvel1(j,y_max+1+k)=xvel1(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          xvel1(1-j,k)=-xvel1(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          xvel1(x_max+1+j,k)=-xvel1(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_YVEL0).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          yvel0(j,1-k)=-yvel0(j,1+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          yvel0(j,y_max+1+k)=-yvel0(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          yvel0(1-j,k)=yvel0(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          yvel0(x_max+1+j,k)=yvel0(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_YVEL1).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          yvel1(j,1-k)=-yvel1(j,1+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          yvel1(j,y_max+1+k)=-yvel1(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          yvel1(1-j,k)=yvel1(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          yvel1(x_max+1+j,k)=yvel1(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          vol_flux_x(j,1-k)=vol_flux_x(j,1+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          vol_flux_x(j,y_max+k)=vol_flux_x(j,y_max-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          vol_flux_x(1-j,k)=-vol_flux_x(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          vol_flux_x(x_max+j+1,k)=-vol_flux_x(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          mass_flux_x(j,1-k)=mass_flux_x(j,1+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+1+depth
          mass_flux_x(j,y_max+k)=mass_flux_x(j,y_max-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          mass_flux_x(1-j,k)=-mass_flux_x(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+depth
        DO j=1,depth
          mass_flux_x(x_max+j+1,k)=-mass_flux_x(x_max+1-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          vol_flux_y(j,1-k)=-vol_flux_y(j,1+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          vol_flux_y(j,y_max+k+1)=-vol_flux_y(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          vol_flux_y(1-j,k)=vol_flux_y(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          vol_flux_y(x_max+j,k)=vol_flux_y(x_max-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    IF(chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          mass_flux_y(j,1-k)=-mass_flux_y(j,1+k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
      DO k=1,depth
!$OMP DO
        DO j=x_min-depth,x_max+depth
          mass_flux_y(j,y_max+k+1)=-mass_flux_y(j,y_max+1-k)
        ENDDO
!$OMP END DO
      ENDDO
    ENDIF
    IF(chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          mass_flux_y(1-j,k)=mass_flux_y(1+j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF(chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-depth,y_max+1+depth
        DO j=1,depth
          mass_flux_y(x_max+j,k)=mass_flux_y(x_max-j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

!$OMP END PARALLEL

END SUBROUTINE update_halo_kernel

END  MODULE update_halo_kernel_module
