/*Crown Copyright 2012 AWE.
*
* This file is part of CloverLeaf.
*
* CloverLeaf is free software: you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the
* Free Software Foundation, either version 3 of the License, or (at your option)
* any later version.
*
* CloverLeaf is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public License along with
* CloverLeaf. If not, see http://www.gnu.org/licenses/. */

/**
 *  @brief Not yet called.
 *  @author Wayne Gaudin
 *  @details Still just a stub.
 */

#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void generate_chunk_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                              double *vertexx,
                              double *vertexy,
                              double *cellx,
                              double *celly,
                              double *density0,
                              double *energy0,
                              double *xvel0,
                              double *yvel0,
                              int *nmbr_f_stts,
                              double *state_density,
                              double *state_energy,
                              double *state_xvel,
                              double *state_yvel,
                              double *state_xmin,
                              double *state_xmax,
                              double *state_ymin,
                              double *state_ymax,
                              double *state_radius,
                              double *state_geometry,
                              g_rct,
                              g_crc)


{
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;
  int number_of_states=*nmbr_f_stts;
  int g_rect=g_rct;
  int g_circ=g_crc;
  double radius;
  int state;

  int j,k,jt,kt;

#pragma omp parallel private(j)
 {
  /* State 1 is always the background state */
#pragma omp for
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
      energy0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=state_energy[1];
    }
  }

#pragma omp for private(j)
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
      density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=state_density[1];
   }
  }

#pragma omp for private(j)
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
      xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=state_xvel[1];
   }
  }

#pragma omp for private(j)
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
      yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=state_yvel[1];
   }
  }

  DO state=2,number_of_states

! Could the velocity setting be thread unsafe?

#pragma omp for private(radius)
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
        IF(state_geometry(state).EQ.g_rect ) THEN
          IF(vertexx(j).GE.state_xmin(state).AND.vertexx(j).LT.state_xmax(state)) THEN
            IF(vertexy(k).GE.state_ymin(state).AND.vertexy(k).LT.state_ymax(state)) THEN
              energy0(j,k)=state_energy(state)
              density0(j,k)=state_density(state)
              DO kt=k,k+1
                DO jt=j,j+1
                  xvel0(jt,kt)=state_xvel(state)
                  yvel0(jt,kt)=state_yvel(state)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ELSEIF(state_geometry(state).EQ.g_circ ) THEN
          radius=SQRT(cellx(j)*cellx(j)+celly(k)*celly(k))
          IF(radius.LE.state_radius(state))THEN
            energy0(j,k)=state_energy(state)
            density0(j,k)=state_density(state)
            DO kt=k,k+1
              DO jt=j,j+1
                xvel0(jt,kt)=state_xvel(state)
                yvel0(jt,kt)=state_yvel(state)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO

  ENDDO
