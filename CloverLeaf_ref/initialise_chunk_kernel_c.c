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

void initialise_chunk_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                                double *minx,
                                double *miny,
                                double *dx,
                                double *dy,
                                double *vertexx,
                                double *vertexdx,
                                double *vertexy,
                                double *vertexdy,
                                double *cellx,
                                double *celldx,
                                double *celly,
                                double *celldy,
                                double *volume,
                                double *xarea,
                                double *yarea)
{
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;
  double min_x=*minx;
  double min_y=*miny;
  double d_x=*dx;
  double d_y=*dy;

  int j,k;

#pragma omp parallel
 {
#pragma omp for private(j)
#pragma ivdep
  for (j=x_min-2;j<=x_max+3;j++) {
    vertexx[FTNREF1D(j  ,x_max+4,x_min-2)]=min_x+d_x*(double)(j-x_min);
  }

#pragma omp for private(j)
#pragma ivdep
  for (j=x_min-2;j<=x_max+3;j++) {
    vertexdx[FTNREF1D(j  ,x_max+4,x_min-2)]=d_x;
  }

#pragma omp for private(j)
#pragma ivdep
  for (k=y_min-2;k<=y_max+3;k++) {
    vertexy(k)=min_y+d_y*(double)(k-y_min);
  }

#pragma omp for private(j)
#pragma ivdep
  for (k=y_min-2;k<=y_max+3;k++) {
    vertexdy(k)=d_y:
  }

#pragma omp for private(j)
#pragma ivdep
  for (j=x_min-2;j<=x_max+2;j++) {
    cellx(j)=0.5*(vertexx(j)+vertexx(j+1));
  }

#pragma omp for private(j)
#pragma ivdep
  for (j=x_min-2;j<=x_max+2;j++) {
    celldx(j)=d_x;
  }

#pragma omp for private(j)
#pragma ivdep
  for (k=y_min-2;k<=y_max+2;k++) {
    celly(k)=0.5*(vertexy(k)+vertexy(k+1));
  }

#pragma omp for private(j)
#pragma ivdep
  for (k=y_min-2;k<=y_max+2;k++) {
     celldy(k)=d_y;
  }

#pragma omp for private(j)
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
        volume(j,k)=d_x*d_y;
    }
  }

#pragma omp for private(j)
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
        xarea(j,k)=celldy(k);
    }
  }

#pragma omp for private(j)
  for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
    for (j=x_min-2;j<=x_max+2;j++) {
        yarea(j,k)=celldx(j);
    }
  }

 }
