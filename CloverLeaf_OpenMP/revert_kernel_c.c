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
 *  @brief C revert kernel.
 *  @author Wayne Gaudin
 *  @details Takes the half step field data used in the predictor and reverts
 *  it to the start of step data, ready for the corrector.
 *  Note that this does not seem necessary in this proxy-app but should be
 *  left in to remain relevant to the full method.
 */

#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void revert_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                      double *density0,
                      double *density1,
                      double *energy0,
                      double *energy1)

{
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;

  int j,k;
  
#pragma omp parallel
 {
#pragma omp for private(j)
  for (k=y_min;k<=y_max;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max;j++) {
      density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
    }
  }
  
#pragma omp for private(j)
  for (k=y_min;k<=y_max;k++) {
#pragma ivdep
    for (j=x_min;j<=x_max;j++) {
      energy1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=energy0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
    }
  }

 }

}
