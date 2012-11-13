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
 *  @brief C PdV kernel.
 *  @author Wayne Gaudin
 *  @details Calculates the change in energy and density in a cell using the
 *  change on cell volume due to the velocity gradients in a cell. The time
 *  level of the velocity data depends on whether it is invoked as the
 *  predictor or corrector.
 */

#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void pdv_kernel_c_(int *prdct,
                   int *xmin,int *xmax,int *ymin,int *ymax,
                double *dtbyt,
                double *xarea,
                double *yarea,
                double *volume,
                double *density0,
                double *density1,
                double *energy0,
                double *energy1,
                double *pressure,
                double *viscosity,
                double *xvel0,
                double *xvel1,
                double *yvel0,
                double *yvel1,
                double *volume_change)
{
  int predict=*prdct;
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;
  double dt=*dtbyt;

  int j,k;
  double recip_volume,energy_change,min_cell_volume,right_flux,left_flux,top_flux,bottom_flux,total_flux;
  
#pragma omp parallel
 {

  if(predict==0) {
    
#pragma omp for private(right_flux,left_flux,top_flux,bottom_flux,total_flux,min_cell_volume,energy_change,recip_volume,j)
    for (k=y_min;k<=y_max;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max;j++) {
        left_flux=  (xarea[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])
                                   *(xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt*0.5;
        right_flux= (xarea[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)])
                                   *(xvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt*0.5;
        bottom_flux=(yarea[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)])
                                   *(yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt*0.5;
        top_flux=   (yarea[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)])
                                   *(yvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt*0.5;

        total_flux=right_flux-left_flux+top_flux-bottom_flux;

        volume_change[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                         /(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+total_flux);

        min_cell_volume=MIN(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+right_flux-left_flux+top_flux-bottom_flux
                           ,MIN(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+right_flux-left_flux
                           ,volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+top_flux-bottom_flux));

        recip_volume=1.0/volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];

        energy_change=(pressure[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]/density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                     +viscosity[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]/density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)])
                      *total_flux*recip_volume;

        energy1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=energy0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]-energy_change;

        density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                           *volume_change[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];

      }
    }
  }
  else{
#pragma omp for private(right_flux,left_flux,top_flux,bottom_flux,total_flux,min_cell_volume,energy_change,recip_volume,j)
    for (k=y_min;k<=y_max;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max;j++) {
        left_flux=  (xarea[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])
                                   *(xvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +xvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel1[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt;
        right_flux= (xarea[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)])
                                   *(xvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)]
                                    +xvel1[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +xvel1[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt;
        bottom_flux=(yarea[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)])
                                   *(yvel0[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                    +yvel1[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt;
        top_flux=   (yarea[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)])
                                   *(yvel0[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel0[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel1[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]
                                    +yvel1[FTNREF2D(j+1,k+1,x_max+5,x_min-2,y_min-2)])
                                    *0.25*dt;

        total_flux=right_flux-left_flux+top_flux-bottom_flux;

        volume_change[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                         /(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+total_flux);

        min_cell_volume=MIN(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+right_flux-left_flux+top_flux-bottom_flux
                           ,MIN(volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+right_flux-left_flux
                           ,volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]+top_flux-bottom_flux));

        recip_volume=1.0/volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];

        energy_change=(pressure[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]/density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                     +viscosity[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]/density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)])
                      *total_flux*recip_volume;

        energy1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=energy0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]-energy_change;

        density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=density0[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                           *volume_change[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];

      }
    }

  }

 }

}

