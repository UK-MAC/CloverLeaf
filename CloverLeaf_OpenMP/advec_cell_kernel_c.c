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
 *  @brief C cell advection kernel.
 *  @author Wayne Gaudin
 *  @details Performs a second order advective remap using van-Leer limiting
 *  with directional splitting.
 */

#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void advec_cell_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                          int *dr,
                          int *swp_nmbr,
                          int *vctr,
                       double *vertexdx,
                       double *vertexdy,
                       double *volume,
                       double *density1,
                       double *energy1,
                       double *mass_flux_x,
                       double *vol_flux_x,
                       double *mass_flux_y,
                       double *vol_flux_y,
                       double *pre_vol,
                       double *post_vol,
                       double *pre_mass,
                       double *post_mass,
                       double *advec_vol,
                       double *post_ener,
                       double *ener_flux)


{
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;
  int sweep_number=*swp_nmbr;
  int dir=*dr;
  int vector=*vctr;

  int j,k,upwind,donor,downwind,dif;

  int g_xdir=1,g_ydir=2;

  double sigma,sigmat,sigmav,sigmam,sigma3,sigma4,diffuw,diffdw,limiter;
  double one_by_six;

  one_by_six=1.0/6.0;

#pragma omp parallel
 {
  if(dir==g_xdir){

    if(sweep_number==1){
#pragma omp for private(j)
      for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
        for (j=x_min-2;j<=x_max+2;j++) {

          pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                           +(vol_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                                            -vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                            +vol_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
                                                            -vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]);
          post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                            -(vol_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                                             -vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]);
        }
      }

    }
    else {
#pragma omp for private(j)
      for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
        for (j=x_min-2;j<=x_max+2;j++) {
          pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                            +vol_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                                            -vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
          post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
        }
      }

    }
#pragma omp for private(upwind,donor,downwind,dif,sigmat,sigma3,sigma4,sigmav,sigma,sigmam,diffuw,diffdw,limiter,j)
    for (k=y_min;k<=y_max;k++) {
      for (j=x_min;j<=x_max+2;j++) {

        if(vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]>0.0){
          upwind   =j-2;
          donor    =j-1;
          downwind =j;
          dif      =donor;
        }
        else {
          upwind   =MIN(j+1,x_max+2);
          donor    =j;
          downwind =j-1;
          dif      =upwind;
        }

        sigmat=fabs(vol_flux_x[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]/pre_vol[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]);
        sigma3=(1.0+sigmat)*(vertexdx[FTNREF1D(j,x_min-2)]/vertexdx[FTNREF1D(dif,x_min-2)]);
        sigma4=2.0-sigmat;

        sigma=sigmat;
        sigmav=sigmat;

        diffuw=density1[FTNREF2D(donor,k  ,x_max+4,x_min-2,y_min-2)]-density1[FTNREF2D(upwind,k  ,x_max+4,x_min-2,y_min-2)];
        diffdw=density1[FTNREF2D(downwind,k  ,x_max+4,x_min-2,y_min-2)]-density1[FTNREF2D(donor,k  ,x_max+4,x_min-2,y_min-2)];
        if(diffuw*diffdw>0.0){
          limiter=(1.0-sigmav)*SIGN(1.0,diffdw)*MIN(fabs(diffuw),MIN(fabs(diffdw)
              ,one_by_six*(sigma3*fabs(diffuw)+sigma4*fabs(diffdw))));
        }
        else{
          limiter=0.0;
        }
        mass_flux_x[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]=vol_flux_x[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]
                                                          *(density1[FTNREF2D(donor,k  ,x_max+4,x_min-2,y_min-2)]+limiter);

        sigmam=fabs(mass_flux_x[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)])/(density1[FTNREF2D(donor,k  ,x_max+4,x_min-2,y_min-2)]
              *pre_vol[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]);
        diffuw=energy1[FTNREF2D(donor,k  ,x_max+4,x_min-2,y_min-2)]-energy1[FTNREF2D(upwind,k  ,x_max+4,x_min-2,y_min-2)];
        diffdw=energy1[FTNREF2D(downwind,k  ,x_max+4,x_min-2,y_min-2)]-energy1[FTNREF2D(donor,k  ,x_max+4,x_min-2,y_min-2)];
        if(diffuw*diffdw>0.0){
          limiter=(1.0-sigmam)*SIGN(1.0,diffdw)*MIN(fabs(diffuw),MIN(fabs(diffdw)
              ,one_by_six*(sigma3*fabs(diffuw)+sigma4*fabs(diffdw))));
        }
        else {
          limiter=0.0;
        }
        ener_flux[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]=mass_flux_x[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]
                                                        *(energy1[FTNREF2D(donor,k  ,x_max+4,x_min-2,y_min-2)]+limiter);

      }
    }
    
#pragma omp for private(j)
    for (k=y_min;k<=y_max;k++) {
#pragma ivdep
     for (j=x_min;j<=x_max;j++) {
        pre_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                           *pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        post_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=pre_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                            +mass_flux_x[FTNREF2D(j  ,k,x_max+5,x_min-2,y_min-2)]
                                                            -mass_flux_x[FTNREF2D(j+1,k,x_max+5,x_min-2,y_min-2)];
        post_ener[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=(energy1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                            *pre_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                            +ener_flux[FTNREF2D(j  ,k,x_max+5,x_min-2,y_min-2)]
                                                            -ener_flux[FTNREF2D(j+1,k,x_max+5,x_min-2,y_min-2)])
                                                            /post_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        advec_vol [FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                             +vol_flux_x[FTNREF2D(j  ,k,x_max+5,x_min-2,y_min-2)]
                                                             -vol_flux_x[FTNREF2D(j+1,k,x_max+5,x_min-2,y_min-2)];

        density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=post_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]/advec_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        energy1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=post_ener[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }

  }
  else if(dir==g_ydir){

    if(sweep_number==1){
      
#pragma omp for private(j)
      for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
        for (j=x_min-2;j<=x_max+2;j++) {

          pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                           +(vol_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
                                                            -vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                            +vol_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                                            -vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]);
          post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                            -(vol_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
                                                             -vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]);
        }
      }

    }
    else {

#pragma omp for private(j)
      for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
        for (j=x_min-2;j<=x_max+2;j++) {
          pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                        +vol_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
                                                        -vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
          post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
        }
      }

    }
    
#pragma omp for private(upwind,donor,downwind,dif,sigmat,sigma3,sigma4,sigmav,sigma,sigmam,diffuw,diffdw,limiter,j)
    for (k=y_min;k<=y_max+2;k++) {
      for (j=x_min;j<=x_max;j++) {

        if(vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]>0.0){
          upwind   =k-2;
          donor    =k-1;
          downwind =k;
          dif      =donor;
        }
        else {
          upwind   =MIN(k+1,y_max+2);
          donor    =k;
          downwind =k-1;
          dif      =upwind;
        }

        sigmat=fabs(vol_flux_y[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)]/pre_vol[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]);
        sigma3=(1.0+sigmat)*(vertexdy[FTNREF1D(k,y_min-2)]/vertexdy[FTNREF1D(dif,y_min-2)]);
        sigma4=2.0-sigmat;

        sigma=sigmat;
        sigmav=sigmat;

        diffuw=density1[FTNREF2D(j  ,donor,x_max+4,x_min-2,y_min-2)]-density1[FTNREF2D(j  ,upwind,x_max+4,x_min-2,y_min-2)];
        diffdw=density1[FTNREF2D(j  ,downwind,x_max+4,x_min-2,y_min-2)]-density1[FTNREF2D(j  ,donor,x_max+4,x_min-2,y_min-2)];

        if(diffuw*diffdw>0.0){
          limiter=(1.0-sigmav)*SIGN(1.0,diffdw)*MIN(fabs(diffuw),MIN(fabs(diffdw)
              ,one_by_six*(sigma3*fabs(diffuw)+sigma4*fabs(diffdw))));
        }
        else{
          limiter=0.0;
        }
        mass_flux_y[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)]=vol_flux_y[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)]
                                                          *(density1[FTNREF2D(j  ,donor,x_max+4,x_min-2,y_min-2)]+limiter);

        sigmam=fabs(mass_flux_y[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)])/(density1[FTNREF2D(j  ,donor,x_max+4,x_min-2,y_min-2)]
              *pre_vol[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]);
        diffuw=energy1[FTNREF2D(j  ,donor,x_max+4,x_min-2,y_min-2)]-energy1[FTNREF2D(j  ,upwind,x_max+4,x_min-2,y_min-2)];
        diffdw=energy1[FTNREF2D(j  ,downwind,x_max+4,x_min-2,y_min-2)]-energy1[FTNREF2D(j  ,donor,x_max+4,x_min-2,y_min-2)];
        if(diffuw*diffdw>0.0){
          limiter=(1.0-sigmam)*SIGN(1.0,diffdw)*MIN(fabs(diffuw),MIN(fabs(diffdw)
              ,one_by_six*(sigma3*fabs(diffuw)+sigma4*fabs(diffdw))));
        }
        else {
          limiter=0.0;
        }
        ener_flux[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]=mass_flux_y[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)]
                                                        *(energy1[FTNREF2D(j  ,donor,x_max+4,x_min-2,y_min-2)]+limiter);

      }
    }
    
#pragma omp for private(j)
    for (k=y_min;k<=y_max;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max;j++) {
        pre_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                          *pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        post_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=pre_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                            +mass_flux_y[FTNREF2D(j,k  ,x_max+4,x_min-2,y_min-2)]
                                                            -mass_flux_y[FTNREF2D(j,k+1,x_max+4,x_min-2,y_min-2)];
        post_ener[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=(energy1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                            *pre_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                            +ener_flux[FTNREF2D(j,k  ,x_max+5,x_min-2,y_min-2)]
                                                            -ener_flux[FTNREF2D(j,k+1,x_max+5,x_min-2,y_min-2)])
                                                            /post_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        advec_vol [FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                             +vol_flux_y[FTNREF2D(j,k  ,x_max+4,x_min-2,y_min-2)]
                                                             -vol_flux_y[FTNREF2D(j,k+1,x_max+4,x_min-2,y_min-2)];

        density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=post_mass[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]/advec_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        energy1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]=post_ener[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }

  }

 }

}

