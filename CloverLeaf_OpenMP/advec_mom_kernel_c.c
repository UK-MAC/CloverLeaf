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
 *  @brief C momentum advection kernel
 *  @author Wayne Gaudin
 *  @details Performs a second order advective remap on the vertex momentum
 *  using van-Leer limiting and directional splitting.
 *  Note that although pre_vol is only set and not used in the update, please
 *  leave it in the method.
 */

#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void advec_mom_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
                      double *xvel1,
                      double *yvel1,
                      double *mass_flux_x,
                      double *vol_flux_x,
                      double *mass_flux_y,
                      double *vol_flux_y,
                      double *volume,
                      double *density1,
                      double *node_flux,
                      double *node_mass_post,
                      double *node_mass_pre,
                      double *advec_vel,
                      double *mom_flux,
                      double *pre_vol,
                      double *post_vol,
                      double *celldx,
                      double *celldy,
                         int *whch_vl,
                         int *swp_nmbr,
                         int *drctn,
                         int *vctr)

{
  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;
  int which_vel=*whch_vl;
  int sweep_number=*swp_nmbr;
  int direction=*drctn;
  int vector=*vctr;

  int j,k,mom_sweep;
  int upwind,donor,downwind,dif;
  double sigma,wind,width;
  double sigma2,wind2;
  double vdiffuw,vdiffdw,auw,adw,limiter;
  double vdiffuw2,vdiffdw2,auw2,limiter2;

  double *vel1;

  if(which_vel==1){
    vel1=xvel1;
  } else{
    vel1=yvel1;
  }

  mom_sweep=direction+2*(sweep_number-1);
  
#pragma omp parallel
 {
  if(mom_sweep==1){
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
      for (j=x_min-2;j<=x_max+2;j++) {
        post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                           +vol_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
                                                           -vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
        pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                          +vol_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                                          -vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
  } else if(mom_sweep==2){
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
      for (j=x_min-2;j<=x_max+2;j++) {
        post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                                +vol_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                                                -vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                               +vol_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
                                                               -vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
      }
    }
  } else if(mom_sweep==3){
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
      for (j=x_min-2;j<=x_max+2;j++) {
        post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
        pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                          +vol_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]
                                                          -vol_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
      }
    }
  } else if(mom_sweep==4){
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
      for (j=x_min-2;j<=x_max+2;j++) {
        post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=volume[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)];
        pre_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                          +vol_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]
                                                          -vol_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
  }

  if(direction==1) {
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
      for (j=x_min-2;j<=x_max+2;j++) {
        node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                            *(mass_flux_x[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                             +mass_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                             +mass_flux_x[FTNREF2D(j+1,k-1,x_max+5,x_min-2,y_min-2)]
                                                             +mass_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
      for (j=x_min-1;j<=x_max+2;j++) {
        node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                                 *(density1[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
      for (j=x_min-1;j<=x_max+2;j++) {
        node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                          -node_flux[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]+node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
    if (vector==1) {
#pragma omp for private(sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind,sigma2,limiter2,vdiffuw2,vdiffdw2,auw2,wind2)
      for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
        for (j=x_min-1;j<=x_max+1;j++) {
          sigma=fabs(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])/(node_mass_pre[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]);
          sigma2=fabs(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])/(node_mass_pre[FTNREF2D(j,k  ,x_max+5,x_min-2,y_min-2)]);
          width=celldx[FTNREF1D(j,x_min-2)];
          vdiffuw=vel1[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j+2,k  ,x_max+5,x_min-2,y_min-2)];
          vdiffdw=vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)];
          vdiffuw2=vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)];
          vdiffdw2=-vdiffdw;
          auw=fabs(vdiffuw);
          adw=fabs(vdiffdw);
          auw2=fabs(vdiffuw2);
          wind=1.0;
          wind2=1.0;
          limiter=0.0;
          if(vdiffdw<=0.0) wind=-1.0;
          if(vdiffdw2<=0.0) wind2=-1.0;
          limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldx[FTNREF1D(j+1,x_min-2)])/6.0,MIN(auw,adw));
          limiter2=wind2*MIN(width*((2.0-sigma2)*adw/width+(1.0+sigma2)*auw2/celldx[FTNREF1D(j-1,x_min-2)])/6.0,MIN(auw2,adw));
          if(vdiffuw*vdiffdw<=0.0) limiter=0.0;
          if(vdiffuw2*vdiffdw2<=0.0) limiter2=0.0;
          if(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]<0.0){
            advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=vel1[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]+(1.0-sigma)*limiter;
          }
          else{
            advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]+(1.0-sigma2)*limiter2;
          }
          mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                             *node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        }
      }
    }else{
#pragma omp for private(upwind,downwind,donor,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind,j)
      for (k=y_min;k<=y_max+1;k++) {
        for (j=x_min-1;j<=x_max+1;j++) {
          if(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]<0.0){
            upwind=j+2;
            donor=j+1;
            downwind=j;
            dif=donor;
          }
          else{
            upwind=j-1;
            donor=j;
            downwind=j+1;
            dif=upwind;
          }
          sigma=fabs(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])/(node_mass_pre[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]);
          width=celldx[FTNREF1D(j,x_min-2)];
          vdiffuw=vel1[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(upwind,k  ,x_max+5,x_min-2,y_min-2)];
          vdiffdw=vel1[FTNREF2D(downwind,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)];
          limiter=0.0;
          if(vdiffuw*vdiffdw>0.0){
            auw=fabs(vdiffuw);
            adw=fabs(vdiffdw);
            wind=1.0;
            if(vdiffdw<=0.0) wind=-1.0;
            limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldx[FTNREF1D(dif,x_min-2)])/6.0,MIN(auw,adw));
          }
          advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=vel1[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]+(1.0-sigma)*limiter;
          mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                             *node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        }
      }
    }

#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max+1;j++) {
        vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=(vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        *node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        +mom_flux[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]
                                                        -mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])
                                                        /node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
  }
  else if(direction==2){
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max+1;j++) {
        node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                            *(mass_flux_y[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                                             +mass_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                             +mass_flux_y[FTNREF2D(j-1,k+1,x_max+4,x_min-2,y_min-2)]
                                                             +mass_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min-1;k<=y_max+2;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max+1;j++) {
        node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                                 *(density1[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min-1;k<=y_max+2;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max+1;j++) {
        node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                          -node_flux[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]+node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
    if (vector==1) {
#pragma omp for private(sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind,sigma2,limiter2,vdiffuw2,vdiffdw2,auw2,wind2)
      for (k=y_min-1;k<=y_max+1;k++) {
#pragma ivdep
        for (j=x_min;j<=x_max+1;j++) {
          sigma=fabs(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])/(node_mass_pre[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]);
          sigma2=fabs(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])/(node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]);
          width=celldy[FTNREF1D(k,x_min-2)];
          vdiffuw=vel1[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j  ,k+2,x_max+5,x_min-2,y_min-2)];
          vdiffdw=vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)];
          vdiffuw2=vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)];
          vdiffdw2=-vdiffdw;
          auw=fabs(vdiffuw);
          adw=fabs(vdiffdw);
          auw2=fabs(vdiffuw2);
          wind=1.0;
          wind2=1.0;
          limiter=0.0;
          if(vdiffdw<=0.0) wind=-1.0;
          if(vdiffdw2<=0.0) wind2=-1.0;
          limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldy[FTNREF1D(k+1,x_min-2)])/6.0,MIN(auw,adw));
          limiter2=wind2*MIN(width*((2.0-sigma2)*adw/width+(1.0+sigma2)*auw2/celldy[FTNREF1D(k-1,x_min-2)])/6.0,MIN(auw2,adw));
          if(vdiffuw*vdiffdw<=0.0) limiter=0.0;
          if(vdiffuw2*vdiffdw2<=0.0) limiter2=0.0;
          if(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]<0.0){
            advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=vel1[FTNREF2D(j  ,k+1,x_max+5,x_min-2,y_min-2)]+(1.0-sigma)*limiter;
          }
          else{
            advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]+(1.0-sigma2)*limiter2;
          }
          mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                             *node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        }
      }
    }else{
#pragma omp for private(upwind,downwind,donor,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind,j)
      for (k=y_min-1;k<=y_max+1;k++) {
        for (j=x_min;j<=x_max+1;j++) {
          if(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]<0.0){
            upwind=k+2;
            donor=k+1;
            downwind=k;
            dif=donor;
          }
          else{
            upwind=k-1;
            donor=k;
            downwind=k+1;
            dif=upwind;
          }
          sigma=fabs(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])/(node_mass_pre[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]);
          width=celldy[FTNREF1D(k,y_min-2)];
          vdiffuw=vel1[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j  ,upwind,x_max+5,x_min-2,y_min-2)];
          vdiffdw=vel1[FTNREF2D(j  ,downwind ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)];
          limiter=0.0;
          if(vdiffuw*vdiffdw>0.0){
            auw=fabs(vdiffuw);
            adw=fabs(vdiffdw);
            wind=1.0;
            if(vdiffdw<=0.0) wind=-1.0;
            limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldy[FTNREF1D(dif,y_min-2)])/6.0,MIN(auw,adw));
          }
          advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=vel1[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]+(1.0-sigma)*limiter;
          mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                             *node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
        }
      }
    }
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max+1;j++) {
        vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=(vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        *node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        +mom_flux[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                        -mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])
                                                        /node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }

  }

 }

}

