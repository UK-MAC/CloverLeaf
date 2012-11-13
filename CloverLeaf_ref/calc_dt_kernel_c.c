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
 *  @brief Not yet called
 *  @author Wayne Gaudin
 *  @details Still just a stub
 */

#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

void calc_dt_kernel(x_min,x_max,y_min,y_max,
                    g_small,g_big,dtmin,
                    dtc_safe,
                    dtu_safe,
                    dtv_safe,
                    dtdiv_safe,
                    xarea,
                    yarea,
                    cellx,
                    celly,
                    celldx,
                    celldy,
                    volume,
                    density0,
                    energy0,
                    pressure,
                    viscosity,
                    soundspeed,
                    xvel0,yvel0,
                    dt_min,
                    dt_min_val,
                    dtl_control,
                    xl_pos,
                    yl_pos,
                    jldt,
                    kldt,
                    small)
{

  int  x_min,x_max,y_min,y_max;
  double g_small,g_big,dtmin,dt_min_val;
  double dtc_safe,dtu_safe,dtv_safe,dtdiv_safe;
  double xarea;
  double yarea;
  double celly;
  double celldy;
  double volume;
  double density0;
  double energy0;
  double pressure;
  double viscosity;
  double soundspeed;
  double xvel0,yvel0;
  double dt_min;

  int dtl_control;
  double xl_pos,yl_pos;
  int jldt,kldt;
  int small;

  int j,k;

  double div,dsx,dsy,dtut,dtvt,dtct,dtdivt,cc,dv1,dv2,jk_control;

!$OMP PARALLEL

  small=0

  dt_min_val = g_big
  jk_control=1.1

!$OMP DO PRIVATE(dsx,dsy,cc,dv1,dv2,div,dtct,dtut,dtvt,dtdivt,j)
  DO k=y_min,y_max
    DO j=x_min,x_max

       dsx=celldx(j)
       dsy=celldy(k)

       cc=soundspeed(j,k)**2
       cc=cc+2.0*viscosity(j,k)/density0(j,k)
       cc=MAX(SQRT(cc),g_small)

       dtct=dtc_safe*MIN(dsx,dsy)/cc

       div=0.0

       dv1=(xvel0(j  ,k)+xvel0(j  ,k+1))*xarea(j  ,k)
       dv2=(xvel0(j+1,k)+xvel0(j+1,k+1))*xarea(j+1,k)

       div=div+dv2-dv1

       dtut=dtu_safe*2.0*volume(j,k)/MAX(ABS(dv1),ABS(dv2),g_small*volume(j,k))

       dv1=(yvel0(j,k  )+yvel0(j+1,k  ))*yarea(j,k  )
       dv2=(yvel0(j,k+1)+yvel0(j+1,k+1))*yarea(j,k+1)

       div=div+dv2-dv1

       dtvt=dtv_safe*2.0*volume(j,k)/MAX(ABS(dv1),ABS(dv2),g_small*volume(j,k))

       div=div/(2.0*volume(j,k))

       IF(div.LT.-g_small)THEN
         dtdivt=dtdiv_safe*(-1.0/div)
       ELSE
         dtdivt=g_big
       ENDIF

       dt_min(j,k)=MIN(dtct,dtut,dtvt,dtdivt)

    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO REDUCTION(MIN : dt_min_val)
  DO k=y_min,y_max
    DO j=x_min,x_max
      IF(dt_min(j,k).LT.dt_min_val) dt_min_val=dt_min(j,k)
    ENDDO
  ENDDO
!$OMP END DO

!$OMP END PARALLEL

  ! Extract the mimimum timestep information
  dtl_control=10.01*(jk_control-INT(jk_control))
  jk_control=jk_control-(jk_control-INT(jk_control))
  jldt=MOD(INT(jk_control),x_max)
  kldt=1+(jk_control/x_max)
  xl_pos=cellx(jldt)
  yl_pos=celly(kldt)

  IF(dt_min_val.LT.dtmin) small=1

  IF(small.NE.0)THEN
    WRITE(0,*) 'Timestep information:'
    WRITE(0,*) 'j, k                 : ',jldt,kldt
    WRITE(0,*) 'x, y                 : ',cellx(jldt),celly(kldt)
    WRITE(0,*) 'timestep : ',dt_min_val
    WRITE(0,*) 'Cell velocities;'
    WRITE(0,*) xvel0(jldt  ,kldt  ),yvel0(jldt  ,kldt  )
    WRITE(0,*) xvel0(jldt+1,kldt  ),yvel0(jldt+1,kldt  )
    WRITE(0,*) xvel0(jldt+1,kldt+1),yvel0(jldt+1,kldt+1)
    WRITE(0,*) xvel0(jldt  ,kldt+1),yvel0(jldt  ,kldt+1)
    WRITE(0,*) 'density, energy, pressure, soundspeed '
    WRITE(0,*) density0(jldt,kldt),energy0(jldt,kldt),pressure(jldt,kldt),soundspeed(jldt,kldt)
  ENDIF

}
