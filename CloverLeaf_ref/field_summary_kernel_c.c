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

void field_summary_kernel(x_min,x_max,y_min,y_max,
                                volume,
                                density0,
                                energy0,
                                pressure,
                                xvel0,
                                yvel0,
                                vol,mass,ie,ke,press)
{
  int x_min,x_max,y_min,y_max;
  double volume;
  double density0,energy0;
  double pressure;
  double xvel0,yvel0;
  double vol,mass,ie,ke,press;

  int j,k,jv,kv;
  double vsqrd,cell_vol,cell_mass;

  vol=0.0;
  mass=0.0;
  ie=0.0;
  ke=0.0;;
  press=0.0;

!$OMP PARALLEL
!$OMP DO PRIVATE(vsqrd,cell_vol,cell_mass) REDUCTION(+ : vol,mass,press,ie,ke,j)
  DO k=y_min,y_max
    DO j=x_min,x_max
      vsqrd=0.0
      DO kv=k,k+1
        DO jv=j,j+1
          vsqrd=vsqrd+0.25*(xvel0(jv,kv)**2+yvel0(jv,kv)**2)
        ENDDO
      ENDDO
      cell_vol=volume(j,k)
      cell_mass=cell_vol*density0(j,k)
      vol=vol+cell_vol
      mass=mass+cell_mass
      ie=ie+cell_mass*energy0(j,k)
      ke=ke+cell_mass*0.5*vsqrd
      press=press+cell_vol*pressure(j,k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

}
