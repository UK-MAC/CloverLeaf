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

!>  @brief Reset field driver
!>  @author Wayne Gaudin
!>  @details Invokes the user specified field reset kernel.

MODULE reset_field_module

CONTAINS

SUBROUTINE reset_field()

  USE clover_module
  USE reset_field_kernel_module

  IMPLICIT NONE

  INTEGER :: c

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels)THEN
        CALL reset_field_kernel(chunks(c)%field%x_min,   &
                              chunks(c)%field%x_max,     &
                              chunks(c)%field%y_min,     &
                              chunks(c)%field%y_max,     &
                              chunks(c)%field%density0,  &
                              chunks(c)%field%density1,  &
                              chunks(c)%field%energy0,   &
                              chunks(c)%field%energy1,   &
                              chunks(c)%field%xvel0,     &
                              chunks(c)%field%xvel1,     &
                              chunks(c)%field%yvel0,     &
                              chunks(c)%field%yvel1      )
      ELSEIF(use_C_kernels)THEN
        CALL reset_field_kernel_c(chunks(c)%field%x_min, &
                              chunks(c)%field%x_max,     &
                              chunks(c)%field%y_min,     &
                              chunks(c)%field%y_max,     &
                              chunks(c)%field%density0,  &
                              chunks(c)%field%density1,  &
                              chunks(c)%field%energy0,   &
                              chunks(c)%field%energy1,   &
                              chunks(c)%field%xvel0,     &
                              chunks(c)%field%xvel1,     &
                              chunks(c)%field%yvel0,     &
                              chunks(c)%field%yvel1      )
      ENDIF

    ENDIF

  ENDDO

END SUBROUTINE reset_field

END MODULE reset_field_module
