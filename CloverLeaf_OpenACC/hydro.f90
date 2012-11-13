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

!>  @brief Controls the main hydro cycle.
!>  @author Wayne Gaudin
!>  @details Controls the top level cycle, invoking all the drivers and checks
!>  for outputs and completion.

SUBROUTINE hydro

  USE clover_module
  USE timestep_module
  USE viscosity_module
  USE PdV_module
  USE accelerate_module
  USE flux_calc_module
  USE advection_module
  USE reset_field_module

  IMPLICIT NONE

  INTEGER         :: cells
  REAL(KIND=8)    :: timer,timerstart
  
  REAL(KIND=8)    :: grind_time
  REAL(KIND=8)    :: step_time,step_grind

  timerstart = timer()

!$ACC DATA &
!$ACC COPYIN(chunks(parallel%task+1)%field%density0)   &
!$ACC COPYIN(chunks(parallel%task+1)%field%density1)   &
!$ACC COPYIN(chunks(parallel%task+1)%field%energy0)    &
!$ACC COPYIN(chunks(parallel%task+1)%field%energy1)    &
!$ACC COPYIN(chunks(parallel%task+1)%field%pressure)   &
!$ACC COPYIN(chunks(parallel%task+1)%field%soundspeed) &
!$ACC COPYIN(chunks(parallel%task+1)%field%viscosity)  &
!$ACC COPYIN(chunks(parallel%task+1)%field%xvel0)      &
!$ACC COPYIN(chunks(parallel%task+1)%field%yvel0)      &
!$ACC COPYIN(chunks(parallel%task+1)%field%xvel1)      &
!$ACC COPYIN(chunks(parallel%task+1)%field%yvel1)      &
!$ACC COPYIN(chunks(parallel%task+1)%field%vol_flux_x) &
!$ACC COPYIN(chunks(parallel%task+1)%field%vol_flux_y) &
!$ACC COPYIN(chunks(parallel%task+1)%field%mass_flux_x)&
!$ACC COPYIN(chunks(parallel%task+1)%field%mass_flux_y)&
!$ACC COPYIN(chunks(parallel%task+1)%field%volume)     &
!$ACC COPYIN(chunks(parallel%task+1)%field%work_array1)&
!$ACC COPYIN(chunks(parallel%task+1)%field%work_array2)&
!$ACC COPYIN(chunks(parallel%task+1)%field%work_array3)&
!$ACC COPYIN(chunks(parallel%task+1)%field%work_array4)&
!$ACC COPYIN(chunks(parallel%task+1)%field%work_array5)&
!$ACC COPYIN(chunks(parallel%task+1)%field%work_array6)&
!$ACC COPYIN(chunks(parallel%task+1)%field%work_array7)&
!$ACC COPYIN(chunks(parallel%task+1)%field%cellx)      &
!$ACC COPYIN(chunks(parallel%task+1)%field%celly)      &
!$ACC COPYIN(chunks(parallel%task+1)%field%celldx)     &
!$ACC COPYIN(chunks(parallel%task+1)%field%celldy)     &
!$ACC COPYIN(chunks(parallel%task+1)%field%vertexx)    &
!$ACC COPYIN(chunks(parallel%task+1)%field%vertexdx)   &
!$ACC COPYIN(chunks(parallel%task+1)%field%vertexy)    &
!$ACC COPYIN(chunks(parallel%task+1)%field%vertexdy)   &
!$ACC COPYIN(chunks(parallel%task+1)%field%xarea)      &
!$ACC COPYIN(chunks(parallel%task+1)%field%yarea)      &
!$ACC COPY(chunks(parallel%task+1)%left_snd_buffer)    &
!$ACC COPY(chunks(parallel%task+1)%left_rcv_buffer)    &
!$ACC COPY(chunks(parallel%task+1)%right_snd_buffer)   &
!$ACC COPY(chunks(parallel%task+1)%right_rcv_buffer)   &
!$ACC COPY(chunks(parallel%task+1)%bottom_snd_buffer)  &
!$ACC COPY(chunks(parallel%task+1)%bottom_rcv_buffer)  &
!$ACC COPY(chunks(parallel%task+1)%top_snd_buffer)     &
!$ACC COPY(chunks(parallel%task+1)%top_rcv_buffer)
  DO

    step_time = timer()

    step = step + 1

    CALL timestep()

    CALL PdV(.TRUE.)

    CALL accelerate()

    CALL PdV(.FALSE.)

    CALL flux_calc()

    CALL advection()

    CALL reset_field()

    advect_x = .NOT. advect_x
  
    time = time + dt

    IF(summary_frequency.NE.0) THEN
      IF(MOD(step, summary_frequency).EQ.0) CALL field_summary()
    ENDIF
    IF(visit_frequency.NE.0) THEN
      IF(MOD(step, visit_frequency).EQ.0) CALL visit()
    ENDIF

    IF(time+g_small.GT.end_time.OR.step.GE.end_step) THEN

      complete=.TRUE.
      CALL field_summary()
      IF(visit_frequency.NE.0) CALL visit()

      IF ( parallel%boss ) THEN
        WRITE(g_out,*)
        WRITE(g_out,*) 'Calculation complete'
        WRITE(g_out,*) 'Clover is finishing'
        WRITE(g_out,*) 'Wall clock ', timer() - timerstart
        WRITE(    0,*) 'Wall clock ', timer() - timerstart
      ENDIF

      CALL clover_finalize

      EXIT

    END IF

    IF (parallel%boss) THEN
      WRITE(g_out,*)"Wall clock ",timer()-timerstart
      WRITE(0    ,*)"Wall clock ",timer()-timerstart
      cells = grid%x_cells * grid%y_cells
      grind_time   = (timer() - timerstart) / (step * cells)
      step_grind   = (timer() - step_time)/cells
      WRITE(0    ,*)"Average time per cell ",grind_time
      WRITE(g_out,*)"Average time per cell ",grind_time
      WRITE(0    ,*)"Step time per cell    ",step_grind
      WRITE(g_out,*)"Step time per cell    ",step_grind

     END IF

  END DO
!$ACC END DATA

END SUBROUTINE hydro
