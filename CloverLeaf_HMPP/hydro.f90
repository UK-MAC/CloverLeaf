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


MODULE hydro_cycle_module

CONTAINS

SUBROUTINE hydro_cycle(c,                 &
                       x_min,             &
                       x_max,             &
                       y_min,             &
                       y_max,             &
                       density0,          &
                       density1,          &
                       energy0,           &
                       energy1,           &
                       pressure,          &
                       soundspeed,        &
                       viscosity,         &
                       xvel0,             &
                       yvel0,             &
                       xvel1,             &
                       yvel1,             &
                       vol_flux_x,        &
                       vol_flux_y,        &
                       mass_flux_x,       &
                       mass_flux_y,       &
                       volume,            &
                       work_array1,       &
                       work_array2,       &
                       work_array3,       &
                       work_array4,       &
                       work_array5,       &
                       work_array6,       &
                       work_array7,       &
                       cellx,             &
                       celly,             &
                       celldx,            &
                       celldy,            &
                       vertexx,           &
                       vertexdx,          &
                       vertexy,           &
                       vertexdy,          &
                       xarea,             &
                       yarea,             &
                       left_snd_buffer,   &
                       left_rcv_buffer,   &
                       right_snd_buffer,  &
                       right_rcv_buffer,  &
                       bottom_snd_buffer, &
                       bottom_rcv_buffer, &
                       top_snd_buffer,    &
                       top_rcv_buffer)

  USE clover_module
  USE timestep_module
  USE viscosity_module
  USE PdV_module
  USE accelerate_module
  USE flux_calc_module
  USE advection_module
  USE reset_field_module

  IMPLICIT NONE
  INTEGER               :: c,x_min,x_max,y_min,y_max

  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: density0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: density1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: energy1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: soundspeed
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: pressure
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: viscosity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: xvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: yvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: xvel1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: yvel1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+2) :: vol_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+3) :: vol_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+2) :: mass_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+3) :: mass_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: work_array1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: work_array2
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: work_array3
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: work_array4
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: work_array5
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: work_array6
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+3) :: work_array7
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: cellx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celly
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celldy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexdx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexy
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexdy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+2) :: xarea
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+3) :: yarea
  REAL(KIND=8) :: left_snd_buffer(:),left_rcv_buffer(:),right_snd_buffer(:),right_rcv_buffer(:)
  REAL(KIND=8) :: bottom_snd_buffer(:),bottom_rcv_buffer(:),top_snd_buffer(:),top_rcv_buffer(:)

  INTEGER         :: cells
  REAL(KIND=8)    :: timer,timerstart
  
  REAL(KIND=8)    :: grind_time
  REAL(KIND=8)    :: step_time,step_grind

  timerstart = timer()

!$ACC DATA &
!$ACC COPYIN(density0)   &
!$ACC COPYIN(density1)   &
!$ACC COPYIN(energy0)    &
!$ACC COPYIN(energy1)    &
!$ACC COPYIN(pressure)   &
!$ACC COPYIN(soundspeed) &
!$ACC COPYIN(viscosity)  &
!$ACC COPYIN(xvel0)      &
!$ACC COPYIN(yvel0)      &
!$ACC COPYIN(xvel1)      &
!$ACC COPYIN(yvel1)      &
!$ACC COPYIN(vol_flux_x) &
!$ACC COPYIN(vol_flux_y) &
!$ACC COPYIN(mass_flux_x)&
!$ACC COPYIN(mass_flux_y)&
!$ACC COPYIN(volume)     &
!$ACC COPYIN(work_array1)&
!$ACC COPYIN(work_array2)&
!$ACC COPYIN(work_array3)&
!$ACC COPYIN(work_array4)&
!$ACC COPYIN(work_array5)&
!$ACC COPYIN(work_array6)&
!$ACC COPYIN(work_array7)&
!$ACC COPYIN(cellx)      &
!$ACC COPYIN(celly)      &
!$ACC COPYIN(celldx)     &
!$ACC COPYIN(celldy)     &
!$ACC COPYIN(vertexx)    &
!$ACC COPYIN(vertexdx)   &
!$ACC COPYIN(vertexy)    &
!$ACC COPYIN(vertexdy)   &
!$ACC COPYIN(xarea)      &
!$ACC COPYIN(yarea)      &
!$ACC COPY(left_snd_buffer)    &
!$ACC COPY(left_rcv_buffer)    &
!$ACC COPY(right_snd_buffer)   &
!$ACC COPY(right_rcv_buffer)   &
!$ACC COPY(bottom_snd_buffer)  &
!$ACC COPY(bottom_rcv_buffer)  &
!$ACC COPY(top_snd_buffer)     &
!$ACC COPY(top_rcv_buffer)

  DO

    step_time = timer()

    step = step + 1

    CALL timestep(c)

    CALL PdV(c,.TRUE.)

    CALL accelerate(c)

    CALL PdV(c,.FALSE.)

    CALL flux_calc(c)

    CALL advection(c)

    CALL reset_field(c)

    advect_x = .NOT. advect_x
  
    time = time + dt

    IF(summary_frequency.NE.0) THEN
      IF(MOD(step, summary_frequency).EQ.0) CALL field_summary(c)
    ENDIF
    IF(visit_frequency.NE.0) THEN
      IF(MOD(step, visit_frequency).EQ.0) CALL visit(c)
    ENDIF

    IF(time+g_small.GT.end_time.OR.step.GE.end_step) THEN

      complete=.TRUE.
      CALL field_summary(c)
      IF(visit_frequency.NE.0) CALL visit(c)

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

END SUBROUTINE hydro_cycle

END MODULE hydro_cycle_module

SUBROUTINE hydro

  USE clover_module
  USE hydro_cycle_module

  IMPLICIT NONE

  INTEGER         :: cells
  REAL(KIND=8)    :: timer,timerstart
  
  REAL(KIND=8)    :: grind_time
  REAL(KIND=8)    :: step_time,step_grind

  timerstart = timer()

  CALL hydro_cycle(parallel%task+1,                          &
                   chunks(parallel%task+1)%field%x_min,      &
                   chunks(parallel%task+1)%field%x_max,      &
                   chunks(parallel%task+1)%field%y_min,      &
                   chunks(parallel%task+1)%field%y_max,      &
                   chunks(parallel%task+1)%field%density0,   &
                   chunks(parallel%task+1)%field%density1,   &
                   chunks(parallel%task+1)%field%energy0,    &
                   chunks(parallel%task+1)%field%energy1,    &
                   chunks(parallel%task+1)%field%pressure,   &
                   chunks(parallel%task+1)%field%soundspeed, &
                   chunks(parallel%task+1)%field%viscosity,  &
                   chunks(parallel%task+1)%field%xvel0,      &
                   chunks(parallel%task+1)%field%yvel0,      &
                   chunks(parallel%task+1)%field%xvel1,      &
                   chunks(parallel%task+1)%field%yvel1,      &
                   chunks(parallel%task+1)%field%vol_flux_x, &
                   chunks(parallel%task+1)%field%vol_flux_y, &
                   chunks(parallel%task+1)%field%mass_flux_x,&
                   chunks(parallel%task+1)%field%mass_flux_y,&
                   chunks(parallel%task+1)%field%volume,     &
                   chunks(parallel%task+1)%field%work_array1,&
                   chunks(parallel%task+1)%field%work_array2,&
                   chunks(parallel%task+1)%field%work_array3,&
                   chunks(parallel%task+1)%field%work_array4,&
                   chunks(parallel%task+1)%field%work_array5,&
                   chunks(parallel%task+1)%field%work_array6,&
                   chunks(parallel%task+1)%field%work_array7,&
                   chunks(parallel%task+1)%field%cellx,      &
                   chunks(parallel%task+1)%field%celly,      &
                   chunks(parallel%task+1)%field%celldx,     &
                   chunks(parallel%task+1)%field%celldy,     &
                   chunks(parallel%task+1)%field%vertexx,    &
                   chunks(parallel%task+1)%field%vertexdx,   &
                   chunks(parallel%task+1)%field%vertexy,    &
                   chunks(parallel%task+1)%field%vertexdy,   &
                   chunks(parallel%task+1)%field%xarea,      &
                   chunks(parallel%task+1)%field%yarea,      &
                   chunks(parallel%task+1)%left_snd_buffer,  &
                   chunks(parallel%task+1)%left_rcv_buffer,  &
                   chunks(parallel%task+1)%right_snd_buffer, &
                   chunks(parallel%task+1)%right_rcv_buffer, &
                   chunks(parallel%task+1)%bottom_snd_buffer,&
                   chunks(parallel%task+1)%bottom_rcv_buffer,&
                   chunks(parallel%task+1)%top_snd_buffer,   &
                   chunks(parallel%task+1)%top_rcv_buffer)

END SUBROUTINE hydro

