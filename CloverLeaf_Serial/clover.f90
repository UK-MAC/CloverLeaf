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

!>  @brief Communication Utilities
!>  @author Wayne Gaudin
!>  @details Contains all utilities required to run CloverLeaf in a distributed
!>  environment, including initialisation, mesh decompostion, reductions and
!>  halo exchange using explicit buffers.
!>
!>  Note the halo exchange is currently coded as simply as possible and no 
!>  optimisations have been implemented, such as post receives before sends or packing
!>  buffers with multiple data fields. This is intentional so the effect of these
!>  optimisations can be measured on large systems, as and when they are added.
!>
!>  Even without these modifications CloverLeaf weak scales well on moderately sized
!>  systems of the order of 10K cores.

MODULE clover_module

  USE data_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE clover_barrier

  INTEGER :: err


END SUBROUTINE clover_barrier

SUBROUTINE clover_abort

  INTEGER :: ierr,err


END SUBROUTINE clover_abort

SUBROUTINE clover_finalize

  INTEGER :: err

  CLOSE(g_out)
  CALL FLUSH(0)
  CALL FLUSH(6)
  CALL FLUSH(g_out)

END SUBROUTINE clover_finalize

SUBROUTINE clover_init_comms

  IMPLICIT NONE

  INTEGER :: err,rank,size

  rank=0
  size=1



  parallel%parallel=.TRUE.
  parallel%task=rank

  IF(rank.EQ.0) THEN
    parallel%boss=.TRUE.
  ENDIF

  parallel%boss_task=0
  parallel%max_task=size

END SUBROUTINE clover_init_comms

SUBROUTINE clover_get_num_chunks(count)

  IMPLICIT NONE

  INTEGER :: count

! Should be changed so there can be more than one chunk per mpi task

  count=parallel%max_task

END SUBROUTINE clover_get_num_chunks

SUBROUTINE clover_decompose(x_cells,y_cells,left,right,bottom,top)

  ! This decomposes the mesh into a number of chunks.
  ! The number of chunks may be a multiple of the number of mpi tasks
  ! Doesn't always return the best split if there are few factors
  ! All factors need to be stored and the best picked. But its ok for now

  IMPLICIT NONE

  INTEGER :: x_cells,y_cells,left(:),right(:),top(:),bottom(:)
  INTEGER :: c,delta_x,delta_y

  REAL(KIND=8) :: mesh_ratio,factor_x,factor_y
  INTEGER  :: chunk_x,chunk_y,mod_x,mod_y,split_found

  INTEGER  :: cx,cy,chunk,add_x,add_y,add_x_prev,add_y_prev

  ! 2D Decomposition of the mesh

  mesh_ratio=real(x_cells)/real(y_cells)

  chunk_x=number_of_chunks
  chunk_y=1

  split_found=0 ! Used to detect 1D decomposition
  DO c=1,number_of_chunks
    IF (MOD(number_of_chunks,c).EQ.0) THEN
      factor_x=number_of_chunks/real(c)
      factor_y=c
      !Compare the factor ratio with the mesh ratio
      IF(factor_x/factor_y.LE.mesh_ratio) THEN
        chunk_y=c
        chunk_x=number_of_chunks/c
        split_found=1
        EXIT
      ENDIF
    ENDIF
  ENDDO

  IF(split_found.EQ.0.OR.chunk_y.EQ.number_of_chunks) THEN ! Prime number or 1D decomp detected
    IF(mesh_ratio.GE.1.0) THEN
      chunk_x=number_of_chunks
      chunk_y=1
    ELSE
      chunk_x=1
      chunk_y=number_of_chunks
    ENDIF
  ENDIF

  delta_x=x_cells/chunk_x
  delta_y=y_cells/chunk_y
  mod_x=MOD(x_cells,chunk_x)
  mod_y=MOD(y_cells,chunk_y)

  ! Set up chunk mesh ranges and chunk connectivity

  add_x_prev=0
  add_y_prev=0
  chunk=1
  DO cy=1,chunk_y
    DO cx=1,chunk_x
      add_x=0
      add_y=0
      IF(cx.LE.mod_x)add_x=1
      IF(cy.LE.mod_y)add_y=1
      left(chunk)=(cx-1)*delta_x+1+add_x_prev
      right(chunk)=left(chunk)+delta_x-1+add_x
      bottom(chunk)=(cy-1)*delta_y+1+add_y_prev
      top(chunk)=bottom(chunk)+delta_y-1+add_y
      chunks(chunk)%chunk_neighbours(chunk_left)=chunk_x*(cy-1)+cx-1
      chunks(chunk)%chunk_neighbours(chunk_right)=chunk_x*(cy-1)+cx+1
      chunks(chunk)%chunk_neighbours(chunk_bottom)=chunk_x*(cy-2)+cx
      chunks(chunk)%chunk_neighbours(chunk_top)=chunk_x*(cy)+cx
      IF(cx.EQ.1)chunks(chunk)%chunk_neighbours(chunk_left)=external_face
      IF(cx.EQ.chunk_x)chunks(chunk)%chunk_neighbours(chunk_right)=external_face
      IF(cy.EQ.1)chunks(chunk)%chunk_neighbours(chunk_bottom)=external_face
      IF(cy.EQ.chunk_y)chunks(chunk)%chunk_neighbours(chunk_top)=external_face
      IF(cx.LE.mod_x)add_x_prev=add_x_prev+1
      chunk=chunk+1
    ENDDO
    add_x_prev=0
    IF(cy.LE.mod_y)add_y_prev=add_y_prev+1
  ENDDO

  IF(parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Mesh ratio of ",mesh_ratio
    WRITE(g_out,*)"Decomposing the mesh into ",chunk_x," by ",chunk_y," chunks"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE clover_decompose

SUBROUTINE clover_allocate_buffers(chunk)

  IMPLICIT NONE

  INTEGER      :: chunk
  
  ! Unallocated buffers for external boundaries caused issues on some systems so they are now
  !  all allocated
  IF(parallel%task.EQ.chunks(chunk)%task)THEN
    !IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
    !ENDIF
  ENDIF

END SUBROUTINE clover_allocate_buffers

SUBROUTINE clover_exchange(fields,depth)

  IMPLICIT NONE

  INTEGER      :: fields(:),depth

  ! Assuming 1 patch per task, this will be changed
  ! Also, not packing all fields for each communication, doing one at a time

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%density0,      &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%density1,      &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%energy0,       &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%energy1,       &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%pressure,      &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%viscosity,     &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%soundspeed,    &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_XVEL0).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%xvel0,         &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_XVEL1).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%xvel1,         &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_YVEL0).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%yvel0,         &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_YVEL1).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%yvel1,         &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%vol_flux_x,    &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,X_FACE_DATA)
  ENDIF

  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%vol_flux_y,    &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,Y_FACE_DATA)
  ENDIF

  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%mass_flux_x,   &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,X_FACE_DATA)
  ENDIF

  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%mass_flux_y,   &
                                 chunks(parallel%task+1)%left_snd_buffer,                      &
                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
                                 chunks(parallel%task+1)%right_snd_buffer,                     &
                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
                                 chunks(parallel%task+1)%top_snd_buffer,                       &
                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
                                 depth,Y_FACE_DATA)
  ENDIF


END SUBROUTINE clover_exchange

SUBROUTINE clover_exchange_message(chunk,field,                            &
                                   left_snd_buffer,                        &
                                   left_rcv_buffer,                        &
                                   right_snd_buffer,                       &
                                   right_rcv_buffer,                       &
                                   bottom_snd_buffer,                      &
                                   bottom_rcv_buffer,                      &
                                   top_snd_buffer,                         &
                                   top_rcv_buffer,                         &
                                   depth,field_type)

  IMPLICIT NONE

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  REAL(KIND=8) :: left_snd_buffer(:),left_rcv_buffer(:),right_snd_buffer(:),right_rcv_buffer(:)
  REAL(KIND=8) :: bottom_snd_buffer(:),bottom_rcv_buffer(:),top_snd_buffer(:),top_rcv_buffer(:)

  INTEGER      :: chunk,depth,field_type

  INTEGER      :: size,err,request(8),tag,message_count,j,k,x_inc,y_inc,index
  INTEGER      :: receiver,sender

  ! Field type will either be cell, vertex, x_face or y_face to get the message limits correct

  !  but this will do for now

  ! I am also sending buffers to chunks with the same task id for now.
  ! This can be improved in the future but at the moment there is just 1 chunk per task anyway

  ! The tag will be a function of the sending chunk and the face it is coming from
  !  like chunk 6 sending the left face

  ! No open mp in here either. May be beneficial will packing and unpacking in the future, though I am not sure.

  ! Change this so it will allow more than 1 chunk per task

  request=0
  message_count=0

  ! Pack and send

  ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
  IF(field_type.EQ.CELL_DATA) THEN
    x_inc=0
    y_inc=0
  ENDIF
  IF(field_type.EQ.VERTEX_DATA) THEN
    x_inc=1
    y_inc=1
  ENDIF
  IF(field_type.EQ.X_FACE_DATA) THEN
    x_inc=1
    y_inc=0
  ENDIF
  IF(field_type.EQ.Y_FACE_DATA) THEN
    x_inc=0
    y_inc=1
  ENDIF

  ! Pack real data into buffers
  IF(parallel%task.EQ.chunks(chunk)%task) THEN
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      size=(1+(chunks(chunk)%field%y_max+y_inc+depth)-(chunks(chunk)%field%y_min-depth))*depth
      DO k=chunks(chunk)%field%y_min-depth,chunks(chunk)%field%y_max+y_inc+depth
        DO j=1,depth
          index=j+(k+depth-1)*depth
          left_snd_buffer(index)=field(chunks(chunk)%field%x_min+x_inc-1+j,k)
        ENDDO
      ENDDO
      tag=4*(chunk)+1 ! 4 because we have 4 faces, 1 because it is leaving the left face
      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task
      tag=4*(chunks(chunk)%chunk_neighbours(chunk_left))+2 ! 4 because we have 4 faces, 1 because it is coming from the right face of the left neighbour
      sender=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task
      message_count=message_count+2
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      size=(1+(chunks(chunk)%field%y_max+y_inc+depth)-(chunks(chunk)%field%y_min-depth))*depth
      DO k=chunks(chunk)%field%y_min-depth,chunks(chunk)%field%y_max+y_inc+depth
        DO j=1,depth
          index=j+(k+depth-1)*depth
          right_snd_buffer(index)=field(chunks(chunk)%field%x_max+1-j,k)
        ENDDO
      ENDDO
      tag=4*chunk+2 ! 4 because we have 4 faces, 2 because it is leaving the right face
      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task
      tag=4*(chunks(chunk)%chunk_neighbours(chunk_right))+1 ! 4 because we have 4 faces, 1 because it is coming from the left face of the right neighbour
      sender=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task
      message_count=message_count+2
    ENDIF
  ENDIF

  ! Wait for the messages


  ! Unpack buffers in halo cells
  IF(parallel%task.EQ.chunks(chunk)%task) THEN
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      DO k=chunks(chunk)%field%y_min-depth,chunks(chunk)%field%y_max+y_inc+depth
        DO j=1,depth
          index=j+(k+depth-1)*depth
          field(chunks(chunk)%field%x_min-j,k)=left_rcv_buffer(index)
        ENDDO
      ENDDO
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      DO k=chunks(chunk)%field%y_min-depth,chunks(chunk)%field%y_max+y_inc+depth
        DO j=1,depth
          index=j+(k+depth-1)*depth
          field(chunks(chunk)%field%x_max+x_inc+j,k)=right_rcv_buffer(index)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  request=0
  message_count=0

  IF(parallel%task.EQ.chunks(chunk)%task) THEN
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      size=(1+(chunks(chunk)%field%x_max+x_inc+depth)-(chunks(chunk)%field%x_min-depth))*depth
      DO k=1,depth
        DO j=chunks(chunk)%field%x_min-depth,chunks(chunk)%field%x_max+x_inc+depth
          index=j+depth+(k-1)*(chunks(chunk)%field%x_max+x_inc+(2*depth))
          bottom_snd_buffer(index)=field(j,chunks(chunk)%field%y_min+y_inc-1+k)
        ENDDO
      ENDDO
      tag=4*(chunk)+3 ! 4 because we have 4 faces, 3 because it is leaving the bottom face
      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task
      tag=4*(chunks(chunk)%chunk_neighbours(chunk_bottom))+4 ! 4 because we have 4 faces, 1 because it is coming from the top face of the bottom neighbour
      sender=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task
      message_count=message_count+2
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      size=(1+(chunks(chunk)%field%x_max+x_inc+depth)-(chunks(chunk)%field%x_min-depth))*depth
      DO k=1,depth
        DO j=chunks(chunk)%field%x_min-depth,chunks(chunk)%field%x_max+x_inc+depth
          index=j+depth+(k-1)*(chunks(chunk)%field%x_max+x_inc+(2*depth))
          top_snd_buffer(index)=field(j,chunks(chunk)%field%y_max+1-k)
        ENDDO
      ENDDO
      tag=4*(chunk)+4 ! 4 because we have 4 faces, 4 because it is leaving the top face
      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task
      tag=4*(chunks(chunk)%chunk_neighbours(chunk_top))+3 ! 4 because we have 4 faces, 4 because it is coming from the left face of the top neighbour
      sender=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task
      message_count=message_count+2
    ENDIF
  ENDIF

  ! Wait for the messages

  ! Unpack buffers in halo cells
  IF(parallel%task.EQ.chunks(chunk)%task) THEN
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      DO k=1,depth
        DO j=chunks(chunk)%field%x_min-depth,chunks(chunk)%field%x_max+x_inc+depth
          index=j+depth+(k-1)*(chunks(chunk)%field%x_max+x_inc+(2*depth))
          field(j,chunks(chunk)%field%y_min-k)=bottom_rcv_buffer(index)
        ENDDO
      ENDDO
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      DO k=1,depth
        DO j=chunks(chunk)%field%x_min-depth,chunks(chunk)%field%x_max+x_inc+depth
          index=j+depth+(k-1)*(chunks(chunk)%field%x_max+x_inc+(2*depth))
          field(j,chunks(chunk)%field%y_max+y_inc+k)=top_rcv_buffer(index)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

END SUBROUTINE clover_exchange_message

SUBROUTINE clover_sum(value)

  ! Only sums to the master

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  INTEGER :: err

  total=value


  value=total

END SUBROUTINE clover_sum

SUBROUTINE clover_min(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: minimum

  INTEGER :: err

  minimum=value


  value=minimum

END SUBROUTINE clover_min

SUBROUTINE clover_check_error(error)

  IMPLICIT NONE

  INTEGER :: error

  INTEGER :: maximum

  INTEGER :: err

  maximum=error


  error=maximum

END SUBROUTINE clover_check_error


END MODULE clover_module
