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

!>  @brief Fortran mpi buffer packing kernel
!>  @author Wayne Gaudin
!>  @details Packs/unpacks mpi send and receive buffers

MODULE pack_kernel_module

CONTAINS

SUBROUTINE pack_left_right_buffers(x_min,x_max,y_min,y_max,              &
                                   chunk_left,chunk_right,external_face, &
                                   x_inc,y_inc,depth,size,               &
                                   field,left_snd_buffer,right_snd_buffer)

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  INTEGER      :: chunk_left,chunk_right,external_face
  INTEGER      :: x_inc,y_inc,depth,size

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  REAL(KIND=8) :: left_snd_buffer(:),right_snd_buffer(:)

  INTEGER      :: j,k,index

  IF(chunk_left.NE.external_face) THEN
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=j+(k+depth-1)*depth
        left_snd_buffer(index)=field(x_min+x_inc-1+j,k)
      ENDDO
    ENDDO
  ENDIF
  IF(chunk_right.NE.external_face) THEN
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=j+(k+depth-1)*depth
        right_snd_buffer(index)=field(x_max+1-j,k)
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE pack_left_right_buffers

SUBROUTINE unpack_left_right_buffers(x_min,x_max,y_min,y_max,              &
                                     chunk_left,chunk_right,external_face, &
                                     x_inc,y_inc,depth,size,               &
                                     field,left_rcv_buffer,right_rcv_buffer)

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  INTEGER      :: chunk_left,chunk_right,external_face
  INTEGER      :: x_inc,y_inc,depth,size

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  REAL(KIND=8) :: left_rcv_buffer(:),right_rcv_buffer(:)

  INTEGER      :: j,k,index

  IF(chunk_left.NE.external_face) THEN
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=j+(k+depth-1)*depth
        field(x_min-j,k)=left_rcv_buffer(index)
      ENDDO
    ENDDO
  ENDIF
  IF(chunk_right.NE.external_face) THEN
    DO k=y_min-depth,y_max+y_inc+depth
      DO j=1,depth
        index=j+(k+depth-1)*depth
        field(x_max+x_inc+j,k)=right_rcv_buffer(index)
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE unpack_left_right_buffers

SUBROUTINE pack_top_bottom_buffers(x_min,x_max,y_min,y_max,              &
                                   chunk_bottom,chunk_top,external_face, &
                                   x_inc,y_inc,depth,size,               &
                                   field,bottom_snd_buffer,top_snd_buffer)

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  INTEGER      :: chunk_bottom,chunk_top,external_face
  INTEGER      :: x_inc,y_inc,depth,size

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  REAL(KIND=8) :: bottom_snd_buffer(:),top_snd_buffer(:)

  INTEGER      :: j,k,index

  IF(chunk_bottom.NE.external_face) THEN
    DO k=1,depth
      DO j=x_min-depth,x_max+x_inc+depth
        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
        bottom_snd_buffer(index)=field(j,y_min+y_inc-1+k)
      ENDDO
    ENDDO
  ENDIF
  IF(chunk_top.NE.external_face) THEN
    DO k=1,depth
      DO j=x_min-depth,x_max+x_inc+depth
        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
        top_snd_buffer(index)=field(j,y_max+1-k)
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE pack_top_bottom_buffers

SUBROUTINE unpack_top_bottom_buffers(x_min,x_max,y_min,y_max,             &
                                    chunk_bottom,chunk_top,external_face, &
                                    x_inc,y_inc,depth,size,               &
                                    field,bottom_rcv_buffer,top_rcv_buffer)

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  INTEGER      :: chunk_bottom,chunk_top,external_face
  INTEGER      :: x_inc,y_inc,depth,size

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  REAL(KIND=8) :: bottom_rcv_buffer(:),top_rcv_buffer(:)

  INTEGER      :: j,k,index

  IF(chunk_bottom.NE.external_face) THEN
    DO k=1,depth
      DO j=x_min-depth,x_max+x_inc+depth
        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
        field(j,y_min-k)=bottom_rcv_buffer(index)
      ENDDO
    ENDDO
  ENDIF
  IF(chunk_top.NE.external_face) THEN
    DO k=1,depth
      DO j=x_min-depth,x_max+x_inc+depth
        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
        field(j,y_max+y_inc+k)=top_rcv_buffer(index)
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE unpack_top_bottom_buffers

END MODULE pack_kernel_module
