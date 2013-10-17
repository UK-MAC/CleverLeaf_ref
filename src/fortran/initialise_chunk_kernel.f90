! Crown Copyright 2012 AWE, Copyright 2013 David Beckingsale
!
! This file is part of CleverLeaf.
!
! CleverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CleverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CleverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran chunk initialisation kernel.
!>  @author Wayne Gaudin, David Beckingsale.
!>  @details Calculates mesh geometry for the mesh chunk based on the mesh size.

SUBROUTINE initialise_chunk_kernel(x_min,x_max,y_min,y_max,       &
                                   xmin,ymin,dx,dy,               &
                                   vertexx,                       &
                                   vertexdx,                      &
                                   vertexy,                       &
                                   vertexdy,                      &
                                   cellx,                         &
                                   celldx,                        &
                                   celly,                         &
                                   celldy,                        &
                                   volume)

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  REAL(KIND=8) :: xmin,ymin,dx,dy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexdx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexy
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexdy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: cellx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celly
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celldy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: volume

  INTEGER      :: j,k

!$OMP PARALLEL
!$OMP DO
  DO j=x_min-2,x_max+3
     vertexx(j)=xmin+dx*float(j-x_min)
  ENDDO
!$OMP END DO

!$OMP DO
  DO j=x_min-2,x_max+3
    vertexdx(j)=dx
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+3
     vertexy(k)=ymin+dy*float(k-y_min)
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+3
    vertexdy(k)=dy
  ENDDO
!$OMP END DO

!$OMP DO
  DO j=x_min-2,x_max+2
     cellx(j)=0.5*(vertexx(j)+vertexx(j+1))
  ENDDO
!$OMP END DO

!$OMP DO
  DO j=x_min-2,x_max+2
     celldx(j)=dx
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+2
     celly(k)=0.5*(vertexy(k)+vertexy(k+1))
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+2
     celldy(k)=dy
  ENDDO
!$OMP END DO

!$OMP DO PRIVATE(j)
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
        volume(j,k)=dx*dy
     ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE initialise_chunk_kernel
