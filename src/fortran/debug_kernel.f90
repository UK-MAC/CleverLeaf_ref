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

!>  @brief Fortran ideal gas kernel.
!>  @author Wayne Gaudin
!>  @details Calculates the pressure and sound speed for the mesh chunk using
!>  the ideal gas equation of state, with a fixed gamma of 1.4.

SUBROUTINE debug_kernel(x_min,x_max,y_min,y_max,                &
                            xvel0,                                &
                            yvel0,                                 &
                            xvel1,                               &
                            yvel1                              )                              

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: yvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: yvel1

  INTEGER :: j,k

  PRINT *, 'xvel1: '
  DO j=x_min,x_max+1
    PRINT *, xvel1(j,y_min+2)
  ENDDO

END SUBROUTINE debug_kernel
