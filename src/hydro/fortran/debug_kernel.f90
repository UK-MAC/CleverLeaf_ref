! 
! Copyright 2013 David Beckingsale.
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
!  
!>  @brief Fortran debug kernel.
!>  @author David Beckingsale
!>  @details This kernel is used for debugging, providing a convenient
!>           breakpoint for inspecting patch variables.
SUBROUTINE debug_kernel(x_min,x_max,y_min,y_max,                &
                            density0,                           &
                            density1,                           &
                            energy0,                            &
                            energy1,                            &
                            pressure,                           &
                            soundspeed,                         &
                            viscosity,                          &
                            xvel0,                              &
                            yvel0,                              &
                            xvel1,                              &
                            yvel1,                              &
                            vol_flux_x,                         &
                            vol_flux_y,                         &
                            mass_flux_x,                        &
                            mass_flux_y,                        &
                            pre_vol,                            &
                            post_vol,                           &
                            pre_mass,                           &
                            post_mass,                          &
                            advec_vol,                          &
                            post_ener,                          &
                            ener_flux)

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: pressure
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: soundspeed
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: viscosity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: yvel0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: yvel1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: mass_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: mass_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_vol
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_vol
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_mass
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_mass
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: advec_vol
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_ener
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: ener_flux

  INTEGER :: j,k

  REAL(KIND=8) :: tmp

  WRITE(*, *) 'density0'
  DO k=y_min-2, y_max+2
    DO j=x_min-2, x_max+2
      WRITE(*, '(7E12.4,X)', advance='no') density0(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'density1'
  DO k=y_min-2, y_max+2
    DO j=x_min-2, x_max+2
      WRITE(*, '(7E12.4,X)', advance='no') density1(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'energy0'
  DO k=y_min-2, y_max+2
    DO j=x_min-2, x_max+2
      WRITE(*, '(7E12.4,X)', advance='no') energy0(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'energy1'
  DO k=y_min-2, y_max+2
    DO j=x_min-2, x_max+2
      WRITE(*, '(7E12.4,X)', advance='no') energy1(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'pressure'
  DO k=y_min-2, y_max+2
    DO j=x_min-2, x_max+2
      WRITE(*, '(7E12.4,X)', advance='no') pressure(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'soundspeed'
  DO k=y_min-2, y_max+2
    DO j=x_min-2, x_max+2
      WRITE(*, '(7E12.4,X)', advance='no') soundspeed(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'viscosity'
  DO k=y_min-2, y_max+2
    DO j=x_min-2, x_max+2
      WRITE(*, '(7E12.4,X)', advance='no') viscosity(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'xvel0'
  DO k=y_min-2, y_max+3
    DO j=x_min-2, x_max+3
      WRITE(*, '(7E12.4,X)', advance='no') xvel0(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'yvel0'
  DO k=y_min-2, y_max+3
    DO j=x_min-2, x_max+3
      WRITE(*, '(7E12.4,X)', advance='no') yvel0(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'xvel1'
  DO k=y_min-2, y_max+3
    DO j=x_min-2, x_max+3
      WRITE(*, '(7E12.4,X)', advance='no') xvel1(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'yvel1'
  DO k=y_min-2, y_max+3
    DO j=x_min-2, x_max+3
      WRITE(*, '(7E12.4,X)', advance='no') yvel1(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'vol_flux_x'
  DO k=y_min-2, y_max+2
    DO j=x_min-2, x_max+3
      WRITE(*, '(7E12.4,X)', advance='no') vol_flux_x(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'vol_flux_y'
  DO k=y_min-2, y_max+3
    DO j=x_min-2, x_max+2
      WRITE(*, '(7E12.4,X)', advance='no') vol_flux_y(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'mass_flux_x'
  DO k=y_min-2, y_max+2
    DO j=x_min-2, x_max+3
      WRITE(*, '(7E12.4,X)', advance='no') mass_flux_x(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO
  WRITE(*, *) 'mass_flux_y'
  DO k=y_min-2, y_max+3
    DO j=x_min-2, x_max+2
      WRITE(*, '(7E12.4,X)', advance='no') mass_flux_y(j,k)
    ENDDO
    WRITE(*, *) ''  ! this gives you the line break
  ENDDO

END SUBROUTINE debug_kernel
