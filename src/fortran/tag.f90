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

!>  @brief Fortran viscosity tagging kernel
!>  @author David Beckingsale
!>  @details Tags cells for refinement based on the absolute viscosity value.
SUBROUTINE tag_q_kernel(x_min,x_max,y_min,y_max,&
                            q_threshold,        &
                            viscosity,          &
                            tags)

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: viscosity
  INTEGER, DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: tags
  REAL(KIND=8) :: q_threshold

  INTEGER :: j,k

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      IF (viscosity(j,k) .GT. q_threshold) THEN
        tags(j,k) = 1
        tags(j+1,k) = 1
        tags(j,k+1) = 1
        tags(j-1,k) = 1
        tags(j,k-1) = 1
        tags(j-1,k-1) = 1
        tags(j+1,k-1) = 1
        tags(j-1,k+1) = 1
        tags(j+1,k+1) = 1
      ENDIF
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tag_q_kernel

!>  @brief Fortran density tagging kernel
!>  @author David Beckingsale
!>  @details Tags cells for refinement based on the density gradient.
SUBROUTINE tag_density_kernel(x_min,x_max,y_min,y_max,&
                            density_gradient_threshold,&
                            density,           &
                            tags)

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density
  INTEGER, DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: tags
  REAL(KIND=8) :: density_gradient_threshold

  INTEGER :: j,k
  REAL(KIND=8) :: d2x,d2y,dxy,dyx,dd

!$OMP PARALLEL
!$OMP DO PRIVATE(d2x,d2y,dxy,dyx,dd)
  DO k=y_min,y_max
    DO j=x_min,x_max
      d2x = abs(density(j+1,k) - 2.0*density(j,k) + density(j-1,k));
      d2y = abs(density(j,k+1) - 2.0*density(j,k) + density(j,k-1));

      dxy = abs(density(j+1,k+1) - 2.0*density(j,k) + density(j-1,k-1));
      dyx = abs(density(j-1,k+1) - 2.0*density(j,k) + density(j+1,k-1));

      dd = max(d2x,max(d2y,max(dxy,dyx)));

      IF (dd .GT. density_gradient_threshold) THEN
        tags(j,k) = 1
        tags(j+1,k) = 1
        tags(j,k+1) = 1
        tags(j-1,k) = 1
        tags(j,k-1) = 1
        tags(j-1,k-1) = 1
        tags(j+1,k-1) = 1
        tags(j-1,k+1) = 1
        tags(j+1,k+1) = 1
      ENDIF
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tag_density_kernel

!>  @brief Fortran energy tagging kernel
!>  @author David Beckingsale
!>  @details Tags cells for refinement based on the energy gradient.
SUBROUTINE tag_energy_kernel(x_min,x_max,y_min,y_max,&
                            energy_gradient_threshold,&
                            energy,           &
                            tags)

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
  INTEGER, DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: tags
  REAL(KIND=8) :: energy_gradient_threshold

  INTEGER :: j,k
  REAL(KIND=8) :: d2x,d2y,dxy,dyx,dd


!$OMP PARALLEL
!$OMP DO PRIVATE(d2x,d2y,dxy,dyx,dd)
  DO k=y_min,y_max
    DO j=x_min,x_max
      d2x = abs(energy(j+1,k) - 2.0*energy(j,k) + energy(j-1,k));
      d2y = abs(energy(j,k+1) - 2.0*energy(j,k) + energy(j,k-1));

      dxy = abs(energy(j+1,k+1) - 2.0*energy(j,k) + energy(j-1,k-1));
      dyx = abs(energy(j-1,k+1) - 2.0*energy(j,k) + energy(j+1,k-1));

      dd = max(d2x,max(d2y,max(dxy,dyx)));

      IF (dd .GT. energy_gradient_threshold) THEN
        tags(j,k) = 1
        tags(j+1,k) = 1
        tags(j,k+1) = 1
        tags(j-1,k) = 1
        tags(j,k-1) = 1
        tags(j-1,k-1) = 1
        tags(j+1,k-1) = 1
        tags(j-1,k+1) = 1
        tags(j+1,k+1) = 1
      ENDIF
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tag_energy_kernel

!>  @brief Fortran "tag all" kernel
!>  @author David Beckingsale
!>  @details Tags all cells for refinement - useful for debugging.
SUBROUTINE tag_all_kernel(x_min,x_max,y_min,y_max,&
                            tags)

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  INTEGER, DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: tags

  INTEGER :: j,k

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
        tags(j,k) = 1
        tags(j+1,k) = 1
        tags(j,k+1) = 1
        tags(j-1,k) = 1
        tags(j,k-1) = 1
        tags(j-1,k-1) = 1
        tags(j+1,k-1) = 1
        tags(j-1,k+1) = 1
        tags(j+1,k+1) = 1
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tag_all_kernel
