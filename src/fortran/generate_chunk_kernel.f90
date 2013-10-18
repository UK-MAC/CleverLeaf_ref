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

!>  @brief Fortran mesh chunk generator
!>  @author Wayne Gaudin
!>  @details Generates the field data on a mesh chunk based on the user specified
!>  input for the states.
!>
!>  Note that state one is always used as the background state, which is then
!>  overwritten by further state definitions.

SUBROUTINE generate_chunk_kernel(x_min,x_max,y_min,y_max, &
                                 vertexx,                 &
                                 vertexy,                 &
                                 cellx,                   &
                                 celly,                   &
                                 density0,                &
                                 energy0,                 &
                                 xvel0,                   &
                                 yvel0,                   &
                                 number_of_states,        &
                                 state_density,           &
                                 state_energy,            &
                                 state_xvel,              &
                                 state_yvel,              &
                                 state_xmin,              &
                                 state_xmax,              &
                                 state_ymin,              &
                                 state_ymax,              &
                                 state_radius,            &
                                 state_geometry,          &
                                 g_rect,                  &
                                 g_circ,                  &
                                 g_point)

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3, y_min-2:y_max+3) :: vertexx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3, y_min-2:y_max+3) :: vertexy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2, y_min-2:y_max+2) :: cellx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2, y_min-2:y_max+2) :: celly
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
  INTEGER      :: number_of_states
  REAL(KIND=8), DIMENSION(1:number_of_states) :: state_density
  REAL(KIND=8), DIMENSION(1:number_of_states) :: state_energy
  REAL(KIND=8), DIMENSION(1:number_of_states) :: state_xvel
  REAL(KIND=8), DIMENSION(1:number_of_states) :: state_yvel
  REAL(KIND=8), DIMENSION(1:number_of_states) :: state_xmin
  REAL(KIND=8), DIMENSION(1:number_of_states) :: state_xmax
  REAL(KIND=8), DIMENSION(1:number_of_states) :: state_ymin
  REAL(KIND=8), DIMENSION(1:number_of_states) :: state_ymax
  REAL(KIND=8), DIMENSION(1:number_of_states) :: state_radius
  INTEGER     , DIMENSION(1:number_of_states) :: state_geometry
  INTEGER      :: g_rect
  INTEGER      :: g_circ
  INTEGER      :: g_point

  REAL(KIND=8) :: radius,x_cent,y_cent
  INTEGER      :: state

  INTEGER      :: j,k,jt,kt

  ! State 1 is always the background state

!$OMP PARALLEL SHARED(x_cent,y_cent)
!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
      energy0(j,k)=state_energy(1)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
      density0(j,k)=state_density(1)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
      xvel0(j,k)=state_xvel(1)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
      yvel0(j,k)=state_yvel(1)
    ENDDO
  ENDDO
!$OMP END DO

  DO state=2,number_of_states

! Could the velocity setting be thread unsafe?
    x_cent=state_xmin(state)
    y_cent=state_ymin(state)

!$OMP DO PRIVATE(radius,jt,kt)
    DO k=y_min-2,y_max+2
      DO j=x_min-2,x_max+2
        IF(state_geometry(state).EQ.g_rect ) THEN
          IF( (vertexx(j+1,k) - state_xmin(state) .GT. (ABS(vertexx(j+1,k))*EPSILON(1.0_8)) &
             .OR. (ABS(vertexx(j+1,k) - state_xmin(state))) .LE. (ABS(vertexx(j+1,k))*EPSILON(1.0_8))) &
              ! check whether vertexx is less than state_xmax, taking into account machine epsilon
              .AND. (state_xmax(state) - vertexx(j,k)) .GT. (abs(vertexx(j,k))*epsilon(1.0_8)) ) THEN
            IF( (vertexy(j,k+1) - state_ymin(state) .GT. (ABS(vertexy(j,k+1))*EPSILON(1.0_8)) &
                .OR. (ABS(vertexy(j,k+1) - state_ymin(state))) .LE.  (ABS(vertexy(j,k+1))*EPSILON(1.0_8))) &
              ! check whether vertexy is less than state_xmax, taking into account machine epsilon
                .AND. (state_ymax(state) - vertexy(j,k)) .GT. (abs(vertexy(j,k))*epsilon(1.0_8)) )  THEN
              energy0(j,k)=state_energy(state)
              density0(j,k)=state_density(state)
              DO kt=k,k+1
                DO jt=j,j+1
                  xvel0(jt,kt)=state_xvel(state)
                  yvel0(jt,kt)=state_yvel(state)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ELSEIF(state_geometry(state).EQ.g_circ ) THEN
          radius=SQRT((cellx(j,k)-x_cent)*(cellx(j,k)-x_cent)+(celly(j,k)-y_cent)*(celly(j,k)-y_cent))
          IF( (state_radius(state) - radius) .GT. (ABS(radius)*EPSILON(1.0_8)) &
              .OR. ABS(radius - state_radius(state)) .LE. (ABS(radius)*EPSILON(1.0_8))) THEN
            energy0(j,k)=state_energy(state)
            density0(j,k)=state_density(state)
            DO kt=k,k+1
              DO jt=j,j+1
                xvel0(jt,kt)=state_xvel(state)
                yvel0(jt,kt)=state_yvel(state)
              ENDDO
            ENDDO
          ENDIF
        ELSEIF(state_geometry(state).EQ.g_point) THEN
            ! check whether vertexx is approximately equal to x_cent and vertexy
            ! is approximately equal to y_cent
          IF ( ABS(vertexx(j,k) - x_cent) .LE. (ABS(vertexx(j,k))*EPSILON(1.0_8)) &
              .AND. ABS(vertexy(j,k) - y_cent) .LE. (ABS(vertexy(j,k))*EPSILON(1.0_8))) THEN
            energy0(j,k)=state_energy(state)
            density0(j,k)=state_density(state)
            DO kt=k,k+1
              DO jt=j,j+1
                xvel0(jt,kt)=state_xvel(state)
                yvel0(jt,kt)=state_yvel(state)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
    ENDDO
!$OMP END DO

  ENDDO

!$OMP END PARALLEL

END SUBROUTINE generate_chunk_kernel
