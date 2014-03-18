! 
! Crown Copyright 2012 AWE, Copyright 2013 David Beckingsale.
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

!>  @brief Fortran cell advection kernel.
!>  @author Wayne Gaudin
!>  @details Performs a second order advective remap using van-Leer limiting
!>  with directional splitting.
SUBROUTINE advec_cell_kernel(x_min,       &
                             x_max,       &
                             y_min,       &
                             y_max,       &
                             dir,         &
                             sweep_number,&
                             vertexdx,    &
                             vertexdy,    &
                             volume,      &
                             density1,    &
                             energy1,     &
                             mass_flux_x, &
                             vol_flux_x,  &
                             mass_flux_y, &
                             vol_flux_y,  &
                             pre_vol,     &
                             post_vol,    &
                             pre_mass,    &
                             post_mass,   &
                             advec_vol,   &
                             post_ener,   &
                             ener_flux    )

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  INTEGER :: sweep_number,dir
  INTEGER :: g_xdir=1,g_ydir=2

  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy1
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

  REAL(KIND=8), DIMENSION(x_min-2:x_max+3, y_min-2:y_max+3) :: vertexdx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3, y_min-2:y_max+3) :: vertexdy

  INTEGER :: j,k,upwind,donor,downwind,dif

  REAL(KIND=8) :: sigma,sigmat,sigmav,sigmam,sigma3,sigma4
  REAL(KIND=8) :: diffuw,diffdw,limiter
  REAL(KIND=8), PARAMETER :: one_by_six=1.0_8/6.0_8

!$OMP PARALLEL

  IF(dir.EQ.g_xdir) THEN

    IF(sweep_number.EQ.1)THEN
!$OMP DO
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          pre_vol(j,k)=volume(j,k)+(vol_flux_x(j+1,k  )-vol_flux_x(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k))
          post_vol(j,k)=pre_vol(j,k)-(vol_flux_x(j+1,k  )-vol_flux_x(j,k))
        ENDDO
      ENDDO 
!$OMP END DO
    ELSE
!$OMP DO
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          pre_vol(j,k)=volume(j,k)+vol_flux_x(j+1,k)-vol_flux_x(j,k)
          post_vol(j,k)=volume(j,k)
        ENDDO
      ENDDO 
!$OMP END DO
    ENDIF

!$OMP DO PRIVATE(upwind,donor,downwind,dif,sigmat,sigma3,sigma4,sigmav,sigma,sigmam, &
!$OMP            diffuw,diffdw,limiter)
    DO k=y_min,y_max
      DO j=x_min,x_max+2

        IF(vol_flux_x(j,k).GT.0.0)THEN
          upwind   =j-2
          donor    =j-1
          downwind =j
          dif      =donor
        ELSE
          upwind   =MIN(j+1,x_max+2)
          donor    =j
          downwind =j-1
          dif      =upwind
        ENDIF

        sigmat=ABS(vol_flux_x(j,k))/pre_vol(donor,k)
        sigma3=(1.0_8+sigmat)*(vertexdx(j,k)/vertexdx(dif,k))
        sigma4=2.0_8-sigmat

        sigma=sigmat
        sigmav=sigmat

        diffuw=density1(donor,k)-density1(upwind,k)
        diffdw=density1(downwind,k)-density1(donor,k)
        IF(diffuw*diffdw.GT.0.0)THEN
          limiter=(1.0_8-sigmav)*SIGN(1.0_8,diffdw)*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
        ELSE
          limiter=0.0
        ENDIF

        IF (j .EQ. 483 .AND. k .EQ. 3) THEN
          WRITE(*,*) 'donor = ', donor
          WRITE(*,*) 'downwind = ', downwind
          WRITE(*,*) 'upwind = ', upwind
          WRITE(*,*) 'dif = ', dif
          WRITE(*,*) 'limiter = ', limiter
          WRITE(*,*) 'sigmat = ', sigmat
          WRITE(*,*) 'sigma3 = ', sigma3
          WRITE(*,*) 'sigma4 = ', sigma4
          WRITE(*,*) 'sigma = ', sigma
          WRITE(*,*) 'sigmav = ', sigmav
          WRITE(*,*) 'diffuw = ', diffuw
          WRITE(*,*) 'diffdw = ', diffdw
          WRITE(*,*) 'density1(donor,k) = ', density1(donor,k)
          WRITE(*,*) 'density1(upwind,k) = ', density1(upwind,k)
          WRITE(*,*) 'density1(downwind,k) = ', density1(downwind,k)
          WRITE(*,*) 'vol_flux_x(j,k)', vol_flux_x(j,k)
          WRITE(*,*) 'vol_flux_x(j+1,k)', vol_flux_x(j+1,k)
          WRITE(*,*) 'pre_vol', pre_vol(j,k)
          WRITE(*,*) 'pre_vol(donor,k)', pre_vol(donor,k)
          WRITE(*,*) 'vol_flux_x(donor,k)', vol_flux_x(donor,k)
          WRITE(*,*) 'vol_flux_x(donor+1,k)', vol_flux_x(donor+1,k)
          WRITE(*,*) 'volume(donor,k)', volume(donor,k)
        ENDIF

        mass_flux_x(j,k)=vol_flux_x(j,k)*(density1(donor,k)+limiter)

        sigmam=ABS(mass_flux_x(j,k))/(density1(donor,k)*pre_vol(donor,k))
        diffuw=energy1(donor,k)-energy1(upwind,k)
        diffdw=energy1(downwind,k)-energy1(donor,k)
        IF(diffuw*diffdw.GT.0.0)THEN
          limiter=(1.0_8-sigmam)*SIGN(1.0_8,diffdw)*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
        ELSE
          limiter=0.0
        ENDIF

        ener_flux(j,k)=mass_flux_x(j,k)*(energy1(donor,k)+limiter)

      ENDDO
    ENDDO
!$OMP END DO

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        pre_mass(j,k)=density1(j,k)*pre_vol(j,k)
        post_mass(j,k)=pre_mass(j,k)+(mass_flux_x(j,k)-mass_flux_x(j+1,k))
        post_ener(j,k)=(energy1(j,k)*pre_mass(j,k)+ener_flux(j,k)-ener_flux(j+1,k))/post_mass(j,k)
        advec_vol(j,k)=pre_vol(j,k)+vol_flux_x(j,k)-vol_flux_x(j+1,k)
        density1(j,k)=post_mass(j,k)/advec_vol(j,k)
        energy1(j,k)=post_ener(j,k)
      IF (j .EQ. 483 .AND. k .EQ. 3 .AND. density1(j,k) .LT. 0.0_8) THEN
        WRITE(*,*) 'advec'
        WRITE(*,*) 'pre_mass', pre_mass(j,k)
        WRITE(*,*) 'pre_vol', pre_vol(j,k)
        WRITE(*,*) 'pre_vol(donor,k)', pre_vol(482,k)
        WRITE(*,*) 'mass_flux_x(j,k)', mass_flux_x(j,k)
        WRITE(*,*) 'mass_flux_x(j+1,k)', mass_flux_x(j+1,k)
        WRITE(*,*) 'vol_flux_x(j,k)', vol_flux_x(j,k)
        WRITE(*,*) 'vol_flux_x(j+1,k)', vol_flux_x(j+1,k)
        WRITE(*,*) 'density1(donor,k)', density1(482,k)
        WRITE(*,*) 'density1(upwind,k)', density1(481,k)
        WRITE(*,*) 'density1(downwind,k)', density1(483,k)
        WRITE(*,*) 'post_mass', post_mass(j,k)
        WRITE(*,*) 'post_ener', post_ener(j,k)
        WRITE(*,*) 'advec_vol', advec_vol(j,k)
        WRITE(*,*) 'density1', density1(j,k)
        WRITE(*,*) 'energy1', energy1(j,k)
        WRITE(*,*) 'vertexdx', vertexdx(j,k)
        WRITE(*,*) 'vertexdx', vertexdx(482,k)
      ENDIF
      ENDDO
    ENDDO
!$OMP END DO

  ELSEIF(dir.EQ.g_ydir) THEN

    IF(sweep_number.EQ.1)THEN
!$OMP DO
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          pre_vol(j,k)=volume(j,k)+(vol_flux_y(j  ,k+1)-vol_flux_y(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k))
          post_vol(j,k)=pre_vol(j,k)-(vol_flux_y(j  ,k+1)-vol_flux_y(j,k))
        ENDDO
      ENDDO
!$OMP END DO
    ELSE
!$OMP DO
      DO k=y_min-2,y_max+2
        DO j=x_min-2,x_max+2
          pre_vol(j,k)=volume(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k)
          post_vol(j,k)=volume(j,k)
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF

!$OMP DO PRIVATE(upwind,donor,downwind,dif,sigmat,sigma3,sigma4,sigmav,sigma,sigmam, &
!$OMP            diffuw,diffdw,limiter)
    DO k=y_min,y_max+2
      DO j=x_min,x_max

        IF(vol_flux_y(j,k).GT.0.0)THEN
          upwind   =k-2
          donor    =k-1
          downwind =k
          dif      =donor
        ELSE
          upwind   =MIN(k+1,y_max+2)
          donor    =k
          downwind =k-1
          dif      =upwind
        ENDIF

        sigmat=ABS(vol_flux_y(j,k))/pre_vol(j,donor)
        sigma3=(1.0_8+sigmat)*(vertexdy(j,k)/vertexdy(j,dif))
        sigma4=2.0_8-sigmat

        sigma=sigmat
        sigmav=sigmat

        diffuw=density1(j,donor)-density1(j,upwind)
        diffdw=density1(j,downwind)-density1(j,donor)
        IF(diffuw*diffdw.GT.0.0)THEN
          limiter=(1.0_8-sigmav)*SIGN(1.0_8,diffdw)*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
        ELSE
          limiter=0.0
        ENDIF
        mass_flux_y(j,k)=vol_flux_y(j,k)*(density1(j,donor)+limiter)

        sigmam=ABS(mass_flux_y(j,k))/(density1(j,donor)*pre_vol(j,donor))
        diffuw=energy1(j,donor)-energy1(j,upwind)
        diffdw=energy1(j,downwind)-energy1(j,donor)
        IF(diffuw*diffdw.GT.0.0)THEN
          limiter=(1.0_8-sigmam)*SIGN(1.0_8,diffdw)*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
        ELSE
          limiter=0.0
        ENDIF
        ener_flux(j,k)=mass_flux_y(j,k)*(energy1(j,donor)+limiter)

      ENDDO
    ENDDO
!$OMP END DO

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        pre_mass(j,k)=density1(j,k)*pre_vol(j,k)
        post_mass(j,k)=pre_mass(j,k)+mass_flux_y(j,k)-mass_flux_y(j,k+1)
        post_ener(j,k)=(energy1(j,k)*pre_mass(j,k)+ener_flux(j,k)-ener_flux(j,k+1))/post_mass(j,k)
        advec_vol(j,k)=pre_vol(j,k)+vol_flux_y(j,k)-vol_flux_y(j,k+1)
        density1(j,k)=post_mass(j,k)/advec_vol(j,k)
        energy1(j,k)=post_ener(j,k)
      IF (j .EQ. 483 .AND. k .EQ. 3 .AND. density1(j,k) .LT. 0.0) THEN
        WRITE(*,*) 'advec'
      ENDIF
      ENDDO
    ENDDO
!$OMP END DO

  ENDIF

!$OMP END PARALLEL

END SUBROUTINE advec_cell_kernel
