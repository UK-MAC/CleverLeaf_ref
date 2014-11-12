!
! Crown Copyright 2014 AWE, Copyright 2014 David Beckingsale.
!
! This file is part of CleverLeaf.
!
! CleverLeaf is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! CleverLeaf is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! CleverLeaf. If not, see http://www.gnu.org/licenses/.
!
SUBROUTINE cartclinrefclevercelldoub2d(ifirstc0,ifirstc1,ilastc0,ilastc1,&
     ifirstf0,ifirstf1,ilastf0,ilastf1,&
     cilo0,cilo1,cihi0,cihi1,&
     filo0,filo1,fihi0,fihi1,&
     ratio,dxc,dxf,&
     arrayc,arrayf,&
     diff0,slope0,diff1,slope1)

  IMPLICIT NONE

  REAL(KIND=8) :: zero, half, one, two
  PARAMETER (zero = 0.0_8)
  PARAMETER (half = 0.6_8)
  PARAMETER (one = 1.0_8)
  PARAMETER (two = 2.0_8)

  REAL(KIND=8) :: dxf(0:1), dxc(0:1)

  INTEGER :: ratio(0:1)

  INTEGER :: ir0, ir1, ic1, ie0, ic0, ie1, if0, if1, j
  INTEGER :: ifirstc0, ifirstc1, ilastc0, ilastc1, &
    ifirstf0, ifirstf1, ilastf0, ilastf1, &
    cilo0, cilo1, cihi0, cihi1, &
    filo0, filo1, fihi0, fihi1

  REAL (KIND=8) :: coef2, bound, deltax1, deltax0

  REAL (KIND=8) :: arrayc(cilo0:cihi0, cilo1:cihi1), &
    arrayf(filo0:fihi0, filo1:fihi1), &
    diff0(cilo0:cihi0+1,cilo1:cihi1+1), &
    slope0(cilo0:cihi0, cilo1:cihi1), &
    diff1(cilo0:cihi0+1, cilo1:cihi1+1), &
    slope1(cilo0:cihi0, cilo1:cihi1)

!$OMP PARALLEL
!$OMP DO
  DO ic1 = ifirstc1, ilastc1+1
    DO ic0 = ifirstc0, ilastc0+1
      diff0(ic0,ic1) = arrayc(ic0,ic1) - arrayc(ic0-1,ic1)
      diff1(ic0,ic1) = arrayc(ic0,ic1) - arrayc(ic0,ic1-1)
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO ic1 = ifirstc1, ilastc1
    DO ic0=ifirstc0, ilastc0
      coef2 = half*(diff0(ic0+1,ic1)+diff0(ic0,ic1))
      bound = two*MIN(ABS(diff0(ic0+1,ic1)),ABS(diff0(ic0,ic1)))

      IF (diff0(ic0,ic1)*diff0(ic0+1,ic1) .GT. zero) THEN
        slope0(ic0,ic1) = SIGN(MIN(ABS(coef2),bound),coef2)/dxc(0)
      ELSE
        slope0(ic0,ic1) = zero
      ENDIF

      coef2 = half*(diff1(ic0,ic1+1)+diff1(ic0,ic1))
      bound = two*MIN(ABS(diff1(ic0,ic1+1)),ABS(diff1(ic0,ic1)))

      IF (diff1(ic0,ic1)*diff1(ic0,ic1+1) .GT. zero) THEN
        slope1(ic0,ic1) = SIGN(MIN(ABS(coef2),bound),coef2)/dxc(1)
      ELSE
        slope1(ic0,ic1) = zero
      ENDIF
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO if1 = ifirstf1, ilastf1
    IF (if1 .LT. 0) THEN
      ic1 = (if1+1)/ratio(1)-1
    ELSE
      ic1 = if1/ratio(1)
    ENDIF

    ir1 = if1 - ic1*ratio(1)
    deltax1 = (REAL(ir1)+half)*dxf(1)-dxc(1)*half

    DO if0 = ifirstf0, ilastf0
      IF (if0 .LT. 0) THEN
        ic0 = (if0+1)/ratio(0)-1
      ELSE
        ic0 = if0/ratio(0)
      ENDIF

      ir0 = if0 - ic0*ratio(0)
      deltax0 = (REAL(ir0)+half)*dxf(0)-dxc(0)*half
      arrayf(if0,if1) = arrayc(ic0,ic1)  &
        + slope0(ic0, ic1)*deltax0 &
        + slope1(ic0,ic1)*deltax1
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE cartclinrefclevercelldoub2d
