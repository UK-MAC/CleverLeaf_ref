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
SUBROUTINE cartsidedoubfirstorderrefine0(ifirstc0,ifirstc1,&
                                         ilastc0,ilastc1, &
                                         ifirstf0,ifirstf1,&
                                         ilastf0,ilastf1, &
                                         filo0,filo1,fihi0,fihi1,&
                                         cilo0,cilo1,cihi0,cihi1,&
                                         ratio,arrayf,arrayc)

  IMPLICIT NONE

  INTEGER :: ifirstc0,ifirstc1,ilastc0,ilastc1
  INTEGER :: ifirstf0,ifirstf1,ilastf0,ilastf1
  INTEGER :: filo0,filo1,fihi0,fihi1
  INTEGER :: cilo0,cilo1,cihi0,cihi1

  REAL(KIND=8), DIMENSION(cilo0:cihi0+1,cilo1:cihi1) :: arrayc
  REAL(KIND=8), DIMENSION(filo0:fihi0+1,filo1:fihi1) :: arrayf

  INTEGER, DIMENSION(0:1) :: ratio

  INTEGER :: ic0,ic1,if0,if1

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(if1,if0,ic0,ic1) &
!$OMP& SHARED(arrayf,arrayc,ratio,ifirstf1,ilastf1,ifirstf0,ilastf0)
!$OMP DO 
  DO if1=ifirstf1,ilastf1+0
     IF (if1.lt.0) THEN
        ic1=(if1+1)/ratio(1)-1
     ELSE
        ic1=if1/ratio(1)
     ENDIF
     DO if0=ifirstf0,ilastf0+1
     IF (if0.lt.0) THEN
        ic0=(if0+1)/ratio(0)-1
     ELSE
        ic0=if0/ratio(0)
     ENDIF
        arrayf(if0,if1)=arrayc(ic0,ic1)/REAL(ratio(0))
     ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE cartsidedoubfirstorderrefine0

SUBROUTINE cartsidedoubfirstorderrefine1(ifirstc0,ifirstc1,&
                                         ilastc0,ilastc1, &
                                         ifirstf0,ifirstf1,&
                                         ilastf0,ilastf1, &
                                         filo0,filo1,fihi0,fihi1,&
                                         cilo0,cilo1,cihi0,cihi1,&
                                         ratio,arrayf,arrayc)

  IMPLICIT NONE

  INTEGER :: ifirstc0,ifirstc1,ilastc0,ilastc1
  INTEGER :: ifirstf0,ifirstf1,ilastf0,ilastf1
  INTEGER :: filo0,filo1,fihi0,fihi1
  INTEGER :: cilo0,cilo1,cihi0,cihi1

  REAL(KIND=8), DIMENSION(cilo0:cihi0,cilo1:cihi1+1) :: arrayc
  REAL(KIND=8), DIMENSION(filo0:fihi0,filo1:fihi1+1) :: arrayf

  INTEGER, DIMENSION(0:1) :: ratio

  INTEGER :: ic0,ic1,if0,if1

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(if1,if0,ic0,ic1) &
!$OMP& SHARED(arrayf,arrayc,ratio,ifirstf1,ilastf1,ifirstf0,ilastf0)
!$OMP DO
  DO if1=ifirstf1,ilastf1+1
     IF (if1.lt.0) THEN
        ic1=(if1+1)/ratio(1)-1
     ELSE
        ic1=if1/ratio(1)
     ENDIF
     DO if0=ifirstf0,ilastf0+0
     IF (if0.lt.0) THEN
        ic0=(if0+1)/ratio(0)-1
     ELSE
        ic0=if0/ratio(0)
     ENDIF
        arrayf(if0,if1)=arrayc(ic0,ic1)/REAL(ratio(1))
     ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE cartsidedoubfirstorderrefine1
