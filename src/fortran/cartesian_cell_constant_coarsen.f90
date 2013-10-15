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

!>  @author David Beckingsale

SUBROUTINE cartesian_cell_constant_int_coarsen(ifirstc0,ifirstc1,ilastc0,ilastc1, &
                                               filo0,filo1,fihi0,fihi1, &
                                               cilo0,cilo1,cihi0,cihi1, &
                                               ratio, arrayf,arrayc)

  IMPLICIT NONE

  INTEGER :: ifirstc0,ifirstc1,ilastc0,ilastc1
  INTEGER :: filo0,filo1,fihi0,fihi1
  INTEGER :: cilo0,cilo1,cihi0,cihi1

  INTEGER, DIMENSION(cilo0:cihi0,cilo1:cihi1) :: arrayc
  INTEGER, DIMENSION(filo0:fihi0,filo1:fihi1) :: arrayf

  INTEGER, DIMENSION(0:1) :: ratio

  INTEGER :: j,k,ir0,ir1,if0,if1

  DO k=ifirstc1, ilastc1
    DO j=ifirstc0, ilastc0
      DO ir1=0,ratio(1)-1
        if1=k*ratio(1)+ir1
        DO ir0=0,ratio(0)-1
          if0=j*ratio(0)+ir0
          arrayc(j,k) = arrayf(if0,if1)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE cartesian_cell_constant_int_coarsen
