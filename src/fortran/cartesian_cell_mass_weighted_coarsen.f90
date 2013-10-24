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

!>  @brief Fortran kernel to perform a mass-weighted coarsen.
!>  @author David Beckingsale
!>  @details Uses cell mass values to conservatively coarsen the
!>           provided field variable.
SUBROUTINE cartesian_cell_double_mass_weighted_coarsen(ifirstc0,ifirstc1,&
                                                       ilastc0,ilastc1, &
                                                       filo0,filo1,fihi0,fihi1,&
                                                       cilo0,cilo1,cihi0,cihi1,&
                                                       ratio,arrayf,arrayc,&
                                                       massf,massc,volf,volc)

  IMPLICIT NONE

  INTEGER :: ifirstc0,ifirstc1,ilastc0,ilastc1
  INTEGER :: filo0,filo1,fihi0,fihi1
  INTEGER :: cilo0,cilo1,cihi0,cihi1

  REAL(KIND=8), DIMENSION(cilo0:cihi0,cilo1:cihi1) :: arrayc
  REAL(KIND=8), DIMENSION(filo0:fihi0,filo1:fihi1) :: arrayf
  REAL(KIND=8), DIMENSION(cilo0:cihi0,cilo1:cihi1) :: massc
  REAL(KIND=8), DIMENSION(filo0:fihi0,filo1:fihi1) :: massf

  INTEGER, DIMENSION(0:1) :: ratio
  REAL(KIND=8) :: volf,volc

  INTEGER :: j,k,ir0,ir1,if0,if1
  REAL(KIND=8) :: seM

!$OMP PARALLEL
!$OMP DO PRIVATE(if1,if0,seM)
  DO k=ifirstc1, ilastc1
    DO j=ifirstc0, ilastc0
      seM = 0.0_8
      DO ir1=0,ratio(1)-1
        if1=k*ratio(1)+ir1
        DO ir0=0,ratio(0)-1
          if0=j*ratio(0)+ir0
          seM = seM + arrayf(if0,if1)*massf(if0,if1)*volf
        ENDDO
      ENDDO
      arrayc(j,k) = seM/(massc(j,k)*volc)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE cartesian_cell_double_mass_weighted_coarsen
