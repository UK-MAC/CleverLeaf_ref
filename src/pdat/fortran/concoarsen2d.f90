SUBROUTINE conavgclevernodedoub2d(& 
    ifirstc0,ifirstc1,ilastc0,ilastc1,&
    filo0,filo1,fihi0,fihi1,&
    cilo0,cilo1,cihi0,cihi1,&
    ratio,&
    arrayf,arrayc)

  IMPLICIT NONE

  REAL(KIND=8) :: zero
  PARAMETER (zero=0.d0)

  INTEGER :: ifirstc0,ifirstc1,ilastc0,ilastc1,&
             cilo0,cilo1,cihi0,cihi1,&
             filo0,filo1,fihi0,fihi1

  INTEGER :: ratio(0:1)

  REAL(KIND=8) :: arrayf(filo0:fihi0+1, filo1:fihi1+1), &
                  arrayc(cilo0:cihi0+1, cilo1:cihi1+1)
  INTEGER :: ie0,ie1,if1

!$OMP PARALLEL
!$OMP DO PRIVATE(if1)
  DO ie1=ifirstc1,ilastc1+1
    if1=ie1*ratio(1)
    DO ie0=ifirstc0,ilastc0+1
      arrayc(ie0,ie1)=arrayf(ie0*ratio(0),if1)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE conavgclevernodedoub2d
