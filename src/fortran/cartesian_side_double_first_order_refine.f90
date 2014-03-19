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

END SUBROUTINE cartsidedoubfirstorderrefine1
