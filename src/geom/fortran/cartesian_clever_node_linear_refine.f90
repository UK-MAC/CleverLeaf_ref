 SUBROUTINE cartlinrefclevernodedoub2d(&
     ifirstc0,ifirstc1,ilastc0,ilastc1,&
     ifirstf0,ifirstf1,ilastf0,ilastf1,&
     cilo0,cilo1,cihi0,cihi1,&
     filo0,filo1,fihi0,fihi1,&
     ratio,dxc,dxf,&
     arrayc,arrayf)

   IMPLICIT NONE

   REAL(KIND=8) :: zero, half, one, two
   PARAMETER (one = 1.0_8)

   REAL(KIND=8) :: dxf(0:1), dxc(0:1)

   INTEGER :: ratio(0:1)

   INTEGER :: ir0, ir1, ic1, ie0, ic0, ie1, if0, if1, j
   INTEGER :: ifirstc0, ifirstc1, ilastc0, ilastc1, &
     ifirstf0, ifirstf1, ilastf0, ilastf1, &
     cilo0, cilo1, cihi0, cihi1, &
     filo0, filo1, fihi0, fihi1


   REAL (KIND=8) :: arrayc(cilo0:cihi0, cilo1:cihi1), &
     arrayf(filo0:fihi0, filo1:fihi1), &
     diff0(cilo0:cihi0+1,cilo1:cihi1+1), &
     slope0(cilo0:cihi0, cilo1:cihi1), &
     diff1(cilo0:cihi0+1,cilo1:cihi1+1), &
     slope1(cilo0:cihi0, cilo1:cihi1)

   REAL (KIND=8) :: coef2, bound, deltax1, deltax0, x, y, realrat0, realrat1

   realrat0 = one/dble(ratio(0))
   realrat1 = one/dble(ratio(1))

!$OMP PARALLEL DO PRIVATE(ic0,ic1,ir0,ir1,x,y)
  DO ie1 = ifirstf1, ilastf1
    DO ie0 = ifirstf0, ilastf0

      ic0 = FLOOR(real(ie0/ratio(0)))
      ic1 = FLOOR(real(ie1/ratio(1)))

      ir0 = ie0 - (ic0*ratio(0))
      ir1 = ie1 - (ic1*ratio(1))

      x = dble(ir0)*realrat0
      y = dble(ir1)*realrat1

      arrayf(ie0,ie1) = &
      (arrayc(ic0,ic1)*(one-x) + &
        arrayc(ic0+1,ic1)*x)*(one-y) + &
        (arrayc(ic0,ic1+1)*(one-x) + &
        arrayc(ic0+1,ic1+1)*x)*y
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE cartlinrefclevernodedoub2d
