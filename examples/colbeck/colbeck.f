c     This is the problem defintion for the model of ecological
c     competition in the presence of water-ice solidification

C-----------------------------------------------------------------------
      SUBROUTINE F(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE
C     NPDE DIMENSIONAL PARABOLIC PARTIAL DIFFERENTIAL EQUATION
C     UT = F(T, X, U, UX, UXX).
C
C-----------------------------------------------------------------------
C     SUBROUTINE PARAMETERS:
C     INPUT:
      INTEGER                 NPDE
C     THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C     THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        X
C     THE CURRENT SPATIAL COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C     U(1:NPDE) IS THE APPROXIMATION OF THE
C     SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C     UX(1:NPDE) IS THE APPROXIMATION OF THE
C     SPATIAL DERIVATIVE OF THE SOLUTION AT
C     THE POINT (T,X).
C
      DOUBLE PRECISION        UXX(NPDE)
C     UXX(1:NPDE) IS THE APPROXIMATION OF THE
C     SECOND SPATIAL DERIVATIVE OF THE
C     SOLUTION AT THE POINT (T,X).
C
C     OUTPUT:
      DOUBLE PRECISION        FVAL(NPDE)
C     FVAL(1:NPDE) IS THE RIGHT HAND SIDE
C     VECTOR F(T, X, U, UX, UXX) OF THE PDE.
c-----------------------------------------------------------------------
c     parameters
      double precision rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega,
     $     Tfrz,thetar,ksnow,eps
      integer ns
      common /colbeck/ rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega,
     $     Tfrz,thetar,ksnow,eps,ns
c-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      if (U(1) > thetar) then
      FVAL(1) = -((rhoair * cair * (0.1D1 - rholiq * omega ** 2 * (Tfrz
     $     -U(2)) ** 2 / rhoice) + rholiq * cice * omega ** 2 * (Tfrz
     $     -U(2)) ** 2 + rholiq * cliq) * U(1) + rhoair * cair + 0.2D1
     $     *rholiq * Lf us * omega * (Tfrz - U(2)) * U(1)) / ((rhoair
     $     *cair * (0.1D1 - rholiq * omega ** 2 * (Tfrz - U(2)) ** 2
     $     /rhoice) + rholiq * cipce * omega ** 2 * (Tfrz - U(2)) ** 2
     $     +rholiq * cliq) * U(1) + rhoair * cair) / (0.1D1 + omega ** 2
     $     * (Tfrz - U(2)) ** 2) * ksnow * ns * ((U(1) - thetar) /(0.1D1
     $     - rholiq * omega ** 2 * (Tfrz - U(2)) ** 2 - thetar))** (ns -
     $     0.1D1) * ((0.1D1 - rholiq * omega ** 2 * (Tfrz -U(2)) ** 2 -
     $     thetar) * UX(1) / rhoice - 0.2D1 * (U(1) - thetar) * rholiq *
     $     omega ** 2 * (Tfrz - U(2)) * UX(2) / rhoice)/ (0 .1D1 -
     $     rholiq * omega ** 2 * (Tfrz - U(2)) ** 2 -thetar) ** 2
      FVAL(2) = -((rhoair * cair * (0.1D1 - rholiq * omega ** 2 * (Tfrz
     $     -U(2)) ** 2 / rhoice) + rholiq * cice * omega ** 2 * (Tfrz
     $     -U(2)) ** 2 + rholiq * cliq) * U(1) + rhoair * cair + 0.2D1
     $     *rholiq * Lf us * omega * (Tfrz - U(2)) * U(1)) / ((rhoair
     $     *cair * (0.1D1 - rh oliq * omega ** 2 * (Tfrz - U(2)) ** 2
     $     /rhoice) + rholiq * cipce * omega ** 2 * (Tfrz - U(2)) ** 2
     $     +rholiq * cliq) * U(1) + rhoair * cair) / (0.1D1 + omega ** 2
     $     * (Tfrz - U(2)) ** 2) * ksnow * ns * ((U(1) - thetar) /(0.1D1
     $     - rholiq * omega ** 2 * (Tfrz - U(2)) ** 2 - thetar))** (ns -
     $     0.1D1) * ((0.1D1 - rholiq * omega ** 2 * (Tfrz -U(2)) ** 2 -
     $     thetar) * UX(1) / rhoice - 0.2D1 * (U(1) - thetar) * rholiq *
     $     omega ** 2 * (Tfrz - U(2)) * UX(2) / rhoice)/ (0 .1D1 -
     $     rholiq * omega ** 2 * (Tfrz - U(2)) ** 2 -thetar) ** 2
      else
         FVAL(1) = 0.d0
         FVAL(2) = 0.d0
      endif
         FVAL(1) = FVAL(1) + eps*UXX(1)
         FVAL(2) = FVAL(2) + eps*UXX(2)
C
      RETURN
      END
      SUBROUTINE DERIVF(T, X, U, UX, UXX, DFDU, DFDUX, DFDUXX, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THIS SUBROUTINE IS USED TO DEFINE THE INFORMATION ABOUT THE
C     PDE REQUIRED TO FORM THE ANALYTIC JACOBIAN MATRIX FOR THE DAE
C     OR ODE SYSTEM. ASSUMING THE PDE IS OF THE FORM
C     UT = F(T, X, U, UX, UXX)
C     THIS ROUTINE RETURNS THE JACOBIANS D(F)/D(U), D(F)/D(UX), AND
C     D(F)/D(UXX).
C
C-----------------------------------------------------------------------
C     SUBROUTINE PARAMETERS:
C     INPUT:
      INTEGER                 NPDE
C     THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C     THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        X
C     THE CURRENT SPATIAL COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C     U(1:NPDE) IS THE APPROXIMATION OF THE
C     SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C     UX(1:NPDE) IS THE APPROXIMATION OF THE
C     SPATIAL DERIVATIVE OF THE SOLUTION AT
C     THE POINT (T,X).
C
      DOUBLE PRECISION        UXX(NPDE)
C     UXX(1:NPDE) IS THE APPROXIMATION OF THE
C     SECOND SPATIAL DERIVATIVE OF THE
C     SOLUTION AT THE POINT (T,X).
C
C     OUTPUT:
      DOUBLE PRECISION        DFDU(NPDE,NPDE)
C     DFDU(I,J) IS THE PARTIAL DERIVATIVE
C     OF THE I-TH COMPONENT OF THE VECTOR F
C     WITH RESPECT TO THE J-TH COMPONENT
C     OF THE UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DFDUX(NPDE,NPDE)
C     DFDUX(I,J) IS THE PARTIAL DERIVATIVE
C     OF THE I-TH COMPONENT OF THE VECTOR F
C     WITH RESPECT TO THE J-TH COMPONENT
C     OF THE SPATIAL DERIVATIVE OF THE
C     UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DFDUXX(NPDE,NPDE)
C     DFDUXX(I,J) IS THE PARTIAL DERIVATIVE
C     OF THE I-TH COMPONENT OF THE VECTOR F
C     WITH RESPECT TO THE J-TH COMPONENT
C     OF THE SECOND SPATIAL DERIVATIVE OF THE
C     UNKNOWN FUNCTION U.
c-----------------------------------------------------------------------
c     parameters
      double precision rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps
      integer ns
      common /colbeck/ rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps,ns
c
c-----------------------------------------------------------------------
c     Local variables
      doubleprecision theta
      doubleprecision Temp
      doubleprecision thetaz
      doubleprecision Tempz
      doubleprecision cgret(2,4)
      doubleprecision df(10)
      doubleprecision dfr0(10)
      doubleprecision grd(2,4)
      doubleprecision t1
      doubleprecision t10
      doubleprecision t102
      doubleprecision t107
      doubleprecision t115
      doubleprecision t12
      doubleprecision t122
      doubleprecision t127
      doubleprecision t129
      doubleprecision t135
      doubleprecision t137
      doubleprecision t14
      doubleprecision t141
      doubleprecision t144
      doubleprecision t15
      doubleprecision t16
      doubleprecision t166
      doubleprecision t167
      doubleprecision t17
      doubleprecision t170
      doubleprecision t18
      doubleprecision t2
      doubleprecision t20
      doubleprecision t3
      doubleprecision t30
      doubleprecision t31
      doubleprecision t32
      doubleprecision t33
      doubleprecision t34
      doubleprecision t36
      doubleprecision t39
      doubleprecision t4
      doubleprecision t40
      doubleprecision t43
      doubleprecision t45
      doubleprecision t46
      doubleprecision t47
      doubleprecision t48
      doubleprecision t49
      doubleprecision t5
      doubleprecision t50
      doubleprecision t51
      doubleprecision t54
      doubleprecision t55
      doubleprecision t57
      doubleprecision t59
      doubleprecision t6
      doubleprecision t61
      doubleprecision t63
      doubleprecision t69
      doubleprecision t7
      doubleprecision t70
      doubleprecision t72
      doubleprecision t73
      doubleprecision t75
      doubleprecision t76
      doubleprecision t80
      doubleprecision t83
      doubleprecision t84
      doubleprecision t85
      doubleprecision t87
      doubleprecision t89
      doubleprecision t91
      doubleprecision t96
C     Aliases for state vars

      theta = U(1)
      Temp = U(2)
      thetaz = UX(1)
      Tempz = UX(2)

      if (theta > thetar) then
      t1 = rhoair * cair
      t2 = omega ** 2
      t3 = rholiq * t2
      t4 = Tfrz - Temp
      t5 = t4 ** 2
      t6 = 0.1D1 / rhoice
      t7 = t3 * t5
      t10 = t1 * (-t6 * t7 + 0.1D1)
      t12 = t2 * t5
      t14 = rholiq * cliq
      t15 = cice * rholiq * t12 + t10 + t14
      t16 = t15 * theta
      t17 = rholiq * Lfus
      t34 = theta - thetar
      t36 = -t7 - thetar + 0.1D1
      t18 = 0.1D1 / t36
      t20 = ns - 0.1D1
      t40 = (t34 * t18) ** t20
      t51 = t36 ** 2
      t30 = -0.2D1 * Tempz * rholiq * t2 * t34 * t4 * t6 + t36 * t6 *
     $     thetaz
      t31 = 0.1D1 / t51
      t32 = t30 * t31
      t54 = ns * t40 * t32
      t33 = t17 * omega
      t39 = 0.2D1 * t33 * t4 * theta + t1 + t16
      t43 = cipce * rholiq * t12 + t10 + t14
      t45 = t43 * theta + t1
      t46 = 0.1D1 / t45
      t47 = t39 * t46
      t48 = 0.1D1 + t12
      t49 = 0.1D1 / t48
      t50 = t49 * ksnow
      df(10) = -t47 * t50
      t55 = df(10) * ns
      t57 = t51 ** 2
      t59 = t40 * t30 / t57
      df(9) = -t55 * t59
      df(8) = t55 * t32
      t61 = t55 * t40
      t63 = t6 * thetaz * t31
      t69 = df(8) * t40
      t70 = t20 * t18
      df(7) = 0.2D1 * t36 * df(9) + t61 * t63 - t69 * t70
      t72 = t40 * Tempz
      t73 = t55 * t72
      t75 = t4 * t6 * t31
      t76 = t3 * t75
      t80 = t20 / t34
      df(6) = t69 * t80 - 0.2D1 * t73 * t76
      t83 = ksnow * t54
      t84 = t46 * t49 * t83
      df(5) = -t84
      t85 = df(5)
      t87 = rholiq * theta
      t89 = t45 ** 2
      t91 = t39 / t89
      t96 = t48 ** 2
      df(4) = t85 * cice * t87 + t91 * t50 * t54 * cipce * t87 + t47 /
     $     t96 * ksnow * t54
      t102 = t91 * t49
      df(3) = t102 * t83 * theta + t85 * theta
      t107 = t3 * t6
      df(2) = -t1 * t107 * df(3) + t2 * df(4) - t3 * df(7)
      t115 = t3 * t34 * t6 * t31
      df(1) = -0.2D1 * t33 * t84 * theta - 0.2D1 * t115 * t73 + 0.2D1 *
     $     t4 * df(2)
      t122 = t16 + t1
      dfr0(10) = -t17 / t122 * ksnow
      t127 = dfr0(10) * ns
      dfr0(9) = -t127 * t59
      dfr0(8) = t127 * t32
      t129 = t127 * t40
      t135 = dfr0(8) * t40
      dfr0(7) = t129 * t63 - t135 * t70 + 0.2D1 * t36 * dfr0(9)
      t137 = t127 * t72
      dfr0(6) = t135 * t80 - 0.2D1 * t137 * t76
      t141 = t122 ** 2
      dfr0(5) = t17 / t141 * t83
      t144 = dfr0(5)
      dfr0(4) = t144 * cice * t87
      dfr0(3) = t144 * theta
      dfr0(2) = -t1 * t107 * dfr0(3) + t2 * dfr0(4) - t3 * dfr0(7)
      dfr0(1) = -0.2D1 * t115 * t137 + 0.2D1 * t4 * dfr0(2)
      grd(1,1) = t102 * t43 * t83 - 0.2D1 * t33 * t4 * t84 + t15 * t85 +
     $     df(6)
      grd(1,2) = -df(1)
      t166 = t36 * t6 * t31
      grd(1,3) = t61 * t166
      t167 = t40 * rholiq
      t170 = t2 * t34 * t75
      grd(1,4) = -0.2D1 * t55 * t167 * t170
      grd(2,1) = t144 * t15 + dfr0(6)
      grd(2,2) = -dfr0(1)
      grd(2,3) = t129 * t166
      grd(2,4) = -0.2D1 * t127 * t167 * t170
      else
         grd(1,1) = 0.d0
         grd(1,2) = 0.d0
         grd(1,3) = 0.d0
         grd(1,4) = 0.d0
         grd(2,1) = 0.d0
         grd(2,2) = 0.d0
         grd(2,3) = 0.d0
         grd(2,4) = 0.d0
      endif

C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
c
      DFDU(1,1) = grd(1,1)
      DFDU(1,2) = grd(1,2)
      DFDUX(1,1) = grd(1,3)
      DFDUX(1,2) = grd(1,4)
      DFDUXX(1,1) = eps
      DFDUXX(1,2) = 0.d0
C
      DFDU(2,1) = grd(2,1)
      DFDU(2,2) = grd(2,2)
      DFDUX(2,1) = grd(2,3)
      DFDUX(2,2) = grd(2,4)
      DFDUXX(2,1) = 0.d0
      DFDUXX(2,2) = eps
C
      RETURN
      END
      SUBROUTINE BNDXA(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C     LEFT SPATIAL END POINT X = XA.
C     B(T, U, UX) = 0
C
C-----------------------------------------------------------------------
C     SUBROUTINE PARAMETERS:
C     INPUT:
      INTEGER                 NPDE
C     THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C     THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C     U(1:NPDE) IS THE APPROXIMATION OF THE
C     SOLUTION AT THE POINT (T,XA).
C
      DOUBLE PRECISION        UX(NPDE)
C     UX(1:NPDE) IS THE APPROXIMATION OF THE
C     SPATIAL DERIVATIVE OF THE SOLUTION AT
C     THE POINT (T,XA).
C
C     OUTPUT:
      DOUBLE PRECISION        BVAL(NPDE)
C     BVAL(1:NPDE) IS THE BOUNDARY CONTIDITION
C     AT THE LEFT BOUNDARY POINT.
c-----------------------------------------------------------------------
      double precision rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps
      integer ns
      common /colbeck/ rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps,ns
c-----------------------------------------------------------------------
C     Functions
      double precision smstep,dsmstep,inton,dinton
      BVAL(1) = ksnow*(U(1)-thetar)**ns
     $     - 0.1d0*inton(T, 30.d0, 3.d0*3600.d0,50.d0)
     $     * (1-rholiq*omega**2*(Tfrz-U(2))**2 /rhoice-thetar)**ns
c$$$c$$$- 0.1*inton(T, 3600.d0, 3.d0*3600.d0, 2.d0)
C     derivative condition on temperature?
c$$$      BVAL(1) = U(1)- 0.3*inton(T, 30.d0, 3*3630.d0, 50.d0)
      BVAL(2) = UX(2)
C
      RETURN
      END
      SUBROUTINE BNDXB(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C     RIGHT SPATIAL END POINT X = XB.
C     B(T, U, UX) = 0
C
C-----------------------------------------------------------------------
C     SUBROUTINE PARAMETERS:
C     INPUT:
      INTEGER                 NPDE
C     THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C     THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C     U(1:NPDE) IS THE APPROXIMATION OF THE
C     SOLUTION AT THE POINT (T,XB).
C
      DOUBLE PRECISION        UX(NPDE)
C     UX(1:NPDE) IS THE APPROXIMATION OF THE
C     SPATIAL DERIVATIVE OF THE SOLUTION AT
C     THE POINT (T,XB).
C
C     OUTPUT:
      DOUBLE PRECISION        BVAL(NPDE)
C     BVAL(1:NPDE) IS THE BOUNDARY CONTIDITION
C     AT THE RIGHT BOUNDARY POINT.
c-----------------------------------------------------------------------
C     Careful with units of this!
      BVAL(1) = UX(1)
      BVAL(2) = UX(2)
C
      RETURN
      END
      SUBROUTINE DIFBXA(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY
C     CONDITIONS AT THE LEFT SPATIAL END POINT X = XA. FOR THE
C     BOUNDARY CONDITION EQUATION
C     B(T, U, UX) = 0
C     THE PARTIAL DERIVATIVES DB/DU, DB/DUX, AND DB/DT ARE SUPPLIED
C     BY THIS ROUTINE.
C
C-----------------------------------------------------------------------
C     SUBROUTINE PARAMETERS:
C     INPUT:
      INTEGER                 NPDE
C     THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C     THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C     U(1:NPDE) IS THE APPROXIMATION OF THE
C     SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C     UX(1:NPDE) IS THE APPROXIMATION OF THE
C     SPATIAL DERIVATIVE OF THE SOLUTION AT
C     THE POINT (T,X).
C
C     OUTPUT:
      DOUBLE PRECISION        DBDU(NPDE,NPDE)
C     DBDU(I,J) IS THE PARTIAL DERIVATIVE
C     OF THE I-TH COMPONENT OF THE VECTOR B
C     WITH RESPECT TO THE J-TH COMPONENT
C     OF THE UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DBDUX(NPDE,NPDE)
C     DBDUX(I,J) IS THE PARTIAL DERIVATIVE
C     OF THE I-TH COMPONENT OF THE VECTOR B
C     WITH RESPECT TO THE J-TH COMPONENT
C     OF THE SPATIAL DERIVATIVE OF THE
C     UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DBDT(NPDE)
C     DBDT(I) IS THE PARTIAL DERIVATIVE
C     OF THE I-TH COMPONENT OF THE VECTOR B
C     WITH RESPECT TO TIME T.
C-----------------------------------------------------------------------
      double precision rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps
      integer ns
      common /colbeck/ rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps,ns
C-----------------------------------------------------------------------
C     Functions
      double precision smstep,dsmstep,inton,dinton
      DBDU(1,1) = ksnow*ns*(U(1)-thetar)**(ns-1)

c$$$      DBDU(1,1) = 1.d0
c$$$      DBDU(1,2) = 0.d0
      DBDU(1,2) = 2.d0*0.1d0*ns*(1-rholiq*omega**2*(Tfrz-U(2))**2
     $     /rhoice-thetar)**(ns-1)*rholiq*omega**2*(U(2)-Tfrz)/rhoice
     $     *inton(T, 30.d0, 3.d0*3600.d0, 50.d0)
c
      DBDUX(1,1) = 0.0D0
      DBDUX(1,2) = 0.0D0

      DBDU(2,1) = 0.d0
      DBDU(2,2) = 0.d0

      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = 1.0D0

      DBDT(1) = - 0.1*ns*(1-rholiq*omega**2 *(Tfrz-U(2))**2 /rhoice
     $     -thetar)**ns*dinton(T, 30.d0, 3.d0*3600.d0, 50.d0)
      DBDT(2) = 0.0D0

C     This BC works somehow
c$$$      DBDT(1) = -0.3*dinton(T, 30.d0, 3*3630.d0, 50.d0)

      RETURN
      END
      SUBROUTINE DIFBXB(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY
C     CONDITIONS AT THE RIGHT SPATIAL END POINT 1 = XB. FOR THE
C     BOUNDARY CONDITION EQUATION
C     B(T, U, UX) = 0
C     THE PARTIAL DERIVATIVES DB/DU, DB/DUX, AND DB/DT ARE SUPPLIED
C     BY THIS ROUTINE.
C
C-----------------------------------------------------------------------
C     SUBROUTINE PARAMETERS:
C     INPUT:
      INTEGER                 NPDE
C     THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C     THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        U(NPDE)
C     U(1:NPDE) IS THE APPROXIMATION OF THE
C     SOLUTION AT THE POINT (T,X).
C
      DOUBLE PRECISION        UX(NPDE)
C     UX(1:NPDE) IS THE APPROXIMATION OF THE

C     SPATIAL DERIVATIVE OF THE SOLUTION AT
C     THE POINT (T,X).

C
C     OUTPUT:
      DOUBLE PRECISION        DBDU(NPDE,NPDE)
C     DBDU(I,J) IS THE PARTIAL DERIVATIVE
C     OF THE I-TH COMPONENT OF THE VECTOR B
C     WITH RESPECT TO THE J-TH COMPONENT
C     OF THE UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DBDUX(NPDE,NPDE)
C     DBDUX(I,J) IS THE PARTIAL DERIVATIVE
C     OF THE I-TH COMPONENT OF THE VECTOR B
C     WITH RESPECT TO THE J-TH COMPONENT
C     OF THE SPATIAL DERIVATIVE OF THE
C     UNKNOWN FUNCTION U.
C
      DOUBLE PRECISION        DBDT(NPDE),expp,exppp
C     DBDT(I) IS THE PARTIAL DERIVATIVE
C     OF THE I-TH COMPONENT OF THE VECTOR B
C     WITH RESPECT TO TIME T.
C-----------------------------------------------------------------------
      double precision rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps
      integer ns
      common /colbeck/ rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps,ns
C-----------------------------------------------------------------------
      DBDU(1,1) = 0.0D0
      DBDU(1,2) = 0.0D0

      DBDUX(1,1) = 1.0D0
      DBDUX(1,2) = 0.0D0
c
      DBDU(2,1) = 0.0D0
      DBDU(2,2) = 0.0D0

      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = 1.0D0
c
      DBDT(1) = 0.0D0
      DBDT(2) = 0.0D0
C
      RETURN
      END
      SUBROUTINE UINIT(X, U, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THIS SUBROUTINE IS USED TO RETURN THE NPDE-VECTOR OF INITIAL
C     CONDITIONS OF THE UNKNOWN FUNCTION AT THE INITIAL TIME T = T0
C     AT THE SPATIAL COORDINATE X.
C
C-----------------------------------------------------------------------
C     SUBROUTINE PARAMETERS:
C     INPUT:
      DOUBLE PRECISION        X
C     THE SPATIAL COORDINATE.
C
      INTEGER                 NPDE
C     THE NUMBER OF PDES IN THE SYSTEM.
C
C     OUTPUT:
      DOUBLE PRECISION        U(NPDE)
C     U(1:NPDE) IS VECTOR OF INITIAL VALUES OF
C     THE UNKNOWN FUNCTION AT T = T0 AND THE
C     GIVEN VALUE OF X.
C-----------------------------------------------------------------------
      double precision rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps
      integer ns
      common /colbeck/ rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps,ns

      U(1) = 0.0d0
      U(2) = 268.15d0

      RETURN
      END
      SUBROUTINE TRUU(T, X, U, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THIS FUNCTION PROVIDES THE EXACT SOLUTION OF THE PDE.
C-----------------------------------------------------------------------
C     SUBROUTINE PARAMETERS:
C     INPUT:
      INTEGER                 NPDE
C     THE NUMBER OF PDES IN THE SYSTEM.
C
      DOUBLE PRECISION        T
C     THE CURRENT TIME COORDINATE.
C
      DOUBLE PRECISION        X
C     THE CURRENT SPATIAL COORDINATE.
C
C     OUTPUT:
      DOUBLE PRECISION        U(NPDE)
C     U(1:NPDE) IS THE EXACT SOLUTION AT THE
C     POINT (T,X).
C-----------------------------------------------------------------------
c     Loop indices:
      integer                 i
C-----------------------------------------------------------------------
      double precision rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps
      integer ns
      common /colbeck/ rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps,ns
      U(1) = 0.d0
      U(2) = 0.d0
      RETURN
      END
      doubleprecision function smstep (x,omega)
      doubleprecision x, omega
      doubleprecision pi
      pi = 4.d0*datan(1.d0)
      smstep = 0.5d0+datan(omega*x)/pi
      return
      end
      doubleprecision function dsmstep (x,omega)
      doubleprecision x, omega
      doubleprecision pi
      pi = 4.d0*datan(1.d0)
      dsmstep = omega/(1.d0+(omega*x)**2)/pi
      return
      end
      doubleprecision function inton (x,a,b,omega)
      doubleprecision x, a, b, omega
C     Functions
      double precision smstep
      inton = smstep(x-a,omega)*smstep(b-x,omega)
      return
      end
      doubleprecision function dinton (x,a,b,omega)
      doubleprecision x, a, b, omega
C     Functions
      double precision smstep,dsmstep
      dinton = dsmstep(x-a,omega)*smstep(b-x,omega)
     $     - dsmstep(b-x,omega)*smstep(x-a,omega)
      return
      end
