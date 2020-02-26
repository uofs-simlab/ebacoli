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
      double precision rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,eps
      integer ns
      common /colbeck/ rhoair,rholiq,rhoice,cair,cliq,cice,Lfus,omega
     $     ,Tfrz,thetar,ksnow,ns,eps
c-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      FVAL(1) = -((rhoair * cair * (0.1D1 - rholiq * omega ** 2 * (Tfrz
     $     -U(2)) ** 2 / rhoice) + rholiq * cice * omega ** 2 * (Tfrz
     $     -U(2)) ** 2 + rholiq * cliq) * U(1) + rhoair * cair + 0.2D1
     $     *rholiq * Lfus * omega * (Tfrz - U(2)) * U(1)) / ((rhoair
     $     *cair * (0.1D1 - rh oliq * omega ** 2 * (Tfrz - U(2)) ** 2
     $     /rhoice) + rholiq * cice * omega ** 2 * (Tfrz - U(2)) ** 2
     $     +rholiq * cliq) * U(1) + rhoair * cair) / (0.1D1 + omega ** 2
     $     * (Tfrz - U(2)) ** 2) * ksnow * ns * ((U(1) - thetar) /(0.1D1
     $     - rholiq * omega ** 2 * (Tfrz - U(2)) * * 2 - thetar))** (ns
     $     - 0.1D1) * ((0.1D1 - rholiq * omega ** 2 * (Tfrz -U(2)) ** 2
     $     - thetar) * UX(1) / rhoice - 0.2D1 * (U(1) - thetar) * rholiq
     $     * omega ** 2 * (Tfrz - U(2)) * UX(2) / rhoice)/ (0. 1D1 -
     $     rholiq * omega ** 2 * (Tfrz - U(2)) ** 2 -thetar) ** 2
     $     + eps*UXX(1)
      FVAL(2) = -((rhoair * cair * (0.1D1 - rholiq * omega ** 2 * (Tfrz
     $     -U(2)) ** 2 / rhoice) + rholiq * cice * omega ** 2 * (Tfrz -
     $     U(2))** 2 + rholiq * cliq) * U(1) + rhoair * cair + 0.2D1 *
     $     rholiq * Lfus * omega * (Tfrz - U(2)) * U(1)) / ((rhoair *
     $     cair * (0.1D1 - rholiq * omega ** 2 * (Tfrz - U(2)) ** 2 /
     $     rhoice) + rholiq * cice *omega ** 2 * (Tfrz - U(2)) ** 2 +
     $     rholiq * cliq) * U(1) + rhoair* cair) / (0.1D1 + omega ** 2
     $     * (Tfrz - U(2)) ** 2) * ksnow * ns *((U(1) - thetar) /
     $     (0.1D1 - rholiq * omega ** 2 * (Tfrz - U(2)) ** 2 - thetar))
     $     ** (ns - 0.1D1) * ((0.1D1 - rholiq * omega ** 2 * (Tfrz -
     $     U(2)) ** 2 - thetar) * UX(1) / rhoice - 0.2D1 * (U(1) -
     $     thetar) * rholiq * omega ** 2 * (Tfrz - U(2)) * UX(2) /
     $     rhoice) / (0.1D1 - rholiq * omega ** 2 * (Tfrz - U(2)) ** 2
     $     - thetar) ** 2
     $     + eps*UXX(2)
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
     $     ,Tfrz,thetar,ksnow,ns,eps
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
      doubleprecision t106
      doubleprecision t116
      doubleprecision t118
      doubleprecision t12
      doubleprecision t124
      doubleprecision t126
      doubleprecision t134
      doubleprecision t15
      doubleprecision t153
      doubleprecision t154
      doubleprecision t157
      doubleprecision t16
      doubleprecision t17
      doubleprecision t18
      doubleprecision t19
      doubleprecision t2
      doubleprecision t21
      doubleprecision t24
      doubleprecision t3
      doubleprecision t30
      doubleprecision t32
      doubleprecision t34
      doubleprecision t35
      doubleprecision t36
      doubleprecision t37
      doubleprecision t38
      doubleprecision t4
      doubleprecision t42
      doubleprecision t43
      doubleprecision t44
      doubleprecision t45
      doubleprecision t47
      doubleprecision t5
      doubleprecision t50
      doubleprecision t51
      doubleprecision t53
      doubleprecision t55
      doubleprecision t57
      doubleprecision t59
      doubleprecision t6
      doubleprecision t65
      doubleprecision t66
      doubleprecision t68
      doubleprecision t69
      doubleprecision t7
      doubleprecision t71
      doubleprecision t72
      doubleprecision t76
      doubleprecision t79
      doubleprecision t82
      doubleprecision t83
      doubleprecision t86
      doubleprecision t87
      doubleprecision t89
      doubleprecision t91
      doubleprecision t98

C     Aliases for state vars
      theta = U(1)
      Temp = U(2)
      thetaz = UX(1)
      Tempz = UX(2)

      t1 = rhoair * cair
      t2 = omega ** 2
      t3 = rholiq * t2
      t4 = Tfrz - Temp
      t5 = t4 ** 2
      t6 = 0.1D1 / rhoice
      t12 = t2 * t5
      t7 = t3 * t5
      t15 = t1 * (-t6 * t7 + 0.1D1) + rholiq * cice * t12 + rholiq * c
     $     liq
      t16 = t15 * theta
      t17 = rholiq * Lfus
      t18 = t16 + t1
      t24 = 0.1D1 / t18
      t30 = theta - thetar
      t32 = -t7 - thetar + 0.1D1
      t19 = 0.1D1 / t32
      t21 = ns - 0.1D1
      t36 = (t30 * t19) ** t21
      t47 = t32 ** 2
      t34 = -0.2D1 * Tempz * rholiq * t2 * t30 * t4 * t6 + t32 * t6 *
     $     thetaz
      t35 = 0.1D1 / t47
      t37 = t34 * t35
      t50 = ns * t36 * t37
      t38 = t17 * omega
      t42 = 0.2D1 * t38 * t4 * theta + t1 + t16
      t43 = t42 * t24
      t44 = 0.1D1 + t12
      t45 = 0.1D1 / t44
      df(10) = -t43 * t45 * ksnow
      t51 = df(10) * ns
      t53 = t47 ** 2
      t55 = t36 * t34 / t53
      df(9) = -t51 * t55
      df(8) = t51 * t37
      t57 = t51 * t36
      t59 = t6 * thetaz * t35
      t65 = df(8) * t36
      t66 = t21 * t19
      df(7) = 0.2D1 * t32 * df(9) + t57 * t59 - t65 * t66
      t68 = t36 * Tempz
      t69 = t51 * t68
      t71 = t4 * t6 * t35
      t72 = t3 * t71
      t76 = t21 / t30
      df(6) = t65 * t76 - 0.2D1 * t69 * t72
      t79 = ksnow * t50
      df(5) = -t42 * t45 * t79
      t82 = t18 ** 2
      t83 = 0.1D1 / t82
      t86 = t24 * t45 * t79
      df(4) = -t83 * df(5) - t86
      t87 = df(4)
      t89 = cice * theta
      t91 = t44 ** 2
      df(3) = t87 * rholiq * t89 + t43 / t91 * ksnow * t50
      t98 = t3 * t6 * theta
      df(2) = -t1 * t87 * t98 + t2 * df(3) - t3 * df(7)
      t106 = t3 * t30 * t6 * t35
      df(1) = -0.2D1 * t38 * t86 * theta - 0.2D1 * t106 * t69 + 0.2D1
     $     * t4 * df(2)
      dfr0(10) = -t17 * t24 * ksnow
      t116 = dfr0(10) * ns
      dfr0(9) = -t116 * t55
      dfr0(8) = t116 * t37
      t118 = t116 * t36
      t124 = dfr0(8) * t36
      dfr0(7) = t118 * t59 - t124 * t66 + 0.2D1 * t32 * dfr0(9)
      t126 = t116 * t68
      dfr0(6) = t124 * t76 - 0.2D1 * t126 * t72
      dfr0(5) = -t17 * ksnow * t50
      dfr0(4) = -dfr0(5) * t83
      t134 = dfr0(4)
      dfr0(3) = t134 * rholiq * t89
      dfr0(2) = -t1 * t134 * t98 + t2 * dfr0(3) - t3 * dfr0(7)
      dfr0(1) = -0.2D1 * t106 * t126 + 0.2D1 * t4 * dfr0(2)
      grd(1,1) = -0.2D1 * t38 * t4 * t86 + t15 * t87 + df(6)
      grd(1,2) = -df(1)
      t153 = t32 * t6 * t35
      grd(1,3) = t57 * t153
      t154 = t36 * rholiq
      t157 = t2 * t30 * t71
      grd(1,4) = -0.2D1 * t51 * t154 * t157
      grd(2,1) = t134 * t15 + dfr0(6)
      grd(2,2) = -dfr0(1)
      grd(2,3) = t118 * t153
      grd(2,4) = -0.2D1 * t116 * t154 * t157
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
     $     ,Tfrz,thetar,ksnow,ns,eps
c-----------------------------------------------------------------------
C     Functions
      double precision smstep,dsmstep,inton,dinton
      BVAL(1) = ksnow*((U(1)-thetar)/(1-rholiq*omega**2*(Tfrz-U(2))**2
     $     /rhoice- thetar))**ns
     $     - 0.1d0*inton(T,3600.d0,3.d0*3600.d0,omega)
c$$$      BVAL(1) = U(1) - 0.5d0*inton(T,3600.d0,3.d0*3600.d0,5.d0)
C     derivative condition on temperature?
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
     $     ,Tfrz,thetar,ksnow,ns,eps
C-----------------------------------------------------------------------
C     Functions
      double precision smstep,dsmstep,inton,dinton
c$$$      DBDU(1,1) = ksnow*ns*(U(1)-thetar)**(ns-1)/(1-rholiq*omega**2
c$$$     $     *(Tfrz-U(2))**2/rhoice-thetar)**ns
      DBDU(1,1) = ksnow*ns/(1-rholiq*omega**2
     $     *(Tfrz-U(2))**2/rhoice-thetar)**ns
      DBDU(1,2) = 2.d0*ns*rholiq/rhoice*ksnow*(U(1)-thetar)**ns*(Tfrz
     $     -U(2))/(1-rholiq*omega**2*(Tfrz-U(2))**2/rhoice-thetar)**(ns
     $     +1)
c$$$      DBDU(1,1) = 1.D0
c$$$      DBDU(1,2) = 0.D0
c
      DBDUX(1,1) = 0.0D0
      DBDUX(1,2) = 0.0D0

      DBDU(2,1) = 0.d0
      DBDU(2,2) = 0.d0

      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = 1.0D0
c
      DBDT(1) = -0.1d0*dinton(T,3600.d0,3.d0*3600.d0,5.d0)
      DBDT(2) = 0.0D0
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
     $     ,Tfrz,thetar,ksnow,ns,eps
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
     $     ,Tfrz,thetar,ksnow,ns,eps

      U(1) = thetar
      U(2) = 273.15d0

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
     $     ,Tfrz,thetar,ksnow,ns,eps
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
     $     - dsmstep(b-x,omega)*smstep(a-x,omega)
      return
      end
