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
         FVAL(1) = -((rhoair * cair * (0.1D1 - rholiq * omega ** 2 *
     $        (Tfrz-U(2)) ** 2 / rhoice) + rholiq * cice * omega ** 2 *
     $        (Tfrz- U(2))** 2 + rholiq * cliq) * U(1) + rhoair * cair
     $        +0.2D1 * rholiq * Lfus * omega * (Tfrz - U(2)) * U(1))
     $        /((rhoair * cair * (0.1D1 - rholiq * omega ** 2 * (Tfrz
     $        -U(2)) ** 2 / rhoice) + rholiq * cice *omega ** 2 * (Tfrz
     $        -U(2)) ** 2 + rholiq * cliq) * U(1) + rhoair* cair)
     $        /(0.1D1 + omega ** 2 * (Tfrz - U(2)) ** 2) * ksnow * ns
     $        *((U(1) - thetar) / (0.1D1 - rholiq * omega ** 2 * (Tfrz
     $        -U(2)) ** 2 * U(1) - thetar)) ** (ns - 0.1D1) * ((0.1D1
     $        -rholiq * omega ** 2 * (Tfrz - U(2)) ** 2 * U(1) / rhoice
     $        -thetar) * UX(1) - (U(1)- thetar) * (0.2D1 * rholiq *
     $        omega** 2 * (Tfrz - U(2)) * UX(2)/ rhoice + (Tfrz - T) **
     $        2 *UX(1))) / (0.1D1 - rholiq * omega **2 * (Tfrz - U(2))
     $        ** 2* U(1) - thetar) ** 2
         FVAL(2) = -((rhoair * cair * (0.1D1 - rholiq * omega ** 2 *
     $        (Tfrz-U(2)) ** 2 / rhoice) + rholiq * cice * omega ** 2 *
     $        (Tfrz- U(2))** 2 + rholiq * cliq) * U(1) + rhoair * cair
     $        +0.2D1 * rholiq * Lfus * omega * (Tfrz - U(2)) * U(1))
     $        /((rhoair * cair * (0.1D1 - rholiq * omega ** 2 * (Tfrz
     $        -U(2)) ** 2 / rhoice) + rholiq * cice *omega ** 2 * (Tfrz
     $        -U(2)) ** 2 + rholiq * cliq) * U(1) + rhoair* cair)
     $        /(0.1D1 + omega ** 2 * (Tfrz - U(2)) ** 2) * ksnow * ns
     $        *((U(1) - thetar) / (0.1D1 - rholiq * omega ** 2 * (Tfrz
     $        -U(2)) ** 2 * U(1) - thetar)) ** (ns - 0.1D1) * ((0.1D1
     $        -rholiq * omega ** 2 * (Tfrz - U(2)) ** 2 * U(1) / rhoice
     $        -thetar) * UX(1) - (U(1)- thetar) * (0.2D1 * rholiq *
     $        omega** 2 * (Tfrz - U(2)) * UX(2)/ rhoice + (Tfrz - T) **
     $        2 *UX(1))) / (0.1D1 - rholiq * omega **2 * (Tfrz - U(2))
     $        ** 2* U(1) - thetar) ** 2
      else
         FVAL(1) = 0.d0
         FVAL(2) = 0.d0
      endif
c
c     And add diffusion
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
        doubleprecision df(11)
        doubleprecision dfr0(11)
        doubleprecision grd(2,4)
        doubleprecision t1
        doubleprecision t103
        doubleprecision t107
        doubleprecision t109
        doubleprecision t110
        doubleprecision t112
        doubleprecision t113
        doubleprecision t12
        doubleprecision t123
        doubleprecision t129
        doubleprecision t141
        doubleprecision t147
        doubleprecision t149
        doubleprecision t15
        doubleprecision t16
        doubleprecision t164
        doubleprecision t165
        doubleprecision t17
        doubleprecision t18
        doubleprecision t19
        doubleprecision t2
        doubleprecision t20
        doubleprecision t22
        doubleprecision t24
        doubleprecision t26
        doubleprecision t29
        doubleprecision t3
        doubleprecision t30
        doubleprecision t33
        doubleprecision t37
        doubleprecision t38
        doubleprecision t4
        doubleprecision t40
        doubleprecision t41
        doubleprecision t42
        doubleprecision t43
        doubleprecision t47
        doubleprecision t48
        doubleprecision t49
        doubleprecision t5
        doubleprecision t50
        doubleprecision t51
        doubleprecision t53
        doubleprecision t56
        doubleprecision t57
        doubleprecision t59
        doubleprecision t6
        doubleprecision t61
        doubleprecision t67
        doubleprecision t68
        doubleprecision t71
        doubleprecision t73
        doubleprecision t79
        doubleprecision t82
        doubleprecision t85
        doubleprecision t88
        doubleprecision t89
        doubleprecision t92
        doubleprecision t93
        doubleprecision t95
        doubleprecision t97
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
        t12 = t2 * t5
        t15 = t1 * (-t3 * t5 * t6 + 0.1D1) + rholiq * cice * t12 + rholi
     $       q * cliq
        t16 = t15 * theta
        t17 = rholiq * Lfus
        t18 = t16 + t1
        t24 = 0.1D1 / t18
        t30 = theta - thetar
        t19 = t3 * t5 * theta
        t33 = -t19 - thetar + 0.1D1
        t20 = 0.1D1 / t33
        t22 = ns - 0.1D1
        t37 = (t30 * t20) ** t22
        t48 = (Tfrz - T) ** 2
        t53 = t33 ** 2
        t26 = t19 * t6
        t29 = t3 * t4
        t38 = 0.2D1 * Tempz * t29 * t6 + t48 * thetaz
        t40 = (-t26 - thetar + 0.1D1) * thetaz - t30 * t38
        t41 = 0.1D1 / t53
        t42 = t40 * t41
        t56 = ns * t37 * t42
        t43 = t17 * omega
        t47 = 0.2D1 * t4 * t43 * theta + t1 + t16
        t49 = t47 * t24
        t50 = 0.1D1 + t12
        t51 = 0.1D1 / t50
        df(11) = -t49 * t51 * ksnow
        t57 = df(11) * ns
        t59 = t53 ** 2
        t61 = t37 * t40 / t59
        df(10) = -t57 * t61
        df(9) = t57 * t42
        t67 = df(9) * t37
        t68 = t22 * t20
        df(8) = 0.2D1 * t33 * df(10) - t67 * t68
        t71 = t3 * t6
        t73 = t71 * thetaz * t41
        df(7) = -t37 * t57 * t73 - t3 * df(8)
        t79 = -t37 * t38 * t41
        t82 = t22 / t30
        df(6) = t57 * t79 + t67 * t82
        t85 = ksnow * t56
        df(5) = -t47 * t51 * t85
        t88 = t18 ** 2
        t89 = 0.1D1 / t88
        t92 = t24 * t51 * t85
        df(4) = -t89 * df(5) - t92
        t93 = df(4)
        t95 = cice * theta
        t97 = t50 ** 2
        df(3) = t93 * rholiq * t95 + t49 / t97 * ksnow * t56
        t103 = t71 * theta
        t107 = df(7)
        df(2) = -t1 * t103 * t93 + t107 * theta + t2 * df(3)
        t109 = t37 * t30
        t110 = t57 * t109
        t112 = t6 * t41
        t113 = Tempz * t3 * t112
        df(1) = -0.2D1 * t43 * t92 * theta - 0.2D1 * t110 * t113 + 0.2D1
     $       * t4 * df(2)
        dfr0(11) = -t17 * t24 * ksnow
        t123 = dfr0(11) * ns
        dfr0(10) = -t123 * t61
        dfr0(9) = t123 * t42
        t129 = dfr0(9) * t37
        dfr0(8) = -t129 * t68 + 0.2D1 * t33 * dfr0(10)
        dfr0(7) = -t123 * t37 * t73 - t3 * dfr0(8)
        dfr0(6) = t123 * t79 + t129 * t82
        dfr0(5) = -t17 * ksnow * t56
        dfr0(4) = -dfr0(5) * t89
        t141 = dfr0(4)
        dfr0(3) = t141 * rholiq * t95
        t147 = dfr0(7)
        dfr0(2) = -t1 * t103 * t141 + t147 * theta + t2 * dfr0(3)
        t149 = t123 * t109
        dfr0(1) = -0.2D1 * t113 * t149 + 0.2D1 * t4 * dfr0(2)
        grd(1,1) = -0.2D1 * t4 * t43 * t92 + t107 * t5 + t15 * t93 +
     $       df(6)
        grd(1,2) = -df(1)
        t164 = t37 * (-t30 * t48 - t26 - thetar + 0.1D1) * t41
        grd(1,3) = t57 * t164
        t165 = t29 * t112
        grd(1,4) = -0.2D1 * t110 * t165
        grd(2,1) = t141 * t15 + t147 * t5 + dfr0(6)
        grd(2,2) = -dfr0(1)
        grd(2,3) = t123 * t164
        grd(2,4) = -0.2D1 * t149 * t165
        cgret(1,1) = grd(1,1)
        cgret(1,2) = grd(1,2)
        cgret(1,3) = grd(1,3)
        cgret(1,4) = grd(1,4)
        cgret(2,1) = grd(2,1)
        cgret(2,2) = grd(2,2)
        cgret(2,3) = grd(2,3)
        cgret(2,4) = grd(2,4)
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
c$$$      BVAL(1) = ksnow*(U(1)-thetar)/(1-rholiq*omega**2*(Tfrz-U(2))
c$$$     $     **2*U(1)/rhoice-thetar)
c$$$     $     - 0.1d0*inton(T, 30.d0, 3.d0 *3600.d0,50.d0)

      BVAL(1) = U(1)-thetar - (1.d0-rholiq*omega**2*(Tfrz-U(2))**2*U(1)
     $     /rhoice-thetar)*(0.1*inton(T,30.d0,3.d0*3600.d0,50.d0)/ksnow)
     $     **(1.d0/ns)
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
      BVAL(1) = U(1)
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
C-----------------------------------------------------------------------
C     Local variables
      double precision combo1,combo2
C-----------------------------------------------------------------------
      combo1 = rholiq*omega**2*(Tfrz-U(2))**2/rhoice
      combo2 = 0.1*inton(T,30.d0,3.d0*3600,50.d0)/ksnow

c$$$      DBDU(1,1) = ksnow*(*denom**ns + rholiq
c$$$     $     *omega**2*(Tfrz-U(2))**2*U(1)**(ns-1)*(U(1)-thetar))/denom
c$$$     $     **(2*n)

      DBDU(1,1) = 1.d0 + combo1*combo2**(1.d0/ns)

c$$$      DBDU(1,1) = 1.d0
c$$$      DBDU(1,2) = 0.d0
c$$$      DBDU(1,2) = 2.d0*ns*(U(1)-thetar)**ns*rholiq*omega**2*(Tfrz-U(2))
c$$$     $     *U(1)/rhoice/denom**(ns+1)
      DBDU(1,2) = -2.D0*rholiq*omega**2*(Tfrz-U(2))/rhoice*U(1)*combo2
     $     **(1.d0/ns)

c
      DBDUX(1,1) = 0.0D0
      DBDUX(1,2) = 0.0D0

      DBDU(2,1) = 0.d0
      DBDU(2,2) = 0.d0

      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = 1.0D0

c$$$      DBDT(1) = -0.1*dinton(T, 30.d0, 3.d0*3600.d0, 50.d0)

      DBDT(1) = - (1.d0 - combo1*U(1) - thetar)*combo2**(1.d0/ns-1.d0)
     $     *(0.1*dinton(T,30.d0,3.d0*3600.d0,50.d0)/ns/ksnow)

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
      DBDU(1,1) = 1.0D0
      DBDU(1,2) = 0.0D0

      DBDUX(1,1) = 0.0D0
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
