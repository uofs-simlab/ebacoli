c     This is the problem defintion for the model of ecological
c     competition in the presence of water-ice solidification

C-----------------------------------------------------------------------
      SUBROUTINE F(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE
C       NPDE DIMENSIONAL PARABOLIC PARTIAL DIFFERENTIAL EQUATION
C                        UT = F(T, X, U, UX, UXX).
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
        DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        FVAL(NPDE)
C                               FVAL(1:NPDE) IS THE RIGHT HAND SIDE
C                               VECTOR F(T, X, U, UX, UXX) OF THE PDE.
c-----------------------------------------------------------------------
c parameters
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
c-----------------------------------------------------------------------
c functions
      double precision volcon, dvolco, ddvolc,
     &     hydcon, dhydco, ddhydc, qs
c
c-----------------------------------------------------------------------
c Local variables
      double precision psi, dpsi, ddpsi, theta, dtheta, k, dk
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      psi = U(1)
      dpsi = UX(1)
      ddpsi = UXX(1)
C
      theta = volcon(psi)
      dtheta = dvolco(psi)
C
      k = hydcon(psi)
      dk = dhydco(psi)
C
      FVAL(1) = ( k*ddpsi + dk*dpsi*(-1.d0+dpsi) + qs(x,t) )
     &     / (ss*theta/phi+dtheta)
C
      RETURN
      END
      SUBROUTINE DERIVF(T, X, U, UX, UXX, DFDU, DFDUX, DFDUXX, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE IS USED TO DEFINE THE INFORMATION ABOUT THE
C       PDE REQUIRED TO FORM THE ANALYTIC JACOBIAN MATRIX FOR THE DAE
C       OR ODE SYSTEM. ASSUMING THE PDE IS OF THE FORM
C                        UT = F(T, X, U, UX, UXX)
C       THIS ROUTINE RETURNS THE JACOBIANS D(F)/D(U), D(F)/D(UX), AND
C       D(F)/D(UXX).
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
        DOUBLE PRECISION        UXX(NPDE)
C                               UXX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SECOND SPATIAL DERIVATIVE OF THE
C                               SOLUTION AT THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        DFDU(NPDE,NPDE)
C                               DFDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DFDUX(NPDE,NPDE)
C                               DFDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE
C                               UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DFDUXX(NPDE,NPDE)
C                               DFDUXX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR F
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SECOND SPATIAL DERIVATIVE OF THE
C                               UNKNOWN FUNCTION U.
c-----------------------------------------------------------------------
c parameters
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
c
c-----------------------------------------------------------------------
c functions
      double precision volcon, dvolco, ddvolc,
     &     hydcon, dhydco, ddhydc, qs
c-----------------------------------------------------------------------
c Local variables
      double precision psi, dpsi, ddpsi, theta, dtheta, ddthet, k, dk,
     $     ddk, denom1
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      psi = U(1)
      dpsi = UX(1)
      ddpsi = UXX(1)
C
      theta = volcon(psi)
      dtheta = dvolco(psi)
      ddthet = ddvolc(psi)
C
      k = hydcon(psi)
      dk = dhydco(psi)
      ddk = ddhydc(psi)
c
      denom1 = ss*theta/phi+dtheta
C
      DFDU(1,1) = -(ss*dtheta/phi+ddthet)*
     &     (k*ddpsi+dk*dpsi*(-1.d0+dpsi)+qs(x,t)) / denom1**2
     &     + (dk*ddpsi+ddk*dpsi*(-1.d0+dpsi))/denom1
c
      DFDUX(1,1) = dk*(-1.d0+2.d0*dpsi)/denom1
c
      DFDUXX(1,1) = k/denom1
C
      RETURN
      END
      SUBROUTINE BNDXA(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C       LEFT SPATIAL END POINT X = XA.
C                           B(T, U, UX) = 0
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,XA).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,XA).
C
C OUTPUT:
        DOUBLE PRECISION        BVAL(NPDE)
C                               BVAL(1:NPDE) IS THE BOUNDARY CONTIDITION
C                               AT THE LEFT BOUNDARY POINT.
c-----------------------------------------------------------------------
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      BVAL(1) = U(1) + 0.75d0
C
      RETURN
      END
      SUBROUTINE BNDXB(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C       RIGHT SPATIAL END POINT X = XB.
C                           B(T, U, UX) = 0
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,XB).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,XB).
C
C OUTPUT:
        DOUBLE PRECISION        BVAL(NPDE)
C                               BVAL(1:NPDE) IS THE BOUNDARY CONTIDITION
C                               AT THE RIGHT BOUNDARY POINT.
c-----------------------------------------------------------------------
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
C     Careful with units of this!
      BVAL(1) = U(1) + 10.0d0
C
      RETURN
      END
      SUBROUTINE DIFBXA(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY
C       CONDITIONS AT THE LEFT SPATIAL END POINT X = XA. FOR THE
C       BOUNDARY CONDITION EQUATION
C                              B(T, U, UX) = 0
C       THE PARTIAL DERIVATIVES DB/DU, DB/DUX, AND DB/DT ARE SUPPLIED
C       BY THIS ROUTINE.
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE
C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).
C
C OUTPUT:
        DOUBLE PRECISION        DBDU(NPDE,NPDE)
C                               DBDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDUX(NPDE,NPDE)
C                               DBDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE
C                               UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDT(NPDE)
C                               DBDT(I) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO TIME T.
C-----------------------------------------------------------------------
      DBDU(1,1) = 1.0D0
c
      DBDUX(1,1) = 0.0D0
c
      DBDT(1) = 0.0D0
      RETURN
      END
      SUBROUTINE DIFBXB(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY
C       CONDITIONS AT THE RIGHT SPATIAL END POINT 1 = XB. FOR THE
C       BOUNDARY CONDITION EQUATION
C                              B(T, U, UX) = 0
C       THE PARTIAL DERIVATIVES DB/DU, DB/DUX, AND DB/DT ARE SUPPLIED
C       BY THIS ROUTINE.
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
C
        DOUBLE PRECISION        UX(NPDE)
C                               UX(1:NPDE) IS THE APPROXIMATION OF THE

C                               SPATIAL DERIVATIVE OF THE SOLUTION AT
C                               THE POINT (T,X).

C
C OUTPUT:
        DOUBLE PRECISION        DBDU(NPDE,NPDE)
C                               DBDU(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDUX(NPDE,NPDE)
C                               DBDUX(I,J) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO THE J-TH COMPONENT
C                               OF THE SPATIAL DERIVATIVE OF THE
C                               UNKNOWN FUNCTION U.
C
        DOUBLE PRECISION        DBDT(NPDE),expp,exppp
C                               DBDT(I) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO TIME T.
C-----------------------------------------------------------------------
      DBDU(1,1) = 1.0D0
c
      DBDUX(1,1) = 0.0D0
c
      DBDT(1) = 0.0D0
C
      RETURN
      END
      SUBROUTINE UINIT(X, U, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE IS USED TO RETURN THE NPDE-VECTOR OF INITIAL
C       CONDITIONS OF THE UNKNOWN FUNCTION AT THE INITIAL TIME T = T0
C       AT THE SPATIAL COORDINATE X.
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      DOUBLE PRECISION        X
C                               THE SPATIAL COORDINATE.
C
      INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
C OUTPUT:
      DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS VECTOR OF INITIAL VALUES OF
C                               THE UNKNOWN FUNCTION AT T = T0 AND THE
C                               GIVEN VALUE OF X.
C-----------------------------------------------------------------------
      if (X .le. 0.6) then
         U(1) = -0.75 - 9.25d0*X/0.6
      else
         U(1) = -10.d0
      end if
C

      RETURN
      END
      SUBROUTINE TRUU(T, X, U, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C     THIS FUNCTION PROVIDES THE EXACT SOLUTION OF THE PDE.
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE NUMBER OF PDES IN THE SYSTEM.
C
        DOUBLE PRECISION        T
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
C OUTPUT:
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE EXACT SOLUTION AT THE
C                               POINT (T,X).
C-----------------------------------------------------------------------
c Loop indices:
        integer                 i
C-----------------------------------------------------------------------
        U(1) = 0.d0
C
      RETURN
      END
      doubleprecision function qs (x,t)
C     Source/sink term
      double precision x, t
      qs = 0.d0
      return
      end
      doubleprecision function hydcon (psi)
      doubleprecision psi
C     Hydraulic conductivity K(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      hydcon = ks * (0.1D1 - (alpha * psi) ** (n - 0.1D1) * (0.1D1 +
     $     (alpha * psi) ** n) ** (-m)) ** 2 / (0.1D1 + (alpha * psi) **
     $     n) ** (m / 0.2D1)
      return
      end
      doubleprecision function dhydco (psi)
      doubleprecision psi
C     Derivative of hydraulic conductivity K'(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      dhydco = 0.2D1 * ks * (0.1D1 - (alpha * psi) ** (n - 0.1D1) * (0
     $     .1D1 + (alpha * psi) ** n) ** (-m)) / (0.1D1 + (alpha *
     $     psi) ** n) ** (m / 0.2D1) * (-(alpha * psi) ** (n - 0.1D1)
     $     * (n - 0.1D1) / p si * (0.1D1 + (alpha * psi) ** n) ** (-m)
     $     + (alpha * psi) ** (n - 0.1D1) * (0.1D1 + (alpha * psi) **
     $     n) ** (-m) * m * (alpha * psi) ** n * n / psi / (0.1D1 +
     $     (alpha * psi) ** n)) - ks * (0.1D1 - (al pha * psi) ** (n -
     $     0.1D1) * (0.1D1 + (alpha * psi) ** n) ** (-m)) ** 2 /
     $     (0.1D1 + (alpha * psi) ** n) ** (m / 0.2D1) * m * (alpha *
     $     psi) ** n * n / psi / (0.1D1 + (alpha * psi) ** n) / 0.2D1
      return
      end
      doubleprecision function ddhydc (psi)
      doubleprecision psi
      doubleprecision t1
      doubleprecision t11
      doubleprecision t12
      doubleprecision t13
      doubleprecision t14
      doubleprecision t18
      doubleprecision t19
      doubleprecision t2
      doubleprecision t22
      doubleprecision t23
      doubleprecision t26
      doubleprecision t27
      doubleprecision t3
      doubleprecision t30
      doubleprecision t36
      doubleprecision t38
      doubleprecision t39
      doubleprecision t4
      doubleprecision t40
      doubleprecision t42
      doubleprecision t47
      doubleprecision t49
      doubleprecision t5
      doubleprecision t50
      doubleprecision t52
      doubleprecision t55
      doubleprecision t58
      doubleprecision t59
      doubleprecision t6
      doubleprecision t61
      doubleprecision t68
      doubleprecision t69
      doubleprecision t7
      doubleprecision t75
      doubleprecision t8
C     2nd derivative of hydraulic conductivity K''(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      t1 = alpha * psi
      t2 = n - 0.1D1
      t3 = t1 ** t2
      t4 = t3 * t2
      t5 = 0.1D1 / psi
      t6 = t1 ** n
      t7 = 0.1D1 + t6
      t8 = t7 ** (-m)
      t11 = t3 * t8
      t12 = t11 * m
      t13 = t6 * n
      t14 = 0.1D1 / t7
      t18 = t12 * t13 * t14 * t5 - t4 * t5 * t8
      t19 = t18 ** 2
      t22 = t7 ** (m / 0.2D1)
      t23 = 0.1D1 / t22
      t26 = 0.1D1 - t11
      t27 = ks * t26
      t30 = m * t6
      t36 = t2 ** 2
      t38 = psi ** 2
      t39 = 0.1D1 / t38
      t40 = t39 * t8
      t42 = t4 * t40
      t47 = m ** 2
      t49 = t6 ** 2
      t50 = n ** 2
      t52 = t7 ** 2
      t55 = t49 * t50 * t39 / t52
      t58 = t39 * t14
      t59 = t6 * t50 * t58
      t61 = t13 * t58
      t68 = t26 ** 2
      t69 = ks * t68
      t75 = t69 * t23 * m
      ddhydc = 0.2D1 * ks * t19 * t23 - 0.2D1 * t27 * t23 * t18 * t30
     $     * n * t5 * t14 + 0.2D1 * t27 * t23 * (0.2D1 * n * t14 * t30
     $     * t42 - t11 * t47 * t55 - t3 * t36 * t40 - t12 * t55 + t12
     $     * t59 - t12 * t61 + t42) + t69 * t23 * t47 * t55 / 0.4D1 -
     $     t75 * t59 / 0.2D1 + t75 * t61 / 0.2D1 + t75 * t55 / 0.2D1
      return
      end
      doubleprecision function volcon (psi)
      doubleprecision psi
C     Volumetric water content theta(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      volcon = (thetas - thetar) / (0.1D1 + (alpha * psi) ** n) ** m +
     &     thetar
        return
      end
      doubleprecision function dvolco (psi)
      doubleprecision psi
C     Derivative of volumetric water content theta'(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      dvolco = -(thetas - thetar) / (0.1D1 + (alpha * psi) ** n) ** m *
     $     m * (alpha * psi) ** n * n / psi / (0.1D1 + (alpha * psi) **
     $     n)
        return
      end
      doubleprecision function ddvolc (psi)
      doubleprecision psi
      doubleprecision t10
      doubleprecision t11
      doubleprecision t13
      doubleprecision t14
      doubleprecision t15
      doubleprecision t18
      doubleprecision t20
      doubleprecision t23
      doubleprecision t3
      doubleprecision t4
      doubleprecision t5
      doubleprecision t7
      doubleprecision t8
C     2nd derivative of Volumetric water content theta''(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
        t3 = (alpha * psi) ** n
        t4 = 0.1D1 + t3
        t5 = t4 ** m
        t7 = (thetas - thetar) / t5
        t8 = m ** 2
        t10 = t3 ** 2
        t11 = n ** 2
        t13 = psi ** 2
        t14 = 0.1D1 / t13
        t15 = t4 ** 2
        t18 = t10 * t11 * t14 / t15
        t20 = t7 * m
        t23 = t14 / t4
        ddvolc = n * t20 * t23 * t3 - t11 * t20 * t23 * t3 + t18 * t7 *
     #t8 + t20 * t18
        return
      end
