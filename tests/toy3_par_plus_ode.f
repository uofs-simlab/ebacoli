c     This is the problem defintion for the toy problem:
c
c     u1_t = x u1_xx + t x u2_x + (1 + t x) e^{-tx} u3
c     u2_t = x u2_xx + t x u1_x + (1 + t x) e^{-tx} u4
c     u3_t = -0.5 e^{tx} u1 - 0.5 \sqrt{1-u4^2}
c     u4_t = -0.5 e^{tx} u2 - 0.5 \sqrt{1-u3^2}
c
c     1 <= x <= 2,   0 <= t <= 1
c
c     BC:
c                  u1(1,t) - e^{-t} \sin{1+t} = 0
c     u1_x(2,t) + t u1(2,t) - e^{-2t} u3(2,t) = 0
c                  u2(1,t) - e^{-t} \sin{1+t} = 0
c     u2_x(2,t) + t u2(2,t) - e^{-2t} u4(2,t) = 0
c
c     IC:
c     u1(x,0) = \sin{x}
c     u2(x,0) = \sin{x}
c     u3(x,0) = \cos{x}
c     u4(x,0) = \cos{x}
c
c     Exact solu1tion:
c     u1(x,t) = e^{-tx} \sin{x+t}
c     u2(x,t) = e^{-tx} \sin{x+t}
c     u3(x,t) = \cos{x+t}
c     u4(x,t) = \cos{x+t}

C-----------------------------------------------------------------------
      SUBROUTINE F(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE
C       NPDE SYSTEM:
C                        UT = F(T, X, U, V, UX, UXX).
C                        VT = G(T, X, U, V)
C
C     Note that U is the first NU elements of the U input vector
C          (Parabolic subsystem) and V are the remaining NPDE-NU
C          elements.
C
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
        INTEGER                 NPDE
C                               THE TOTAL NUMBER OF EQUATIONS IN THE SYSTEM.
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
C                               VECTOR [F(T, X, U, V, UX, UXX)
C                                       G(T, X, U, V)         ] OF THE SYSTEM.
C-----------------------------------------------------------------------
c Loop indices:
      integer                 i
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
c     The parabolic equations
      FVAL(1)= X*UXX(1) + T*X*UX(2) + (1+T*X)*DEXP(-T*X)*U(3)
      FVAL(2)= X*UXX(2) + T*X*UX(1) + (1+T*X)*DEXP(-T*X)*U(4)

c     The ODEs
      FVAL(3) = -0.5D0*(DEXP(T*X)*U(1) + DSQRT(1-U(4)*U(4)))
      FVAL(4) = -0.5D0*(DEXP(T*X)*U(2) + DSQRT(1-U(3)*U(3)))
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
c Loop indices:
        integer                 i, j
C-----------------------------------------------------------------------
C
C     ASSIGN DFDU(1:NPDE,1:NPDE), DFDUX(1:NPDE,1:NPDE), AND
C     DFDUXX(1:NPDE,1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      DFDU(1,1) = 0.D0
      DFDU(1,2) = 0.D0
      DFDU(1,3) = (1+T*X)*DEXP(-T*X)
      DFDU(1,4) = 0.D0

      DFDU(2,1) = 0.D0
      DFDU(2,2) = 0.D0
      DFDU(2,3) = 0.D0
      DFDU(2,4) = (1+T*X)*DEXP(-T*X)

      DFDU(3,1) = -0.5D0*DEXP(T*X)
      DFDU(3,2) = 0.D0
      DFDU(3,3) = 0.D0
      DFDU(3,4) = -0.5D0*DSQRT(1-U(4)*U(4))*U(4)

      DFDU(4,1) = 0.D0
      DFDU(4,2) = -0.5D0*DEXP(T*X)
      DFDU(4,3) = -0.5D0*DSQRT(1-U(3)*U(3))*U(3)
      DFDU(4,4) = 0.D0
c
      DFDUX(1,1) = 0.D0
      DFDUX(1,2) = T*X
      DFDUX(1,3) = 0.D0
      DFDUX(1,4) = 0.D0

      DFDUX(2,1) = T*X
      DFDUX(2,2) = 0.D0
      DFDUX(2,3) = 0.D0
      DFDUX(2,4) = 0.D0

      DFDUX(3,1) = 0.D0
      DFDUX(3,2) = 0.D0
      DFDUX(3,3) = 0.D0
      DFDUX(3,4) = 0.D0

      DFDUX(4,1) = 0.D0
      DFDUX(4,2) = 0.D0
      DFDUX(4,3) = 0.D0
      DFDUX(4,4) = 0.D0
c
      DFDUXX(1,1) = X
      DFDUXX(1,2) = 0.d0
      DFDUXX(1,3) = 0.d0
      DFDUXX(1,4) = 0.d0

      DFDUXX(2,1) = 0.d0
      DFDUXX(2,2) = X
      DFDUXX(2,3) = 0.d0
      DFDUXX(2,4) = 0.d0

      DFDUXX(3,1) = 0.d0
      DFDUXX(3,2) = 0.D0
      DFDUXX(3,3) = 0.d0
      DFDUXX(3,4) = 0.d0

      DFDUXX(4,1) = 0.d0
      DFDUXX(4,2) = 0.D0
      DFDUXX(4,3) = 0.d0
      DFDUXX(4,4) = 0.d0
C
      RETURN
      END
      SUBROUTINE BNDXA(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C       LEFT SPATIAL END POINT X = XA.
C                           B(T, U, V, UX) = 0
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
c Loop indices:
        integer                 i
C-----------------------------------------------------------------------
      BVAL(1) = U(1) - DEXP(-T)*DSIN(1+T)
      BVAL(2) = U(2) - DEXP(-T)*DSIN(1+T)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE BNDXB(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C       RIGHT SPATIAL END POINT X = XB.
C                           B(T, U, V, UX) = 0
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
c Loop indices:
        integer                 i
C-----------------------------------------------------------------------
      BVAL(1) = UX(1) - T*DEXP(-2.D0*T)*UX(3) - DEXP(-2.D0*T)*U(3)
      BVAL(2) = UX(2) - T*DEXP(-2.D0*T)*UX(4) - DEXP(-2.D0*T)*U(4)
C
      RETURN
      END
      SUBROUTINE DIFBXA(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C PURPOSE:
C       THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY
C       CONDITIONS AT THE LEFT SPATIAL END POINT X = XA. FOR THE
C       BOUNDARY CONDITION EQUATION
C                              B(T, U, V, UX) = 0
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
c Loop indices:
        integer                 i, j
C-----------------------------------------------------------------------
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF
C     U(1:NPDE), UX(1:NPDE),
C
      DBDU(1,1) = 1.0D0
      DBDU(1,2) = 0.0D0
      DBDU(1,3) = 0.0D0
      DBDU(1,4) = 0.0D0

      DBDU(2,1) = 0.0D0
      DBDU(2,2) = 1.0D0
      DBDU(2,3) = 0.0D0
      DBDU(2,4) = 0.0D0

      DBDUX(1,1) = 0.0D0
      DBDUX(1,2) = 0.0D0
      DBDUX(1,3) = 0.0D0
      DBDUX(1,4) = 0.0D0

      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = 0.0D0
      DBDUX(2,3) = 0.0D0
      DBDUX(2,4) = 0.0D0

      DBDT(1) = -DEXP(-T)*(DCOS(1+T)-DSIN(1+T))
      DBDT(2) = -DEXP(-T)*(DCOS(1+T)-DSIN(1+T))
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
        DOUBLE PRECISION        DBDT(NPDE)
C                               DBDT(I) IS THE PARTIAL DERIVATIVE
C                               OF THE I-TH COMPONENT OF THE VECTOR B
C                               WITH RESPECT TO TIME T.
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, j
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF
C     U(1:NPDE), UX(1:NPDE)
C
      DBDU(1,1) = 0.D0
      DBDU(1,2) = 0.D0
      DBDU(1,3) = -DEXP(-2*T)
      DBDU(1,4) = 0.D0

      DBDU(2,1) = 0.D0
      DBDU(2,2) = 0.D0
      DBDU(2,3) = 0.D0
      DBDU(2,4) = -DEXP(-2*T)

      DBDUX(1,1) = 1.0D0
      DBDUX(1,2) = 0.0D0
      DBDUX(1,3) = -T*DEXP(-2*T)
      DBDUX(1,4) = 0.0D0

      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = 1.0D0
      DBDUX(2,3) = 0.0D0
      DBDUX(2,4) = -T*DEXP(-2*T)

      DBDT(1)=(2.D0*T-1.D0)*DEXP(-2.D0*T)*UX(3)+2.D0*DEXP(-2.D0*T)*U(3)
      DBDT(2)=(2.D0*T-1.D0)*DEXP(-2.D0*T)*UX(4)+2.D0*DEXP(-2.D0*T)*U(4)

      RETURN
      END
C-----------------------------------------------------------------------
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
c-----------------------------------------------------------------------
c Loop indices:
      integer                 i
C-----------------------------------------------------------------------
C
C     ASSIGN U(1:NPDE) THE INITIAL VALUES OF U(T0,X).
      U(1) = DSIN(X)
      U(2) = DSIN(X)
      U(3) = DCOS(X)
      U(4) = DCOS(X)
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
      U(1) = DEXP(-T*X)*DSIN(X+T)
      U(2) = DEXP(-T*X)*DSIN(X+T)
      U(3) = DCOS(X+T)
      U(4) = DCOS(X+T)
C
      RETURN
      END
