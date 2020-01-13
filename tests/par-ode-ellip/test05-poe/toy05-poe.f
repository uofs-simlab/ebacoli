c     Problem description in ./toy05_system.{tex,pdf}
c

C-----------------------------------------------------------------------
      SUBROUTINE F(T, X, U, UX, UXX, FVAL, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THIS SUBROUTINE DEFINES THE RIGHT HAND SIDE VECTOR OF THE
C     NPDE SYSTEM:
C     UT = F(T, X, U, V, UX, UXX).
C     VT = G(T, X, U, V)
C
C-----------------------------------------------------------------------
C     SUBROUTINE PARAMETERS:
C     INPUT:
      INTEGER                 NPDE
C     THE TOTAL NUMBER OF EQUATIONS IN THE SYSTEM.
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
C     VECTOR [F(T, X, U, V, UX, UXX)
C     G(T, X, U, V)         ] OF THE SYSTEM.
C-----------------------------------------------------------------------
c     Loop indices:
      integer                 i
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
c     The parabolic equations
      FVAL(1) = -X*U(1) + U(3) + UXX(1) - 2.D0*T*U(3) - (T*T-1.D0)*U(1)
     $     + (UX(5) + U(5)**2 -1.D0) + DACOS(U(2)) -2.D0*T - U(4)

      FVAL(2) = UX(2) + (UXX(2)+UX(6)) + (UX(1)-T*U(1)+U(3))
     $     + (DSQRT(1.D0-U(2)**2) - U(6))

C     The ODE equations
      FVAL(3) = -X*U(3) - U(1) + U(2) - DSQRT(1.D0-U(6)**2)

      FVAL(4) = -2.D0*DATANH(U(5)) + DSIN(DEXP(X*T)*U(1)-U(6))
     $     + (DEXP(X*T)*U(3)-U(2)) + U(4)**(1.5D0) - (X-T)**3

c     The elliptic constraints
      FVAL(5) = -UXX(5) - 2.D0*U(5)*UX(1)**2 + DSIN(U(4))
     $     + DEXP(2.D0*X*T)*(U(1)**2+U(3)**2 - 1.D0) - DSIN((X-T)**2)

      FVAL(6) = -UXX(6) + DEXP(X*T)*U(1) + (U(2)**2+U(6)**2)
     $     - (UX(6)**2+UX(2)**2) + (UX(1) - T*U(1) + U(3))
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
c     Loop indices:
      integer                 I, J
C-----------------------------------------------------------------------
C
C     ASSIGN DFDU(1:NPDE,1:NPDE), DFDUX(1:NPDE,1:NPDE), AND
C     DFDUXX(1:NPDE,1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
C     ZERO OUT ARRAYS
c$$$      DO J = 1, NPDE
c$$$         DO I = 1,NPDE
c$$$            DFDU(I,J) = 0.D0
c$$$            DFDUX(I,J) = 0.D0
c$$$            DFDUXX(I,J) = 0.D0
c$$$         END DO
c$$$      END DO

      DFDU(1,1) = -X - (T*T-1.D0)
      DFDU(1,2) = -1.D0/(1.D0-U(2)**2)
      DFDU(1,3) = 1.D0 - 2.D0*T
      DFDU(1,4) = -1.D0
      DFDU(1,5) = 2.D0*U(5)
      DFDU(1,6) = 0.D0

      DFDU(2,1) = -T
      DFDU(2,2) = -U(2)/(DSQRT(1-U(2)**2))
      DFDU(2,3) = 1.D0
      DFDU(2,4) = 0.D0
      DFDU(2,5) = 0.D0
      DFDU(2,6) = -1.D0

      DFDU(3,1) = -1.D0
      DFDU(3,2) = 1.D0
      DFDU(3,3) = -X
      DFDU(3,4) = 0.D0
      DFDU(3,5) = 0.D0
      DFDU(3,6) = -U(6)/(DSQRT(1.D0-U(6)**2))

      DFDU(4,1) = DEXP(X*T)*DCOS(DEXP(X*T)*U(1)-U(6))
      DFDU(4,2) = -1.D0
      DFDU(4,3) = DEXP(X*T)
      DFDU(4,4) = 1.5D0*DSQRT(U(4))
      DFDU(4,5) = -2.D0/(1.D0-U(5)**2)
      DFDU(4,6) = -DCOS(DEXP(X*T)*U(1)-U(6))

      DFDU(5,1) = 2.D0*DEXP(2.D0*X*T)*U(1)
      DFDU(5,2) = 0.D0
      DFDU(5,3) = 2.D0*DEXP(2.D0*X*T)*U(3)
      DFDU(5,4) = DCOS(U(4))
      DFDU(5,5) = -2.D0*UX(1)**2
      DFDU(5,6) = 0.D0

      DFDU(6,1) = DEXP(X*T) - T
      DFDU(6,2) = 2*U(2)
      DFDU(6,3) = 1.D0
      DFDU(6,4) = 0.D0
      DFDU(6,5) = 0.D0
      DFDU(6,6) = 2.D0*U(6)
c
      DFDUX(1,1) = 0.D0
      DFDUX(1,2) = 0.D0
      DFDUX(1,3) = 0.D0
      DFDUX(1,4) = 0.D0
      DFDUX(1,5) = 1.D0
      DFDUX(1,6) = 0.D0

      DFDUX(2,1) = 1.D0
      DFDUX(2,2) = 1.D0
      DFDUX(2,3) = 0.D0
      DFDUX(2,4) = 0.D0
      DFDUX(2,5) = 0.D0
      DFDUX(2,6) = 1.D0

      DFDUX(3,1) = 0.D0
      DFDUX(3,2) = 0.D0
      DFDUX(3,3) = 0.D0
      DFDUX(3,4) = 0.D0
      DFDUX(3,5) = 0.D0
      DFDUX(3,6) = 0.D0

      DFDUX(4,1) = 0.D0
      DFDUX(4,2) = 0.D0
      DFDUX(4,3) = 0.D0
      DFDUX(4,4) = 0.D0
      DFDUX(4,5) = 0.D0
      DFDUX(4,6) = 0.D0

      DFDUX(5,1) = -4.D0*U(5)*UX(1)
      DFDUX(5,2) = 0.D0
      DFDUX(5,3) = 0.D0
      DFDUX(5,4) = 0.D0
      DFDUX(5,5) = 0.D0
      DFDUX(5,6) = 0.D0

      DFDUX(6,1) = 1.D0
      DFDUX(6,2) = -2.D0*UX(2)
      DFDUX(6,3) = 0.D0
      DFDUX(6,4) = 0.D0
      DFDUX(6,5) = 0.D0
      DFDUX(6,6) = -2.D0*UX(6)
c
      DFDUXX(1,1) = 1.D0
      DFDUXX(1,2) = 0.D0
      DFDUXX(1,3) = 0.D0
      DFDUXX(1,4) = 0.D0
      DFDUXX(1,5) = 0.D0
      DFDUXX(1,6) = 0.D0

      DFDUXX(2,1) = 0.D0
      DFDUXX(2,2) = 1.D0
      DFDUXX(2,3) = 0.D0
      DFDUXX(2,4) = 0.D0
      DFDUXX(2,5) = 0.D0
      DFDUXX(2,6) = 0.D0

      DFDUXX(3,1) = 0.D0
      DFDUXX(3,2) = 0.D0
      DFDUXX(3,3) = 0.D0
      DFDUXX(3,4) = 0.D0
      DFDUXX(3,5) = 0.D0
      DFDUXX(3,6) = 0.D0

      DFDUXX(4,1) = 0.D0
      DFDUXX(4,2) = 0.D0
      DFDUXX(4,3) = 0.D0
      DFDUXX(4,4) = 0.D0
      DFDUXX(4,5) = 0.D0
      DFDUXX(4,6) = 0.D0

      DFDUXX(5,1) = 0.D0
      DFDUXX(5,2) = 0.D0
      DFDUXX(5,3) = 0.D0
      DFDUXX(5,4) = 0.D0
      DFDUXX(5,5) = -1.D0
      DFDUXX(5,6) = 0.D0

      DFDUXX(6,1) = 0.D0
      DFDUXX(6,2) = 0.D0
      DFDUXX(6,3) = 0.D0
      DFDUXX(6,4) = 0.D0
      DFDUXX(6,5) = 0.D0
      DFDUXX(6,6) = -1.D0

      DFDUXX(2,1) = 0.D0
      DFDUXX(2,2) = 0.D0
      DFDUXX(2,3) = 0.D0

      DFDUXX(3,1) = 0.D0
      DFDUXX(3,2) = 0.D0
      DFDUXX(3,3) = 1.D0
C
      RETURN
      END
      SUBROUTINE BNDXA(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C     LEFT SPATIAL END POINT X = XA.
C     B(T, U, V, UX) = 0
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
c     Loop indices:
      integer                 i
C-----------------------------------------------------------------------
C
C     DO NOT SET ODE COMPONENTS OF BVAL
C
      BVAL(1) = U(1) - DEXP(-T)*DSIN(1.D0+T)
      BVAL(2) = U(2) - DCOS(1.D0+T)
      BVAL(5) = U(5) - DTANH(1.D0-T)
      BVAL(6) = U(6) - DSIN(1.D0+T)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE BNDXB(T, U, UX, BVAL, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THE SUBROUTINE IS USED TO DEFINE THE BOUNDARY CONDITIONS AT THE
C     RIGHT SPATIAL END POINT X = XB.
C     B(T, U, V, UX) = 0
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
c     Loop indices:
      integer                 i
C-----------------------------------------------------------------------
c     DO NOT SET ODE COMPONENTS OF BVAL
      BVAL(1) = -UX(1) - T*U(1) + U(3)
      BVAL(2) = -UX(2) + DSQRT(1.D0+UX(6)**2)
      BVAL(5) = -UX(5) + 1.D0 - U(5)**2 + U(4) - (2.D0-T)**2
      BVAL(6) = -UX(6) + U(2) + UX(2) + U(6)
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DIFBXA(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
C-----------------------------------------------------------------------
C     PURPOSE:
C     THE SUBROUTINE IS USED TO DEFINE THE DIFFERENTIATED BOUNDARY
C     CONDITIONS AT THE LEFT SPATIAL END POINT X = XA. FOR THE
C     BOUNDARY CONDITION EQUATION
C     B(T, U, V, UX) = 0
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
c     Loop indices:
      integer                 i, j
C-----------------------------------------------------------------------
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF
C     U(1:NPDE), UX(1:NPDE),
C

c     DO NOT SET ODE ROWS OF DERIVATIVE
      DBDU(1,1) = 1.0D0
      DBDU(1,2) = 0.0D0
      DBDU(1,3) = 0.0D0
      DBDU(1,4) = 0.0D0
      DBDU(1,5) = 0.0D0
      DBDU(1,6) = 0.0D0

      DBDU(2,1) = 0.0D0
      DBDU(2,2) = 1.0D0
      DBDU(2,3) = 0.0D0
      DBDU(2,4) = 0.0D0
      DBDU(2,5) = 0.0D0
      DBDU(2,6) = 0.0D0

      DBDU(5,1) = 0.0D0
      DBDU(5,2) = 0.0D0
      DBDU(5,3) = 0.0D0
      DBDU(5,4) = 0.0D0
      DBDU(5,5) = 1.0D0
      DBDU(5,6) = 0.0D0

      DBDU(6,1) = 0.0D0
      DBDU(6,2) = 0.0D0
      DBDU(6,3) = 0.0D0
      DBDU(6,4) = 0.0D0
      DBDU(6,5) = 0.0D0
      DBDU(6,6) = 1.0D0

      DBDUX(1,1) = 0.0D0
      DBDUX(1,2) = 0.0D0
      DBDUX(1,3) = 0.0D0
      DBDUX(1,4) = 0.0D0
      DBDUX(1,5) = 0.0D0
      DBDUX(1,6) = 0.0D0

      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = 0.0D0
      DBDUX(2,3) = 0.0D0
      DBDUX(2,4) = 0.0D0
      DBDUX(2,5) = 0.0D0
      DBDUX(2,6) = 0.0D0

      DBDUX(5,1) = 0.0D0
      DBDUX(5,2) = 0.0D0
      DBDUX(5,3) = 0.0D0
      DBDUX(5,4) = 0.0D0
      DBDUX(5,5) = 0.0D0
      DBDUX(5,6) = 0.0D0

      DBDUX(6,1) = 0.0D0
      DBDUX(6,2) = 0.0D0
      DBDUX(6,3) = 0.0D0
      DBDUX(6,4) = 0.0D0
      DBDUX(6,5) = 0.0D0
      DBDUX(6,6) = 0.0D0

      DBDT(1) = DEXP(-T)*DSIN(1.D0+T)+DEXP(-T)*DCOS(T)
      DBDT(2) = DSIN(1.D0+T)
      DBDT(5) = 1.D0/DCOSH(1.D0-T)**2
      DBDT(6) = -DCOS(1.D0+T)
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
      DOUBLE PRECISION        DBDT(NPDE)
C     DBDT(I) IS THE PARTIAL DERIVATIVE
C     OF THE I-TH COMPONENT OF THE VECTOR B
C     WITH RESPECT TO TIME T.
c-----------------------------------------------------------------------
c     Loop indices:
      integer                 i, j
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF
C     U(1:NPDE), UX(1:NPDE)
C
c     DO NOT SET ODE ROWS OF DERIVATIVE
      DBDU(1,1) = -T
      DBDU(1,2) = 0.0D0
      DBDU(1,3) = 1.0D0
      DBDU(1,4) = 0.0D0
      DBDU(1,5) = 0.0D0
      DBDU(1,6) = 0.0D0

      DBDU(2,1) = 0.0D0
      DBDU(2,2) = 1.0D0
      DBDU(2,3) = 0.0D0
      DBDU(2,4) = 0.0D0
      DBDU(2,5) = 0.0D0
      DBDU(2,6) = 0.0D0

      DBDU(5,1) = 0.0D0
      DBDU(5,2) = 0.0D0
      DBDU(5,3) = 0.0D0
      DBDU(5,4) = 1.0D0
      DBDU(5,5) = -2.0D0*U(5)
      DBDU(5,6) = 0.0D0

      DBDU(6,1) = 0.0D0
      DBDU(6,2) = 1.0D0
      DBDU(6,3) = 0.0D0
      DBDU(6,4) = 0.0D0
      DBDU(6,5) = 0.0D0
      DBDU(6,6) = 1.0D0

      DBDUX(1,1) = -1.0D0
      DBDUX(1,2) = 0.0D0
      DBDUX(1,3) = 0.0D0
      DBDUX(1,4) = 0.0D0
      DBDUX(1,5) = 0.0D0
      DBDUX(1,6) = 0.0D0

      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = -1.0D0
      DBDUX(2,3) = 0.0D0
      DBDUX(2,4) = 0.0D0
      DBDUX(2,5) = 0.0D0
      DBDUX(2,6) = UX(6)/(DSQRT(1.D0+UX(6)**2))

      DBDUX(5,1) = 0.0D0
      DBDUX(5,2) = 0.0D0
      DBDUX(5,3) = 0.0D0
      DBDUX(5,4) = 0.0D0
      DBDUX(5,5) = -1.0D0
      DBDUX(5,6) = 0.0D0

      DBDUX(6,1) = 0.0D0
      DBDUX(6,2) = 1.0D0
      DBDUX(6,3) = 0.0D0
      DBDUX(6,4) = 0.0D0
      DBDUX(6,5) = 0.0D0
      DBDUX(6,6) = -1.0D0

      DBDT(1) = U(1)
      DBDT(2) = 0.D0
      DBDT(5) = 2.D0*(2.D0-T)
      DBDT(6) = 0.D0
      RETURN
      END
C-----------------------------------------------------------------------
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
c-----------------------------------------------------------------------
c     Loop indices:
      integer                 i
C-----------------------------------------------------------------------
C
C     ASSIGN U(1:NPDE) THE INITIAL VALUES OF U(T0,X).
      U(1) = DSIN(X)
      U(2) = DCOS(X)
      U(3) = DCOS(X)
      U(4) = X*X
      U(5) = DTANH(X)
      U(6) = DSIN(X)
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
      U(1) = DEXP(-X*T)*DSIN(X+T)
      U(2) = DCOS(X+T)
      U(3) = DEXP(-X*T)*DCOS(X+T)
      U(4) = (X-T)**2
      U(5) = DTANH(X-T)
      U(6) = DSIN(X+T)
C
      RETURN
      END
