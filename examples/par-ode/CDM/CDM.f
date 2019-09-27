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
        double precision    d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
        common /ecosolid/   d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
        double precision    FUNC1,FUNC2
c
c-----------------------------------------------------------------------
c Loop indices:
      integer                 i
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
c     The parabolic equation
      FVAL(1)=d1*UXX(1) + FUNC1(U(1)) -
     &     U(1)*(s1*U(2)+k1*(1.D0-U(3)))/epsilon
      FVAL(2)=d2*UXX(2) + FUNC2(U(2)) -
     &     U(2)*(s2*U(1)+k2*U(3))/epsilon
      FVAL(3)=((1.D0-U(3))*U(1)-U(2)*U(3))/epsilon
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
        double precision    d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
        common /ecosolid/   d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
        double precision    DFUNC1,DFUNC2
c
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, j
C-----------------------------------------------------------------------
C
C     ASSIGN DFDU(1:NPDE,1:NPDE), DFDUX(1:NPDE,1:NPDE), AND
C     DFDUXX(1:NPDE,1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      DFDU(1,1) = DFUNC1(U(1)) - (s1*U(2)+k1*(1.D0-U(3)))/epsilon
      DFDU(1,2) = -s1*U(1)/epsilon
      DFDU(1,3) =  k1*U(1)/epsilon

      DFDU(2,1) = -s2*U(2)/epsilon
      DFDU(2,2) = DFUNC2(U(2)) - (s2*U(1)+k2*U(3))/epsilon
      DFDU(2,3) =  -k2*U(2)/epsilon

      DFDU(3,1) = (1.D0-U(3))/epsilon
      DFDU(3,2) = -U(3)/epsilon
      DFDU(3,3) = -(U(1)+U(2))/epsilon
c
      DFDUX(1,1) = 0.D0
      DFDUX(1,2) = 0.D0
      DFDUX(1,3) = 0.D0

      DFDUX(2,1) = 0.D0
      DFDUX(2,2) = 0.D0
      DFDUX(2,3) = 0.D0

      DFDUX(3,1) = 0.D0
      DFDUX(3,2) = 0.D0
      DFDUX(3,3) = 0.D0
c
      DFDUXX(1,1) = d1
      DFDUXX(1,2) = 0.d0
      DFDUXX(1,3) = 0.d0

      DFDUXX(2,1) = 0.d0
      DFDUXX(2,2) = d2
      DFDUXX(2,3) = 0.d0

      DFDUXX(3,1) = 0.d0
      DFDUXX(3,2) = 0.d0
      DFDUXX(3,3) = 0.d0
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
c Initial separataion info
        double precision      U1left, U1right, U2left, U2right,
     &       U3left, U3right, INSEP1(2), INSEP2(2), INSEP3(2)
        common /insep/       U1left, U1right, U2left, U2right,
     &       U3left, U3right, INSEP1, INSEP2, INSEP3
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
C-----------------------------------------------------------------------
c$$$      BVAL(1) = U(1) - U1left
c$$$      BVAL(2) = U(2) - U2left
      BVAL(1) = UX(1)
      BVAL(2) = UX(2)
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
c Initial separataion info
        double precision      U1left, U1right, U2left, U2right,
     &       U3left, U3right, INSEP1(2), INSEP2(2), INSEP3(2)
        common /insep/       U1left, U1right, U2left, U2right,
     &       U3left, U3right, INSEP1, INSEP2, INSEP3
C-----------------------------------------------------------------------
c Loop indices:
        integer                 i
C-----------------------------------------------------------------------
c$$$      BVAL(1) = U(1) - U1right
c$$$      BVAL(2) = U(2) - U2right
      BVAL(1) = UX(1)
      BVAL(2) = UX(2)
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
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, j
C-----------------------------------------------------------------------
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF
C     U(1:NPDE), UX(1:NPDE).
C
      DBDU(1,1) = 0.0D0
      DBDU(1,2) = 0.0D0
      DBDU(1,3) = 0.0D0

      DBDU(2,1) = 0.0D0
      DBDU(2,2) = 0.0D0
      DBDU(2,3) = 0.0D0
c
      DBDUX(1,1) = 1.0D0
      DBDUX(1,2) = 0.0D0
      DBDUX(1,3) = 0.0D0

      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = 1.0D0
      DBDUX(2,3) = 0.0D0
c
      DBDT(1) = 0.0D0
      DBDT(2) = 0.0D0
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
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, j
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF
C     U(1:NPDE), UX(1:NPDE).
C
      DBDU(1,1) = 0.0D0
      DBDU(1,2) = 0.0D0
      DBDU(1,3) = 0.0D0

      DBDU(2,1) = 0.0D0
      DBDU(2,2) = 0.0D0
      DBDU(2,3) = 0.0D0
c
      DBDUX(1,1) = 1.0D0
      DBDUX(1,2) = 0.0D0
      DBDUX(1,3) = 0.0D0

      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = 1.0D0
      DBDUX(2,3) = 0.0D0
c
      DBDT(1) = 0.0D0
      DBDT(2) = 0.0D0
c
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
c-----------------------------------------------------------------------
c Initial separataion info
        double precision      U1left, U1right, U2left, U2right,
     &       U3left, U3right, INSEP1(2), INSEP2(2), INSEP3(2)
        common /insep/       U1left, U1right, U2left, U2right,
     &       U3left, U3right, INSEP1, INSEP2, INSEP3
c-----------------------------------------------------------------------
c Loop indices:
      integer                 i
C-----------------------------------------------------------------------

C
C     ASSIGN U(1:NPDE) THE INITIAL VALUES OF U(T0,X).
C
C     Note these are values for:
C            epsilon = 0.1
C            beta = 1
C           gamma = 0.5
C           lambda = sigma_i = 1
      IF ( X .le. INSEP1(1) ) THEN
         U(1) = U1left
      ELSE IF ( X .le. INSEP1(2)) THEN
         U(1) = (U1right-U1left)*(X-INSEP1(1))/(INSEP1(2) -
     &             INSEP1(1))+U1left
      ELSE
         U(1) = U1right
      ENDIF

      IF ( X .le. INSEP2(1) ) THEN
         U(2) = U2left
      ELSE IF ( X .le. INSEP2(2)) THEN
         U(2) = (U2right-U2left)*(X-INSEP2(1))/(INSEP2(2) -
     &             INSEP2(1))+U2left
      ELSE
         U(2) = U2right
      ENDIF

      IF ( X .le. INSEP3(1) ) THEN
         U(3) = U3left
      ELSE IF ( X .le. INSEP3(2)) THEN
         U(3) = (U3right-U3left)*(X-INSEP3(1))/(INSEP3(2) -
     &             INSEP3(1))+U3left
      ELSE
         U(3) = U3right
      ENDIF

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

C
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNC1(U)
C-----------------------------------------------------------------------
C PURPOSE:
C     THIS FUNCTION is the reaction term for the 1st unknown
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      DOUBLE PRECISION       U

      double precision    d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
      common /ecosolid/   d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
C
C-----------------------------------------------------------------------
      FUNC1 = r1*(U - U**2/C1)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DFUNC1(U)
C-----------------------------------------------------------------------
C PURPOSE: THIS FUNCTION is the derivative of the reaction term for the
C     1st unknown
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      DOUBLE PRECISION       U

      double precision    d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
      common /ecosolid/   d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
C
C-----------------------------------------------------------------------
      DFUNC1 = r1*(1.D0 - 2.D0*U/C1)
      RETURN
      END
      DOUBLE PRECISION FUNCTION FUNC2(U)
C-----------------------------------------------------------------------
C PURPOSE:
C     THIS FUNCTION is the reaction term for the 1st unknown
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      DOUBLE PRECISION       U

      double precision    d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
      common /ecosolid/   d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
C
C-----------------------------------------------------------------------
      FUNC2 = r2*(U - U**2/C2)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DFUNC2(U)
C-----------------------------------------------------------------------
C PURPOSE: THIS FUNCTION is the derivative of the reaction term for the
C     1st unknown
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS:
C INPUT:
      DOUBLE PRECISION       U

      double precision    d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
      common /ecosolid/   d1,d2,s1,s2,k1,k2,epsilon,r1,r2,C1,C2
C
C-----------------------------------------------------------------------
      DFUNC2 = r2*(1.D0 - 2.D0*U/C2)
      RETURN
      END
