c     This is the problem defintion for the monodomain equation with
c     the Fitzhugh--Nagumo cell model.
c
c     ** We should insert parameters into this. **

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
        double precision    epsilon, beta, gamma, sigmai, sigmae
        common /fhn/       epsilon, beta, gamma, sigmai, sigmae
        double precision   istim
c
c-----------------------------------------------------------------------
c Loop indices:
      integer                 i
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      if (X .le. 35 .and. X .ge. 30) then
         if (T .le. 5) then
            istim = -0.5d0
         else if (T .le. 10) then
            istim = -1.d0
         else if (T .le. 11) then
            istim = -2.0d0
         else
            istim = 0.d0
         end if
      else
         istim = 0
      end if

c     The parabolic equation
      FVAL(1)= sigmai*( UXX(1) + UXX(5) ) +
     &     ( U(1) - U(1)*U(1)*U(1)/3.D0 - U(3) ) / epsilon - istim

      FVAL(2)= sigmai*( UXX(2) + UXX(6) ) +
     &     ( U(2) - U(2)*U(2)*U(2)/3.D0 - U(4) ) / epsilon - istim

C     The cell model equation
      FVAL(3) = epsilon * ( U(1) + beta - gamma*U(3) )
C
      FVAL(4) = epsilon * ( U(2) + beta - gamma*U(4) )
C
      FVAL(5) = sigmai*UXX(1) +(sigmai+sigmae)*UXX(5)
C
      FVAL(6) = sigmai*UXX(2) +(sigmai+sigmae)*UXX(6)

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
        double precision    epsilon, beta, gamma, sigmai, sigmae
        common /fhn/       epsilon, beta, gamma, sigmai, sigmae
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
C     TODO index arrays in column major
      DFDU   = 0.0D0
      DFDUX  = 0.0D0
      DFDUXX = 0.0D0

      DFDU(1,1) = ( 1.D0 - U(1)*U(1) ) / epsilon
      DFDU(1,3) = -1.D0/epsilon

      DFDU(2,2) = ( 1.D0 - U(2)*U(2) ) / epsilon
      DFDU(2,4) = -1.D0/epsilon

      DFDU(3,1) = epsilon
      DFDU(3,3) = -gamma*epsilon

      DFDU(4,2) = epsilon
      DFDU(4,4) = -gamma*epsilon

      DFDUXX(1,1) = sigmai
      DFDUXX(1,5) = sigmai

      DFDUXX(2,2) = sigmai
      DFDUXX(2,6) = sigmai

      DFDUXX(5,1) = sigmai
      DFDUXX(5,5) = sigmai+sigmae

      DFDUXX(6,2) = sigmai
      DFDUXX(6,6) = sigmai+sigmae
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
c Loop indices:
        integer                 i
C-----------------------------------------------------------------------
      BVAL(1) = UX(1) - UX(5)
      BVAL(5) = UX(5)
C
      BVAL(2) = UX(2) - UX(6)
      BVAL(6) = UX(6)
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
C-----------------------------------------------------------------------
c Loop indices:
        integer                 i
C-----------------------------------------------------------------------
      BVAL(1) = UX(1) - UX(5)
      BVAL(5) = U(5)
C
      BVAL(2) = UX(2) - UX(6)
      BVAL(6) = U(6)

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
      DBDU  = 0.0D0
      DBDUX = 0.0D0
      DBDT  = 0.0D0
c
      DBDUX(1,1) = 1.0D0
      DBDUX(2,2) = 1.0D0
      DBDUX(1,5) = -1.0D0
      DBDUX(2,6) = -1.0D0
c
      DBDUX(5,5) = 1.0D0
      DBDUX(6,6) = 1.0D0
c
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
      DBDU  = 0.0D0
      DBDUX = 0.0D0
      DBDT  = 0.0D0
c
      DBDUX(1,1) = 1.0D0
      DBDUX(2,2) = 1.0D0
      DBDUX(1,5) = -1.0D0
      DBDUX(2,6) = -1.0D0
c
      DBDU(5,5) = 1.0D0
      DBDU(6,6) = 1.0D0

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
        double precision      Ulow, Uhigh, U2low
c-----------------------------------------------------------------------
c Loop indices:
      integer                 i
C-----------------------------------------------------------------------
      Ulow  = -1.2879118919372559D0
      Uhigh = 2.D0
      U2low = -0.5758181214332581D0

C
C     ASSIGN U(1:NPDE) THE INITIAL VALUES OF U(T0,X).
C
C     Note these are values for:
C            epsilon = 0.1
C            beta = 1
C           gamma = 0.5
C           lambda = sigmai = 1
c$$$      IF ( x .le. 10.D0 ) THEN
c$$$         U(1) = (Ulow-Uhigh)*x/10.D0 + Uhigh
c$$$      ELSE
c$$$         U(1) = Ulow
c$$$      ENDIF
      U(1) = Ulow
      U(2) = Ulow
      U(3) = U2low
      U(4) = U2low
      U(5) = 0.0D0
      U(6) = 0.0D0
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
