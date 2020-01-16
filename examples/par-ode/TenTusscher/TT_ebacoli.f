C       Problem definition for ten Tusscher Monodomain model
C
C       To make executable:
c >> gfortran -o TT driver_mesh.f f_d_COL1.f f_d_COL2.f f_d_COL3.f f_d_COL4.f f_d_COL5.f  f_d_COL6.f f_d_COL7.f f_d_COL8.f
c f_d_COL9.f f_d_COL10.f f_d_COL11.f f_d_COL12.f f_d_COL13.f  f_d_COL14.f f_d_COL15.f f_d_COL16.f  f_d_COL17.f f_d_COL18.f f_d_COL19.f
c TT.f ebacoli2.f <<
c       To run the executable:
c       ./TT
C-----------------------------------------------------------------------
      SUBROUTINE F(T,X,U,UX,UXX,FVAL,NPDE)
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
c External subroutines (from CellML)
        external computeRates
c and dummy time arg
        double precision voi
c-----------------------------------------------------------------------
C parameters
        double precision CONSTS(53), RATES(19), ALGBRC(70)
        common /TT/ CONSTS, RATES, ALGBRC
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
C-----------------------------------------------------------------------
C
C     ASSIGN FVAL(1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
C Compute rates from CellML generated code
        call computeRates(voi, CONSTS, FVAL, U, ALGBRC)

C and add diffusion term
        FVAL(1) = (1.75D0/1400.0D0) * UXX(1) + FVAL(1)
      RETURN
      END
C-----------------------------------------------------------------------
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
        DOUBLE PRECISION        T,td
C                               THE CURRENT TIME COORDINATE.
C
        DOUBLE PRECISION        X
C                               THE CURRENT SPATIAL COORDINATE.
C
        DOUBLE PRECISION        U(NPDE)
C                               U(1:NPDE) IS THE APPROXIMATION OF THE
C                               SOLUTION AT THE POINT (T,X).
        DOUBLE PRECISION ud(npde)
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
        DOUBLE PRECISION uxxd(npde)
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
      double precision          fval(NPDE)
      DOUBLE PRECISION  col1(NPDE),col2(NPDE),col3(NPDE),col4(NPDE)
      DOUBLE PRECISION  col5(NPDE),col6(NPDE),col7(NPDE),col8(NPDE)
      DOUBLE PRECISION  col9(NPDE),col10(NPDE),col11(NPDE),col12(NPDE)
      DOUBLE PRECISION  col13(NPDE),col14(NPDE),col15(NPDE),col16(NPDE)
      DOUBLE PRECISION  col17(NPDE),col18(NPDE),col19(NPDE)
c-----------------------------------------------------------------------
c External functions
      external f_d_COL1
      external f_d_COL2
      external f_d_COL3
      external f_d_COL4
      external f_d_COL5
      external f_d_COL6
      external f_d_COL7
      external f_d_COL8
      external f_d_COL9
      external f_d_COL10
      external f_d_COL11
      external f_d_COL12
      external f_d_COL13
      external f_d_COL14
      external f_d_COL15
      external f_d_COL16
      external f_d_COL17
      external f_d_COL18
      external f_d_COL19

c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, j
C-----------------------------------------------------------------------
C
C     ASSIGN DFDU(1:NPDE,1:NPDE), DFDUX(1:NPDE,1:NPDE), AND
C     DFDUXX(1:NPDE,1:NPDE) ACCORDING TO THE RIGHT HAND SIDE OF THE PDE
C     IN TERMS OF U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c this line returns D fval(i) / D U(1)  whichis the first column of the Jacobian matrix
      call f_d_COL1(t,x,u,ud,ux,uxx,uxxd,fval,col1,npde)
c this line returns D fval(i) / D U(2)  whichis the 2nd column of the Jacobian matrix
      call f_d_COL2(t,x,u,ud,ux,uxx,uxxd,fval,col2,npde)
c this line returns D fval(i) / D U(3)  whichis the 3rd column of the Jacobian matrix
      call f_d_COL3(t,x,u,ud,ux,uxx,uxxd,fval,col3,npde)
c this line returns D fval(i) / D U(4)  whichis the 4th column of the Jacobian matrix
      call f_d_COL4(t,x,u,ud,ux,uxx,uxxd,fval,col4,npde)
c this line returns D fval(i) / D U(5)  whichis the 5th column of the Jacobian matrix
      call f_d_COL5(t,x,u,ud,ux,uxx,uxxd,fval,col5,npde)
c this line returns D fval(i) / D U(6)  whichis the 6th column of the Jacobian matrix
      call f_d_COL6(t,x,u,ud,ux,uxx,uxxd,fval,col6,npde)
c this line returns D fval(i) / D U(7)  whichis the 7th column of the Jacobian matrix
      call f_d_COL7(t,x,u,ud,ux,uxx,uxxd,fval,col7,npde)
c this line returns D fval(i) / D U(8)  whichis the 8th column of the Jacobian matrix
      call f_d_COL8(t,x,u,ud,ux,uxx,uxxd,fval,col8,npde)
c this line returns D fval(i) / D U(9)  whichis the 9th column of the Jacobian matrix
      call f_d_COL9(t,x,u,ud,ux,uxx,uxxd,fval,col9,npde)
c this line returns D fval(i) / D U(10)  whichis the 10th column of the Jacobian matrix
      call f_d_COL10(t,x,u,ud,ux,uxx,uxxd,fval,col10,npde)
c this line returns D fval(i) / D U(11)  whichis the 11th column of the Jacobian matrix
      call f_d_COL11(t,x,u,ud,ux,uxx,uxxd,fval,col11,npde)
c this line returns D fval(i) / D U(12)  whichis the 12th column of the Jacobian matrix
      call f_d_COL12(t,x,u,ud,ux,uxx,uxxd,fval,col12,npde)
c this line returns D fval(i) / D U(13)  whichis the 13th column of the Jacobian matrix
      call f_d_COL13(t,x,u,ud,ux,uxx,uxxd,fval,col13,npde)
c this line returns D fval(i) / D U(14)  whichis the 14th column of the Jacobian matrix
      call f_d_COL14(t,x,u,ud,ux,uxx,uxxd,fval,col14,npde)
c this line returns D fval(i) / D U(15)  whichis the 15th column of the Jacobian matrix
      call f_d_COL15(t,x,u,ud,ux,uxx,uxxd,fval,col15,npde)
c this line returns D fval(i) / D U(16)  whichis the 16th column of the Jacobian matrix
      call f_d_COL16(t,x,u,ud,ux,uxx,uxxd,fval,col16,npde)
c this line returns D fval(i) / D U(17)  whichis the 17th column of the Jacobian matrix
      call f_d_COL17(t,x,u,ud,ux,uxx,uxxd,fval,col17,npde)
c this line returns D fval(i) / D U(18)  whichis the 18th column of the Jacobian matrix
      call f_d_COL18(t,x,u,ud,ux,uxx,uxxd,fval,col18,npde)
c this line returns D fval(i) / D U(19)  whichis the 19th column of the Jacobian matrix
      call f_d_COL19(t,x,u,ud,ux,uxx,uxxd,fval,col19,npde)
c fill the 1st col of the jacobian
      do i=1,npde
          DFDU(i,1) = col1(i)
      enddo
c fill the 2nd col of the jacobian
      do i=1,npde
          DFDU(i,2) = col2(i)
      enddo
c fill the 3rd col of the jacobian
      do i=1,npde
          DFDU(i,3) = col3(i)
      enddo
c fill the 4th col of the jacobian
      do i=1,npde
          DFDU(i,4) = col4(i)
      enddo
c fill the 5th col of the jacobian
      do i=1,npde
          DFDU(i,5) = col5(i)
      enddo
c fill the 6th col of the jacobian
      do i=1,npde
          DFDU(i,6) = col6(i)
      enddo
c fill the 7th col of the jacobian
      do i=1,npde
          DFDU(i,7) = col7(i)
      enddo
c fill the 8th col of the jacobian
      do i=1,npde
          DFDU(i,8) = col8(i)
      enddo
c fill the 9th col of the jacobian
      do i=1,npde
          DFDU(i,9) = col9(i)
      enddo
c fill the 10th col of the jacobian
      do i=1,npde
          DFDU(i,10) = col10(i)
      enddo
c fill the 11th col of the jacobian
      do i=1,npde
          DFDU(i,11) = col11(i)
      enddo
c fill the 12th col of the jacobian
      do i=1,npde
          DFDU(i,12) = col12(i)
      enddo
c fill the 13th col of the jacobian
      do i=1,npde
          DFDU(i,13) = col13(i)
      enddo
c fill the 14th col of the jacobian
      do i=1,npde
          DFDU(i,14) = col14(i)
      enddo
c fill the 15th col of the jacobian
      do i=1,npde
          DFDU(i,15) = col15(i)
      enddo
c fill the 16th col of the jacobian
      do i=1,npde
          DFDU(i,16) = col16(i)
      enddo
c fill the 17th col of the jacobian
      do i=1,npde
          DFDU(i,17) = col17(i)
      enddo
c fill the 18th col of the jacobian
      do i=1,npde
          DFDU(i,18) = col18(i)
      enddo
c fill the 19th col of the jacobian
      do i=1,npde
          DFDU(i,19) = col19(i)
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,npde
        do j=1,npde
          DFDUX(i,j) = 0.0D0
        enddo
      enddo
c first row of the matrix DFDUX is also zero, because in the First equation (Monodomain) we only have a U_xx term and I_ion
c where we know that I_ion is a rate variable for one of the ODEs and for sure does not have any spatial derivative.
C
      do i=1,npde
        do j=1,npde
          DFDUXX(i,j) = 0.0D0
        enddo
      enddo
c only one element is not zero, and that is the element corresponding to the first equation (PDE: Monodomain)
      DFDUXX(1,1) = 1.75D0/1400.D0
      RETURN
      END
C-----------------------------------------------------------------------
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
      BVAL(1) = UX(1)
C
      RETURN
      END
C-----------------------------------------------------------------------
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
      DOUBLE PRECISION coeff1
      COMMON /BURGER/ coeff1
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i
C-----------------------------------------------------------------------
      BVAL(1) = UX(1)
C
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
C-----------------------------------------------------------------------
C External functions
        external initConsts
C-----------------------------------------------------------------------
C parameters
        double precision CONSTS(53), RATES(19), ALGBRC(70)
        common /TT/ CONSTS, RATES, ALGBRC
C-----------------------------------------------------------------------
        double precision      Ulow, Uhigh, U2low
C-----------------------------------------------------------------------
c Loop indices:
        integer                 i
      INTRINSIC SIN
C-----------------------------------------------------------------------
C
C     ASSIGN U(1:NPDE) THE INITIAL VALUES OF U(T0,X).
C

C     set default states/constants
      call initConsts(CONSTS,U)

C     Low state set to equilibrium value
      Ulow  = U(1)
      Uhigh = -50

      IF ( x .le. 10.D0 ) THEN
         U(1) = (Ulow-Uhigh)*x/10.D0 + Uhigh
      ENDIF
C
      RETURN
      END
C-----------------------------------------------------------------------
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
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF
C     U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      DBDU(1,1) = 0.0D0
      DBDU(1,2) = 0.0D0
      DBDU(1,3) = 0.0D0
      DBDU(1,4) = 0.0D0
      DBDU(1,5) = 0.0D0
      DBDU(1,6) = 0.0D0
      DBDU(1,7) = 0.0D0
      DBDU(1,8) = 0.0D0
      DBDU(1,9) = 0.0D0
      DBDU(1,10) = 0.0D0
      DBDU(1,11) = 0.0D0
      DBDU(1,12) = 0.0D0
      DBDU(1,13) = 0.0D0
      DBDU(1,14) = 0.0D0
      DBDU(1,15) = 0.0D0
      DBDU(1,16) = 0.0D0
      DBDU(1,17) = 0.0D0
      DBDU(1,18) = 0.0D0
      DBDU(1,19) = 0.0D0

      DBDUX(1,1) = 1.0D0
      DBDUX(1,2) = 0.0D0
      DBDUX(1,3) = 0.0D0
      DBDUX(1,4) = 0.0D0
      DBDUX(1,5) = 0.0D0
      DBDUX(1,6) = 0.0D0
      DBDUX(1,7) = 0.0D0
      DBDUX(1,8) = 0.0D0
      DBDUX(1,9) = 0.0D0
      DBDUX(1,10) = 0.0D0
      DBDUX(1,11) = 0.0D0
      DBDUX(1,12) = 0.0D0
      DBDUX(1,13) = 0.0D0
      DBDUX(1,14) = 0.0D0
      DBDUX(1,15) = 0.0D0
      DBDUX(1,16) = 0.0D0
      DBDUX(1,17) = 0.0D0
      DBDUX(1,18) = 0.0D0
      DBDUX(1,19) = 0.0D0
C
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
c-----------------------------------------------------------------------
c Loop indices:
        integer                 i, j
C-----------------------------------------------------------------------
C
C     ASSIGN DBDU(1:NPDE,1:NPDE), DBDU(1:NPDE,1:NPDE), AND DBDT(1:NPDE)
C     ACCORDING TO THE RIGHT BOUNDARY CONDITION EQUATION IN TERMS OF
C     U(1:NPDE), UX(1:NPDE), UXX(1:NPDE).
C
      DBDU(1,1) = 0.0D0
      DBDU(1,2) = 0.0D0
      DBDU(1,3) = 0.0D0
      DBDU(1,4) = 0.0D0
      DBDU(1,5) = 0.0D0
      DBDU(1,6) = 0.0D0
      DBDU(1,7) = 0.0D0
      DBDU(1,8) = 0.0D0
      DBDU(1,9) = 0.0D0
      DBDU(1,10) = 0.0D0
      DBDU(1,11) = 0.0D0
      DBDU(1,12) = 0.0D0
      DBDU(1,13) = 0.0D0
      DBDU(1,14) = 0.0D0
      DBDU(1,15) = 0.0D0
      DBDU(1,16) = 0.0D0
      DBDU(1,17) = 0.0D0
      DBDU(1,18) = 0.0D0
      DBDU(1,19) = 0.0D0

      DBDUX(1,1) = 1.0D0
      DBDUX(1,2) = 0.0D0
      DBDUX(1,3) = 0.0D0
      DBDUX(1,4) = 0.0D0
      DBDUX(1,5) = 0.0D0
      DBDUX(1,6) = 0.0D0
      DBDUX(1,7) = 0.0D0
      DBDUX(1,8) = 0.0D0
      DBDUX(1,9) = 0.0D0
      DBDUX(1,10) = 0.0D0
      DBDUX(1,11) = 0.0D0
      DBDUX(1,12) = 0.0D0
      DBDUX(1,13) = 0.0D0
      DBDUX(1,14) = 0.0D0
      DBDUX(1,15) = 0.0D0
      DBDUX(1,16) = 0.0D0
      DBDUX(1,17) = 0.0D0
      DBDUX(1,18) = 0.0D0
      DBDUX(1,19) = 0.0D0

      DBDT(1) = 0.0D0

      RETURN
      END
      real FUNCTION TERNRY(TEST, VALA, VALB)
      LOGICAL TEST
      double precision VALA, VALB
      IF (TEST) THEN
        TERNRY = VALA
      ELSE
        TERNRY = VALB
      ENDIF
      RETURN
      END
