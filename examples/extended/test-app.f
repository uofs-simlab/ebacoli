      SUBROUTINE F(T, X, U, UX, UXX, FVAL, NPDE)
      INTEGER          NPDE
      DOUBLE PRECISION T, X, U(NPDE), UX(NPDE), UXX(NPDE)
      DOUBLE PRECISION FVAL(NPDE)
      FVAL(1)= X*UXX(1) + T*X*UX(2) + (1+T*X)*DEXP(-T*X)*U(3)
      FVAL(2)= X*UXX(2) + T*X*UX(1) + (1+T*X)*DEXP(-T*X)*U(4)
      FVAL(3) = -0.5D0*(DEXP(T*X)*U(1) + DSQRT(1-U(4)*U(4)))
      FVAL(4) = -0.5D0*(DEXP(T*X)*U(2) + DSQRT(1-U(3)*U(3)))
      RETURN
      END
      SUBROUTINE DERIVF(T, X, U, UX, UXX, DFDU, DFDUX, DFDUXX, NPDE)
      INTEGER          NPDE
      DOUBLE PRECISION T, X, U(NPDE), UX(NPDE), UXX(NPDE)
      DOUBLE PRECISION DFDU(NPDE,NPDE), DFDUX(NPDE,NPDE),
     &                 DFDUXX(NPDE,NPDE)
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
      DFDUXX(1,1) = X
      DFDUXX(1,2) = 0.D0
      DFDUXX(1,3) = 0.D0
      DFDUXX(1,4) = 0.D0
      DFDUXX(2,1) = 0.D0
      DFDUXX(2,2) = X
      DFDUXX(2,3) = 0.D0
      DFDUXX(2,4) = 0.D0
      DFDUXX(3,1) = 0.D0
      DFDUXX(3,2) = 0.D0
      DFDUXX(3,3) = 0.D0
      DFDUXX(3,4) = 0.D0
      DFDUXX(4,1) = 0.D0
      DFDUXX(4,2) = 0.D0
      DFDUXX(4,3) = 0.D0
      DFDUXX(4,4) = 0.D0
      RETURN
      END

      SUBROUTINE BNDXL(T, U, UX, BVAL, NPDE)
      INTEGER          NPDE
      DOUBLE PRECISION T, U(NPDE), UX(NPDE)
      DOUBLE PRECISION BVAL(NPDE)
      BVAL(1) = U(1) - DEXP(-T)*DSIN(1+T)
      BVAL(2) = U(2) - DEXP(-T)*DSIN(1+T)
      RETURN
      END
      SUBROUTINE BNDXR(T, U, UX, BVAL, NPDE)
      INTEGER          NPDE
      DOUBLE PRECISION T, U(NPDE), UX(NPDE)
      DOUBLE PRECISION BVAL(NPDE)
      BVAL(1) = UX(1) + T*U(1) - DEXP(-2.D0*T)*U(3)
      BVAL(2) = UX(2) + T*U(2) - DEXP(-2.D0*T)*U(4)
      RETURN
      END
      SUBROUTINE DIFBXL(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
      INTEGER          NPDE
      DOUBLE PRECISION T, U(NPDE), UX(NPDE)
      DOUBLE PRECISION DBDU(NPDE,NPDE), DBDUX(NPDE,NPDE),
     &                 DBDT(NPDE)
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
      SUBROUTINE DIFBXR(T, U, UX, DBDU, DBDUX, DBDT, NPDE)
      INTEGER          NPDE
      DOUBLE PRECISION T, U(NPDE), UX(NPDE)
      DOUBLE PRECISION DBDU(NPDE,NPDE), DBDUX(NPDE,NPDE),
     &                 DBDT(NPDE)
      DBDU(1,1) = T
      DBDU(1,2) = 0.D0
      DBDU(1,3) = -DEXP(-2*T)
      DBDU(1,4) = 0.D0
      DBDU(2,1) = 0.D0
      DBDU(2,2) = T
      DBDU(2,3) = 0.D0
      DBDU(2,4) = -DEXP(-2*T)
      DBDUX(1,1) = 1.0D0
      DBDUX(1,2) = 0.0D0
      DBDUX(1,3) = 0.0D0
      DBDUX(1,4) = 0.0D0
      DBDUX(2,1) = 0.0D0
      DBDUX(2,2) = 1.0D0
      DBDUX(2,3) = 0.0D0
      DBDUX(2,4) = 0.0D0
      DBDT(1) = U(1) + 2.D0*DEXP(-2*T)*U(3)
      DBDT(2) = U(2) + 2.D0*DEXP(-2*T)*U(4)
      RETURN
      END

      SUBROUTINE UINIT(X, U, NPDE)
      INTEGER                 NPDE
      DOUBLE PRECISION        X, U(NPDE)
      U(1) = DSIN(X)
      U(2) = DSIN(X)
      U(3) = DCOS(X)
      U(4) = DCOS(X)
      RETURN
      END
      SUBROUTINE TRUU(T, X, U, NPDE)
      INTEGER                 NPDE
      DOUBLE PRECISION        T, X, U(NPDE)
      U(1) = DEXP(-T*X)*DSIN(X+T)
      U(2) = DEXP(-T*X)*DSIN(X+T)
      U(3) = DCOS(X+T)
      U(4) = DCOS(X+T)
      RETURN
      END
