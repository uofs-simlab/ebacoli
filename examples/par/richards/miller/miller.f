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
      BVAL(1) = U(1) - 0.1d0
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
      BVAL(1) = U(1)
c$$$      BVAL(1) = UX(1)
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
c$$$      DBDU(1,1) = 0.0D0
c
      DBDUX(1,1) = 0.0D0
c$$$      DBDUX(1,1) = 1.0D0
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
      double precision         xsplit
C-----------------------------------------------------------------------
      xsplit= 0.1
      if (X > xsplit) then
         U(1) = X - 10.d0
      else
         U(1) = 0.1d0 + (xsplit-10.d0-0.1d0)*X/xsplit
      endif
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
      if (psi .lt. 0) then
         hydcon = ks * (0.1D1 - (alpha * dabs(psi)) ** (n - 0.1D1) *
     $        (0.1D1 + (alpha * dabs(psi)) ** n) ** (-m)) ** 2 / (0.1D1
     $        + (alpha* dabs (psi)) ** n) ** (m / 0.2D1)
      else
         hydcon = thetas
      endif
      return
      end
      doubleprecision function dhydco (psi)
      doubleprecision psi
C     Derivative of hydraulic conductivity K'(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      if (psi .lt. 0) then
         dhydco = 0.2D1 * ks * (0.1D1 - (alpha * dabs(psi)) ** (n -
     $        0.1D1)* (0.1D1 + (alpha * dabs(psi)) ** n) ** (-m)) /
     $        (0.1D1 +(alpha * dabs(psi)) ** n) ** (m / 0.2D1) * (
     $        -(alpha *dabs(psi)) ** (n - 0.1D 1) * (n - 0.1D1) *
     $        dabs(psi) / psi /dabs(psi) * (0.1D1 + (alpha * a bs(psi))
     $        ** n) ** (-m) +(alpha * dabs(psi)) ** (n - 0.1D1) * (0.1D
     $        1 + (alpha *dabs(psi)) ** n) ** (-m) * m * (alpha *
     $        dabs(psi)) ** n * n *dabs(psi) / psi / dabs(psi) / (0.1D1
     $        + (alpha * dabs(psi)) **n)) - ks * (0.1D1 - (alpha *
     $        dabs(psi)) ** (n - 0.1D1) *(0.1D1 + ( alpha * dabs(psi))
     $        ** n) ** (-m)) ** 2 / (0.1D1 +(alpha * dabs(psi)) ** n) **
     $        (m / 0.2D1) * m * (alpha *dabs(psi)) ** n * n * dabs(psi)
     $        / psi / dabs(psi) / (0.1D1 +(alpha * dabs(psi)) ** n) /
     $        0.2D1
      else
         dhydco = 0.d0
      endif
      return
      end
      doubleprecision function ddhydc (psi)
        doubleprecision psi
        doubleprecision s0
        doubleprecision s1
        doubleprecision s10
        doubleprecision s100
        doubleprecision s101
        doubleprecision s102
        doubleprecision s103
        doubleprecision s104
        doubleprecision s105
        doubleprecision s106
        doubleprecision s107
        doubleprecision s108
        doubleprecision s109
        doubleprecision s11
        doubleprecision s110
        doubleprecision s111
        doubleprecision s112
        doubleprecision s113
        doubleprecision s114
        doubleprecision s115
        doubleprecision s116
        doubleprecision s117
        doubleprecision s118
        doubleprecision s119
        doubleprecision s12
        doubleprecision s120
        doubleprecision s121
        doubleprecision s122
        doubleprecision s123
        doubleprecision s124
        doubleprecision s125
        doubleprecision s126
        doubleprecision s127
        doubleprecision s128
        doubleprecision s129
        doubleprecision s13
        doubleprecision s130
        doubleprecision s131
        doubleprecision s132
        doubleprecision s133
        doubleprecision s134
        doubleprecision s135
        doubleprecision s136
        doubleprecision s137
        doubleprecision s138
        doubleprecision s139
        doubleprecision s14
        doubleprecision s140
        doubleprecision s141
        doubleprecision s142
        doubleprecision s143
        doubleprecision s144
        doubleprecision s145
        doubleprecision s146
        doubleprecision s147
        doubleprecision s148
        doubleprecision s149
        doubleprecision s15
        doubleprecision s150
        doubleprecision s151
        doubleprecision s152
        doubleprecision s153
        doubleprecision s154
        doubleprecision s155
        doubleprecision s156
        doubleprecision s157
        doubleprecision s158
        doubleprecision s159
        doubleprecision s16
        doubleprecision s160
        doubleprecision s161
        doubleprecision s162
        doubleprecision s163
        doubleprecision s164
        doubleprecision s165
        doubleprecision s166
        doubleprecision s167
        doubleprecision s168
        doubleprecision s169
        doubleprecision s17
        doubleprecision s170
        doubleprecision s171
        doubleprecision s172
        doubleprecision s173
        doubleprecision s174
        doubleprecision s175
        doubleprecision s176
        doubleprecision s177
        doubleprecision s178
        doubleprecision s179
        doubleprecision s18
        doubleprecision s180
        doubleprecision s181
        doubleprecision s182
        doubleprecision s183
        doubleprecision s184
        doubleprecision s185
        doubleprecision s186
        doubleprecision s187
        doubleprecision s188
        doubleprecision s189
        doubleprecision s19
        doubleprecision s190
        doubleprecision s191
        doubleprecision s192
        doubleprecision s193
        doubleprecision s194
        doubleprecision s195
        doubleprecision s196
        doubleprecision s197
        doubleprecision s198
        doubleprecision s199
        doubleprecision s2
        doubleprecision s20
        doubleprecision s200
        doubleprecision s201
        doubleprecision s202
        doubleprecision s203
        doubleprecision s204
        doubleprecision s205
        doubleprecision s206
        doubleprecision s207
        doubleprecision s208
        doubleprecision s209
        doubleprecision s21
        doubleprecision s210
        doubleprecision s211
        doubleprecision s212
        doubleprecision s213
        doubleprecision s214
        doubleprecision s215
        doubleprecision s216
        doubleprecision s217
        doubleprecision s218
        doubleprecision s219
        doubleprecision s22
        doubleprecision s220
        doubleprecision s221
        doubleprecision s222
        doubleprecision s223
        doubleprecision s224
        doubleprecision s225
        doubleprecision s226
        doubleprecision s227
        doubleprecision s228
        doubleprecision s229
        doubleprecision s23
        doubleprecision s230
        doubleprecision s231
        doubleprecision s232
        doubleprecision s233
        doubleprecision s234
        doubleprecision s235
        doubleprecision s236
        doubleprecision s237
        doubleprecision s24
        doubleprecision s25
        doubleprecision s26
        doubleprecision s27
        doubleprecision s28
        doubleprecision s29
        doubleprecision s3
        doubleprecision s30
        doubleprecision s31
        doubleprecision s32
        doubleprecision s33
        doubleprecision s34
        doubleprecision s35
        doubleprecision s36
        doubleprecision s37
        doubleprecision s38
        doubleprecision s39
        doubleprecision s4
        doubleprecision s40
        doubleprecision s41
        doubleprecision s42
        doubleprecision s43
        doubleprecision s44
        doubleprecision s45
        doubleprecision s46
        doubleprecision s47
        doubleprecision s48
        doubleprecision s49
        doubleprecision s5
        doubleprecision s50
        doubleprecision s51
        doubleprecision s52
        doubleprecision s53
        doubleprecision s54
        doubleprecision s55
        doubleprecision s56
        doubleprecision s57
        doubleprecision s58
        doubleprecision s59
        doubleprecision s6
        doubleprecision s60
        doubleprecision s61
        doubleprecision s62
        doubleprecision s63
        doubleprecision s64
        doubleprecision s65
        doubleprecision s66
        doubleprecision s67
        doubleprecision s68
        doubleprecision s69
        doubleprecision s7
        doubleprecision s70
        doubleprecision s71
        doubleprecision s72
        doubleprecision s73
        doubleprecision s74
        doubleprecision s75
        doubleprecision s76
        doubleprecision s77
        doubleprecision s78
        doubleprecision s79
        doubleprecision s8
        doubleprecision s80
        doubleprecision s81
        doubleprecision s82
        doubleprecision s83
        doubleprecision s84
        doubleprecision s85
        doubleprecision s86
        doubleprecision s87
        doubleprecision s88
        doubleprecision s89
        doubleprecision s9
        doubleprecision s90
        doubleprecision s91
        doubleprecision s92
        doubleprecision s93
        doubleprecision s94
        doubleprecision s95
        doubleprecision s96
        doubleprecision s97
        doubleprecision s98
        doubleprecision s99
C     2nd derivative of hydraulic conductivity K''(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      if (psi .lt. 0) then
         s0 = alpha * dabs(psi)
         s1 = n - 0.1D1
         s2 = dabs(psi)
         s3 = alpha * s2
         s4 = 0.1D1 + s3 ** n
         s5 = -m
         s6 = s0 ** s1 * s4 ** s5
         s7 = 0.1D1 / s2 * dabs(psi) / psi
         s8 = alpha * s2
         s9 = n - 0.1D1
         s10 = alpha * s2
         s11 = 0.1D1 + s10 ** n
         s12 = -m
         s13 = alpha * s2
         s14 = alpha * s2
         s15 = 0.1D1 + s14 ** n
         s16 = s7 * s11 ** s12
         s17 = s13 ** n * s8 ** s9
         s18 = s16 * s17
         s19 = (-n + 0.1D1) * s6 * s7 + m * n * s18 / s15
         s20 = alpha * s2
         s21 = 0.1D1 + s20 ** n
         s22 = m / 0.2D1
         s23 = s21 ** s22
         s24 = alpha * s2
         s25 = n - 0.1D1
         s26 = alpha * s2
         s27 = 0.1D1 + s26 ** n
         s28 = -m
         s29 = alpha * s2
         s30 = 0.1D1 + s29 ** n
         s31 = m / 0.2D1
         s32 = s30 ** s31
         s33 = alpha * s2
         s34 = n - 0.1D1
         s35 = alpha * s2
         s36 = 0.1D1 + s35 ** n
         s37 = -m
         s38 = s7 * s33 ** s34
         s39 = alpha * s2
         s40 = n - 0.1D1
         s41 = alpha * s2
         s42 = 0.1D1 + s41 ** n
         s43 = -m
         s44 = alpha * s2
         s45 = alpha * s2
         s46 = 0.1D1 + s45 ** n
         s47 = s7 * s39 ** s40
         s48 = s42 ** s43 * s44 ** n
         s49 = s47 * s48
         s50 = alpha * s2
         s51 = alpha * s2
         s52 = 0.1D1 + s51 ** n
         s53 = s7 * s50 ** n
         s54 = 0.1D1 / s32 / s52
         s55 = (0.1D1 - s24 ** s25 * s27 ** s28) * ((-n + 0.1D1) * s38 *
     $        s36 ** s37 + m * n * s49 / s46)
         s56 = s53 * s54
         s57 = alpha * s2
         s58 = n - 0.1D1
         s59 = alpha * s2
         s60 = 0.1D1 + s59 ** n
         s61 = -m
         s62 = alpha * s2
         s63 = 0.1D1 + s62 ** n
         s64 = m / 0.2D1
         s65 = s63 ** s64
         s66 = alpha * s2
         s67 = n - 0.1D1
         s68 = dabs(psi) / psi
         s69 = alpha * s2
         s70 = 0.1D1 + s69 ** n
         s71 = -m
         s72 = s66 ** s67 * s70 ** s71
         s73 = 0.1D1 / s2 ** 2 * s68 ** 2
         s74 = alpha * s2
         s75 = n - 0.1D1
         s76 = alpha * s2
         s77 = 0.1D1 + s76 ** n
         s78 = -m
         s79 = s74 ** s75 * s77 ** s78
         s80 = 0.1D1 / s2 * 0.d0
         s81 = alpha * s2
         s82 = n - 0.1D1
         s83 = alpha * s2
         s84 = 0.1D1 + s83 ** n
         s85 = -m
         s86 = s73 * s81 ** s82
         s87 = alpha * s2
         s88 = n - 0.1D1
         s89 = alpha * s2
         s90 = 0.1D1 + s89 ** n
         s91 = -m
         s92 = alpha * s2
         s93 = alpha * s2
         s94 = 0.1D1 + s93 ** n
         s95 = s73 * s87 ** s88
         s96 = s90 ** s91 * s92 ** n
         s97 = s95 * s96
         s98 = alpha * s2
         s99 = n - 0.1D1
         s100 = alpha * s2
         s101 = 0.1D1 + s100 ** n
         s102 = -m
         s103 = alpha * s2
         s104 = s103 ** n
         s105 = alpha * s2
         s106 = 0.1D1 + s105 ** n
         s107 = s73 * s101 ** s102
         s108 = s98 ** s99 * s104 ** 2
         s109 = s107 * s108
         s110 = alpha * s2
         s111 = n - 0.1D1
         s112 = alpha * s2
         s113 = 0.1D1 + s112 ** n
         s114 = -m
         s115 = alpha * s2
         s116 = alpha * s2
         s117 = 0.1D1 + s116 ** n
         s118 = s73 * s110 ** s111
         s119 = s113 ** s114 * s115 ** n
         s120 = s118 * s119
         s121 = alpha * s2
         s122 = n - 0.1D1
         s123 = alpha * s2
         s124 = 0.1D1 + s123 ** n
         s125 = -m
         s126 = alpha * s2
         s127 = alpha * s2
         s128 = 0.1D1 + s127 ** n
         s129 = s80 * s121 ** s122
         s130 = s124 ** s125 * s126 ** n
         s131 = s129 * s130
         s132 = alpha * s2
         s133 = n - 0.1D1
         s134 = alpha * s2
         s135 = 0.1D1 + s134 ** n
         s136 = -m
         s137 = alpha * s2
         s138 = alpha * s2
         s139 = 0.1D1 + s138 ** n
         s140 = s73 * s132 ** s133
         s141 = s135 ** s136 * s137 ** n
         s142 = s140 * s141
         s143 = alpha * s2
         s144 = n - 0.1D1
         s145 = alpha * s2
         s146 = 0.1D1 + s145 ** n
         s147 = -m
         s148 = alpha * s2
         s149 = s148 ** n
         s150 = alpha * s2
         s151 = 0.1D1 + s150 ** n
         s152 = s73 * s143 ** s144
         s153 = s146 ** s147 * s149 ** 2
         s154 = s152 * s153
         s155 = 0.1D1 / s65 * (0.1D1 - s57 ** s58 * s60 ** s61)
         s156 = alpha * s2
         s157 = n - 0.1D1
         s158 = alpha * s2
         s159 = 0.1D1 + s158 ** n
         s160 = -m
         s161 = 0.1D1 - s156 ** s157 * s159 ** s160
         s162 = alpha * s2
         s163 = 0.1D1 + s162 ** n
         s164 = m / 0.2D1
         s165 = s163 ** s164
         s166 = alpha * s2
         s167 = s166 ** n
         s168 = alpha * s2
         s169 = 0.1D1 + s168 ** n
         s170 = s73 * s161 ** 2
         s171 = 0.1D1 / s165 * s167 ** 2
         s172 = s170 * s171
         s173 = alpha * s2
         s174 = n - 0.1D1
         s175 = alpha * s2
         s176 = 0.1D1 + s175 ** n
         s177 = -m
         s178 = 0.1D1 - s173 ** s174 * s176 ** s177
         s179 = alpha * s2
         s180 = 0.1D1 + s179 ** n
         s181 = m / 0.2D1
         s182 = s180 ** s181
         s183 = alpha * s2
         s184 = alpha * s2
         s185 = 0.1D1 + s184 ** n
         s186 = s73 * s183 ** n
         s187 = s178 ** 2 / s182
         s188 = s186 * s187
         s189 = alpha * s2
         s190 = n - 0.1D1
         s191 = alpha * s2
         s192 = 0.1D1 + s191 ** n
         s193 = -m
         s194 = 0.1D1 - s189 ** s190 * s192 ** s193
         s195 = alpha * s2
         s196 = 0.1D1 + s195 ** n
         s197 = m / 0.2D1
         s198 = s196 ** s197
         s199 = alpha * s2
         s200 = alpha * s2
         s201 = 0.1D1 + s200 ** n
         s202 = s80 * s199 ** n
         s203 = s194 ** 2 / s198
         s204 = s202 * s203
         s205 = alpha * s2
         s206 = n - 0.1D1
         s207 = alpha * s2
         s208 = 0.1D1 + s207 ** n
         s209 = -m
         s210 = 0.1D1 - s205 ** s206 * s208 ** s209
         s211 = alpha * s2
         s212 = 0.1D1 + s211 ** n
         s213 = m / 0.2D1
         s214 = s212 ** s213
         s215 = alpha * s2
         s216 = alpha * s2
         s217 = 0.1D1 + s216 ** n
         s218 = s73 * s215 ** n
         s219 = s210 ** 2 / s214
         s220 = s218 * s219
         s221 = alpha * s2
         s222 = n - 0.1D1
         s223 = alpha * s2
         s224 = 0.1D1 + s223 ** n
         s225 = -m
         s226 = 0.1D1 - s221 ** s222 * s224 ** s225
         s227 = alpha * s2
         s228 = 0.1D1 + s227 ** n
         s229 = m / 0.2D1
         s230 = s228 ** s229
         s231 = alpha * s2
         s232 = s231 ** n
         s233 = alpha * s2
         s234 = 0.1D1 + s233 ** n
         s235 = s73 * s226 ** 2
         s236 = 0.1D1 / s230 * s232 ** 2
         s237 = s235 * s236
         ddhydc = 0.2D1 * ks * s19 ** 2 / s23 - 0.2D1 * ks * m * n * s55
     $        * s5 6 + 0.2D1 * ks * s155 * (-(n - 0.1D1) ** 2 * s72 *
     $        s73 + (-n + 0.1 D1) * s79 * s80 + (n - 0.1D1) * s86 * s84
     $        ** s85 + (0.2D1 * n - 0. 2D1) * m * n * s97 / s94 - m ** 2
     $        * n ** 2 * s109 / s106 ** 2 + m * n ** 2 * s120 / s117 + m
     $        * n * s131 / s128 - m * n * s142 / s139 - m * n ** 2 *
     $        s154 / s151 ** 2) + ks * m ** 2 * n ** 2 * s172 / s169 **
     $        2 / 0.4D1 - ks * m * n ** 2 * s188 / s185 / 0.2D1 - ks * m
     $        * n * s204 / s201 / 0.2D1 + ks * m * n * s220 / s217 /
     $        0.2D1 + ks * m * n ** 2 * s237 / s234 ** 2 / 0.2D1
      else
         dhydco = 0.d0
      endif
      return
      end
      doubleprecision function volcon (psi)
      doubleprecision psi
C     Volumetric water content theta(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      if (psi .lt. 0) then
         volcon = (thetas - thetar) / (0.1D1 + (alpha * dabs(psi)) ** n)
     $        ** m + thetar
      else
         volcon = ks
      endif
      return
      end
      doubleprecision function dvolco (psi)
      doubleprecision psi
C     Derivative of volumetric water content theta'(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      if (psi .lt. 0) then
         dvolco = -(thetas - thetar) / (0.1D1 + (alpha * dabs(psi)) **
     $        n)** m * m * (alpha * dabs(psi)) ** n * n * dabs(psi) /
     $        psi /dabs(psi) / (0.1D1 + (alpha * dabs(psi)) ** n)
      else
         dvolco = 0.d0
      endif
      return
      end
      doubleprecision function ddvolc (psi)
      doubleprecision psi
C     2nd derivative of Volumetric water content theta''(psi)
      double precision ss,phi,thetas,thetar,alpha,n,m,ks
      common /richards/ ss,phi,thetas,thetar,alpha,n,m,ks
      if (psi .lt. 0) then
         ddvolc = (thetas - thetar) / (0.1D1 + (alpha * dabs(psi)) ** n)
     $        ** m * m ** 2 * ((alpha * dabs(psi)) ** n) ** 2 * n ** 2
     $        *dabs(psi) / psi ** 2 / dabs(psi) ** 2 / (0.1D1 + (alpha
     $        *dabs(psi)) ** n) ** 2 - (thetas - thetar) / (0.1D1 +
     $        (alpha *dabs(psi)) ** n) ** m * m * (alpha * dabs(psi)) **
     $        n * n ** 2* dabs(psi) / psi ** 2 / dabs(ps i) ** 2 /
     $        (0.1D1 + (alpha *dabs(psi)) ** n) - (thetas - thetar) /
     $        (0.1D1 + (alpha *dabs(psi)) ** n) ** m * m * (alpha *
     $        dabs(psi)) ** n * n *0.d0 / dabs(psi) / (0.1D1 + (alpha *
     $        dabs(psi)) * * n) +(thetas - thetar) / (0.1D1 + (alpha *
     $        dabs(psi)) ** n) ** m *m * (alpha * dabs(psi)) ** n * n *
     $        dabs(psi) / psi ** 2 /dabs(psi) ** 2 / (0.1D1 + (alpha *
     $        dabs(psi)) ** n) + (thetas- thetar) / (0 .1D1 + (alpha *
     $        dabs(psi)) ** n) ** m * m *((alpha * dabs(psi)) ** n) ** 2
     $        * n ** 2 * dabs(psi) / psi **2 / dabs(psi) ** 2 / (0.1D1 +
     $        (alpha * dabs(psi)) ** n) ** 2
      else
         ddvolc = 0.d0
      endif
      return
      end
