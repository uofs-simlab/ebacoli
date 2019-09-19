      SUBROUTINE BSPCNDMAX(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,
     *                  NBLOKS,BOTBLK,NRWBOT,EST,V,ISIGN,WORK,
     *                  PIVOT,IFLAG)

C
C***************************************************************
C
C  THIS PROGRAM  COMPUTES AN ESTIMATE OF THE CONDITION NUMBER 
C  IN THE MAX NORM OF AN ALMOST BLOCK DIAGONAL MATRIX OF THE FORM:
C  
C
C               TOPBLK
C               ARRAY(1)
C                     ARRAY(2)
C                          .
C                             .
C                                .
C                                   .
C                                    ARRAY(NBLOKS)
C                                           BOTBLK
C
C  WHERE
C           TOPBLK IS  NRWTOP  BY NOVRLP
C           ARRAY(K), K=1,NBLOKS, ARE NRWBLK BY NRWBLK+NOVRLP
C           BOTBLK IS NRWBOT BY NOVRLP,
C  AND
C           NOVRLP = NRWTOP + NRWBOT
C  WITH
C           NOVRLP.LE.NRWBLK .
C
C  THE LINEAR SYSTEM IS OF ORDER  N = NBLOKS*NRWBLK + NOVRLP.
C
C  THE METHOD IMPLEMENTED IS BASED ON GAUSS ELIMINATION WITH
C  ALTERNATE ROW AND COLUMN ELIMINATION WITH PARTIAL PIVOTING,
C  WHICH PRODUCES A STABLE DECOMPOSITION OF THE MATRIX  A
C  WITHOUT INTRODUCING FILL-IN. SEE COLROW DOCUMENTATION FOR
C  SAMPLE DRIVING PROGRAM ETC.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C               ISIGN  - INTEGER(N)
C                         WORK SPACE
C
C               V      - DOUBLE PRECISION (N), = INV(A)*W, WHERE
C                        NORM(V)/NORM(W) ESTIMATES NORM OF INV(A).
C
C               WORK   - DOUBLE PRECISION(N), WORK SPACE.
C
C       *** ON RETURN  ...
C
C               EST    - ESTIMATED MAX NORM CONDITION NUMBER.
C               V      - APPROXIMATE NULL VECTOR IF A IS SINGULAR.
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  AUXILIARY PROGRAMS  *****
C
C       CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
C    *     BOTBLK,NRWBOT,PIVOT,IFLAG)
C            - DECOMPOSES THE MATRIX  A  USING MODIFIED
C              ALTERNATE ROW AND COLUMN ELIMINATON WITH
C              PARTIAL PIVOTING, AND IS USED FOR THIS
C              PURPOSE IN C O L R O W.
C              THE ARGUMENTS ARE AS IN C O L R O W.
C
C
C       CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
C    *     BOTBLK,NRWBOT,PIVOT,B,JOB)
C            - SOLVES THE SYSTEM A*X = B ONCE A IS DECOMPOSED.
C              THE ARGUMENTS ARE ALL AS IN C O L R O W.
C
C       ABDNRMMAX(NBLOKS,NTOP,NBOT,NOVRLP,
C    *         NRWBLK,NCLBLK,TOP,A,BOT)
C              COMPUTES THE MAX NORM OF THE ABD MATRIX TOP, A, BOT.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*),
     *          IFLAG,ISIGN(*)
        DOUBLE PRECISION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),EST,V(*),WORK(*)
        DOUBLE PRECISION ABDNRMMAX,NORM
        INTEGER KASE,ISOLVE,NOUT
        NOUT = 6
C
C       FIRST, COMPUTE THE NORM OF THE MATRIX ITSELF:
C
        NORM = ABDNRMMAX(NBLOKS,NRWTOP,NRWBOT,NOVRLP,NRWBLK,NCLBLK,
     *                TOPBLK,ARRAY,BOTBLK)
C
C       THEN, DO THE FACTORIZATION USING CRDCMP:
C
        CALL CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,IFLAG)

        IF ( IFLAG .NE. 0 ) RETURN
C
C       MAKING JOB=2-KASE GIVES DONEST A WHEN IT ASKS FOR A TRANSPOSE
C       AND A TRANSPOSE WHEN IT ASKS FOR A.  THIS CAUSES DONEST TO RETURN
C       THE MAX NORM OF A INVERSE.  THIS METHOD TAKES ADVANTAGE OF THE 
C       IDEA THAT THE ONE NORM OF A TRANSPOSE IS EQUAL TO THE MAX NORM OF
C       A.
C
        EST = 0.0D0
        ISOLVE = 0
        KASE = 0
  55    CALL DONEST(N,V,WORK,ISIGN,EST,KASE)
        IF (KASE .NE. 0) THEN
           ISOLVE = ISOLVE+1
           CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,
     *                 WORK,2-KASE)

           GOTO 55
        END IF
        EST = EST*NORM
        RETURN
        END
C
C==================================================================
C
        DOUBLE PRECISION FUNCTION ABDNRMMAX(NBLOKS,NTOP,NBOT,NOVRLP,
     *                            NRWBLK,NCLBLK,TOP,A,BOT)

C******************************************************************
C       ABDNRMMAX IS USED WITH DONEST TO ESTIMATE THE MAX NORM
C       OF AN ALMOST BLOCK DIAGONAL MATRIX LIKE THE ONES HANDLED BY       
C       COLROW.
C******************************************************************

        INTRINSIC DABS
        INTEGER NBLOKS,NTOP,NBOT,NOVRLP,NRWBLK,NCLBLK
        DOUBLE PRECISION TOP(NTOP,*),A(NRWBLK,NCLBLK,*),BOT(NBOT,*)
        INTEGER I,J,K
        DOUBLE PRECISION TEMP,RWSUM,CELL
        DOUBLE PRECISION DABS
        TEMP = 0.0D0
C
C       FIRST, CHECK THE TOP BLOCK:
C
        DO 10 I=1,NTOP
          RWSUM=0.0D0
          DO 20 J=1,NOVRLP
            CELL = TOP(I,J)
            RWSUM=RWSUM+DABS(CELL)
20        CONTINUE
          IF (RWSUM .GT. TEMP) THEN
            TEMP = RWSUM
          ENDIF
10      CONTINUE

C       NEXT CHECK THE BLOCKS IN THE ARRAY

        DO 30 I=1,NBLOKS
          DO 40 J=1,NRWBLK
            RWSUM = 0.0D0
            DO 50 K=1,NCLBLK
              CELL = A(J,K,I)
              RWSUM=RWSUM+DABS(CELL)
50          CONTINUE
            IF (RWSUM .GT. TEMP) THEN
              TEMP = RWSUM
            ENDIF
40        CONTINUE
30      CONTINUE

C       LAST, CHECK THE BOTTOM BLOCK:
C
        DO 60 I=1,NBOT
          RWSUM=0.0D0
          DO 70 J=1,NOVRLP
            CELL = BOT(I,J)
            RWSUM=RWSUM+DABS(CELL)
70        CONTINUE
          IF (RWSUM .GT. TEMP) THEN
            TEMP = RWSUM
          ENDIF
60      CONTINUE                                                
C

        ABDNRMMAX = TEMP
        RETURN
        END
        
      SUBROUTINE DONEST (N, V, X, ISGN, EST, KASE)
      INTEGER N, ISGN(N), KASE
      DOUBLE PRECISION V(N), X(N), EST

C
C     DONEST ESTIMATES THE 1-NORM OF A SQUARE, DOUBLE PRECISION MATRIX  A.
C     REVERSE COMMUNICATION IS USED FOR EVALUATING
C     MATRIX-VECTOR PRODUCTS. 
C
C     ON ENTRY
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX.  N .GE. 1.
C
C        ISGN    INTEGER(N)
C                USED AS WORKSPACE.
C
C        KASE    INTEGER
C                = 0.
C
C     ON INTERMEDIATE RETURNS 
C
C        KASE    = 1 OR 2.
C
C        X       DOUBLE PRECISION(N)
C                MUST BE OVERWRITTEN BY 
C
C                     A*X,             IF KASE=1, 
C                     TRANSPOSE(A)*X,  IF KASE=2, 
C
C                AND DONEST MUST BE RE-CALLED, WITH ALL THE OTHER
C                PARAMETERS UNCHANGED.
C
C     ON FINAL RETURN
C
C        KASE    = 0.
C
C        EST     DOUBLE PRECISION
C                CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C        V       DOUBLE PRECISION(N)
C                = A*W,   WHERE  EST = NORM(V)/NORM(W)
C                         (W  IS NOT RETURNED).
C
C     THIS VERSION DATED MARCH 16, 1988.
C     NICK HIGHAM, UNIVERSITY OF MANCHESTER.
C
C     MODIFIED FOR DOUBLE PRECISION ON JUNE 11, 1996.
C
C     REFERENCE
C     N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C     THE ONE-NORM OF A REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C     TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C     UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C     SUBROUTINES AND FUNCTIONS
C     BLAS     IDAMAX, DASUM, DCOPY
C     GENERIC  ABS, NINT, FLOAT, SIGN
C
        INTRINSIC FLOAT
        DOUBLE PRECISION FLOAT

        INTRINSIC ABS
        DOUBLE PRECISION ABS

        INTRINSIC SIGN
        DOUBLE PRECISION SIGN


      DOUBLE PRECISION DASUM
      INTEGER IDAMAX

      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO = 0.0D0)
      PARAMETER (ONE = 1.0D0)
      PARAMETER (TWO = 2.0D0)
C
C     INTERNAL VARIABLES
      INTEGER I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION ALTSGN, ESTOLD, TEMP
C
      SAVE
C
      IF (KASE .EQ. 0) THEN
         DO 10,I = 1,N
            X(I) = ONE/FLOAT(N)
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      ENDIF
C
      GOTO (100, 200, 300, 400, 500) JUMP
C
C     ................ ENTRY   (JUMP = 1)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
C
  100 CONTINUE
      IF (N .EQ. 1) THEN
         V(1) = X(1)
         EST = ABS(V(1))
C        ... QUIT
         GOTO 510
      ENDIF
      EST = DASUM(N,X,1)
C
      DO 110,I = 1,N
         X(I) = SIGN(ONE,X(I))
         ISGN(I) = NINT(X(I)) 
  110 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
C
C     ................ ENTRY   (JUMP = 2)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
  200 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
C
C     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
  220 CONTINUE
      DO 230,I = 1,N
         X(I) = ZERO 
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      RETURN
C
C     ................ ENTRY   (JUMP = 3)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
  300 CONTINUE
      CALL DCOPY(N,X,1,V,1)
      ESTOLD = EST
      EST = DASUM(N,V,1)
      DO 310,I = 1,N
         IF ( NINT( SIGN(ONE,X(I)) ) .NE. ISGN(I) ) GOTO 320
  310 CONTINUE
C     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GOTO 410
C
  320 CONTINUE
C     TEST FOR CYCLING.
      IF (EST .LE. ESTOLD) GOTO 410
C      
      DO 330,I = 1,N
         X(I) = SIGN(ONE,X(I))
         ISGN(I) = NINT(X(I)) 
  330 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
C
C     ................ ENTRY   (JUMP = 4)
C     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF (   (  X(JLAST) .NE. ABS(X(J))  ) .AND.
     +       (ITER .LT. ITMAX)   ) THEN
         ITER = ITER + 1
         GOTO 220
      ENDIF
C
C     ITERATION COMPLETE.  FINAL STAGE. 
C
  410 CONTINUE
      ALTSGN = ONE
      DO 420,I = 1,N
         X(I) = ALTSGN * (ONE + FLOAT(I-1)/FLOAT(N-1))
         ALTSGN = -ALTSGN
  420 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
C
C     ................ ENTRY   (JUMP = 5)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
  500 CONTINUE
      TEMP = TWO*DASUM(N,X,1)/FLOAT(3*N) 
      IF (TEMP. GT. EST) THEN 
         CALL DCOPY(N,X,1,V,1)
         EST = TEMP 
      ENDIF
C
  510 KASE = 0
      RETURN
C
      END
