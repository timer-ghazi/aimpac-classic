! Generated from gridv.f
! Contains: RDPSI, MAKNAME, MAKNAME_DIRECT, NUMBER


      SUBROUTINE RDPSI
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  RDPSI READS AIMPAC WAVEFUNCTION FILE.  
C
C   INPUT CONSISTS OF;
C     MODE   - WAVEFUNCTION TYPE (SLATER OR GAUSSIAN)
C     NMO    - NO. OF MOLECULAR ORBITALS
C     NPRIMS - NO. OF (PRIMITIVE) BASIS FUNCTIONS
C     NCENT  - NO. OF NUCLEI
C
C  THEN FOR EACH NUCLEUS;
C     XC/YC/ZC   - COORDINATES
C     CHARG      - ATOMIC NUMBER
C
C  AND FOR EACH BASIS FUNCTION;
C     ICENT  - THE NO. OF THE NUCLEUS ON WHICH IT IS CENTERED
C     ITYPE  - FUNCTION TYPE (SEE ARRAYS LABELS AND LABELG)
C     EX     - EXPONENT
C
C  AND FOR EACH MOLECULAR ORBITAL;
C     PO    - OCCUPATION NUMBER
C     EORB  - ORBITAL ENERGY
C     CO    - COEFFICIENTS OF PRIMITIVE BASIS FUNCTIONS  
C             INCLUDING ALL NORMALIZATION AND CONTRACTION COEFFICIENTS
C
C  FOLOWED BY;
C     TOTE  - MOLECULAR ENERGY
C     GAMMA - VIRIAL RATIO (V/T)
C
      CHARACTER*80 WFNTTL,JOBTTL,LINE
      CHARACTER*8 ATNAM
      PARAMETER(MCENT=200,MMO=500,MPRIMS=2000,MPTS=2000,NTYPE=20)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),NMO
      COMMON /PRIMS/ DIV(MPRIMS),COO(MPRIMS,MMO),EXX(MPRIMS),
     +ICT(MPRIMS),SUM(MPRIMS),ITP(NTYPE),NPRIMS
      COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(MCENT),NAT
      COMMON /VALUES/ THRESH1,THRESH2,GAMMA,TOTE
      COMMON /UNITS/ INPT,IOUT,IWFN,IDBG
      DIMENSION ITYPE(MPRIMS),CO(MPRIMS,MMO),ICENT(MPRIMS),EX(MPRIMS)
C
C    THE ITYPE ARRAY REPRESENTS THE FOLLOWING GAUSSIAN ORBITAL TYPES:
C     S
C     PX, PY, PZ 
C     DXX, DYY, DZZ, DXY, DXZ, DYZ
C     FXXX, FYYY, FZZZ, FXXY, FXXZ, FYYZ, FXYY, FXZZ, FYZZ, FXYZ
C
C    READ IN WAVEFUNCTION TITLE
C
      READ(IWFN,101) WFNTTL
C
C    READ IN MODE, NUMBER OF MOLECULAR ORBITALS, NUMBER OF PRIMITIVES,
C    AND NUMBER OF ATOMIC CENTERS
C
      READ(IWFN,102) MODE,NMO,NPRIMS,NCENT
C
      IF (NMO .GT. MMO) STOP 'Too many MOs'
      IF (NPRIMS .GT. MPRIMS) STOP 'Too many primitives'
      IF (NCENT .GT. MCENT) STOP 'Too many nuclei'
C    READ IN ATOM NAMES, CARTESIAN COORDINATES, AND NUMBER OF
C    ELECTRONS PER ATOM
C
      DO 100 I = 1,NCENT
        READ(IWFN,103) ATNAM(I),J,XC(J),YC(J),ZC(J),CHARG(J)
100   CONTINUE
C
C    READ IN FUNCTION CENTERS, FUNCTION TYPES, AND EXPONENTS
C
      READ(IWFN,104) (ICENT(I),I=1,NPRIMS)
      READ(IWFN,104) (ITYPE(I),I=1,NPRIMS)
      READ(IWFN,105) (EX(I),I=1,NPRIMS)
C
C    LOOP OVER MOLECULAR ORBITALS READING ORBITAL OCCUPANICES,
C    ORBITAL ENERGIES, AND COEFFICIENTS
C
      DO 110 I = 1,NMO
        READ(IWFN,106) PO(I),EORB(I)
        READ(IWFN,107) (CO(J,I),J=1,NPRIMS)
110   CONTINUE
C
      READ(IWFN,101) LINE
C
C    READ IN ENERGY AND -(V/T), THE VIRIAL RATIO
C
      READ(IWFN,109) TOTE,GAMMA
C
      N=0
      DO 160 K=1,NTYPE
      DO 170 J=1,NPRIMS
      IF(ITYPE(J).EQ.K) THEN
      N=N+1
      EXX(N)=EX(J)
      ICT(N)=ICENT(J)
      SUM(N)=-2.0*EX(J)
      DIV(N)=1./SUM(N)
      DO 180 L=1,NMO
      COO(N,L)=CO(J,L)
180   CONTINUE
      ENDIF
170   CONTINUE
      ITP(K)=N
160   CONTINUE 
C
      RETURN
C    FORMATS
C
101   FORMAT (A80)
102   FORMAT (4X,A4,12X,3(I3,17X))
103   FORMAT (A8,11X,I3,2X,3F12.8,10X,F5.1)
104   FORMAT (20X,20I3)
105   FORMAT (10X,5E14.7)
106   FORMAT (35X,F12.8,15X,F12.8)
107   FORMAT (5E16.8)
109   FORMAT (17X,F20.12,18X,F13.8)
      END

      SUBROUTINE MAKNAME(I,STRING,L,EXT)
      CHARACTER*(*) STRING,EXT
      INTEGER I,J,EXTLEN,DOT_POS
      CALL GETARG(I,STRING)
      J = LEN(STRING)
      EXTLEN = LEN(EXT)
      
C     Find the length of the actual argument (before spaces)
      L = 0
      DO 10 N = 1,J
        IF(STRING(N:N) .EQ. ' ') THEN
          L = N - 1
          GOTO 20
        ENDIF
10    CONTINUE
      L = J
      
20    CONTINUE
      IF (L .EQ. 0) RETURN
      
C     Check if filename already has the desired extension
      IF (L .GE. EXTLEN) THEN
        IF (STRING(L-EXTLEN+1:L) .EQ. EXT) THEN
C         Already has correct extension, just trim spaces
          STRING = STRING(1:L)
          RETURN
        ENDIF
      ENDIF
      
C     Check if filename has any extension (contains a dot)
      DOT_POS = 0
      DO 30 N = L,1,-1
        IF (STRING(N:N) .EQ. '.') THEN
          DOT_POS = N
          GOTO 40
        ENDIF
30    CONTINUE
      
40    CONTINUE
      IF (DOT_POS .GT. 0) THEN
C       Has an extension, keep the full filename as-is
        STRING = STRING(1:L)
      ELSE
C       No extension, add the default one
        STRING = STRING(1:L)//EXT
        L = L + EXTLEN
      ENDIF
      
      RETURN
      END

      SUBROUTINE MAKNAME_DIRECT(STRING,L,EXT)
      CHARACTER*(*) STRING,EXT
      INTEGER J,EXTLEN,DOT_POS
      J = LEN(STRING)
      EXTLEN = LEN(EXT)
      
C     Find the length of the actual argument (before spaces)
      L = 0
      DO 10 N = 1,J
        IF(STRING(N:N) .EQ. ' ') THEN
          L = N - 1
          GOTO 20
        ENDIF
10    CONTINUE
      L = J
      
20    CONTINUE
      IF (L .EQ. 0) RETURN
      
C     Check if filename already has the desired extension
      IF (L .GE. EXTLEN) THEN
        IF (STRING(L-EXTLEN+1:L) .EQ. EXT) THEN
C         Already has correct extension, just trim spaces
          STRING = STRING(1:L)
          RETURN
        ENDIF
      ENDIF
      
C     Check if filename has any extension (contains a dot)
      DOT_POS = 0
      DO 30 N = L,1,-1
        IF (STRING(N:N) .EQ. '.') THEN
          DOT_POS = N
          GOTO 40
        ENDIF
30    CONTINUE
      
40    CONTINUE
      IF (DOT_POS .GT. 0) THEN
C       Replace existing extension with new one
        STRING = STRING(1:DOT_POS-1)//EXT
        L = DOT_POS - 1 + EXTLEN
      ELSE
C       No extension, add the default one
        STRING = STRING(1:L)//EXT
        L = L + EXTLEN
      ENDIF
      
      RETURN
      END

        FUNCTION        NUMBER  (LINE, LPST, NUM, DEC)
C
C CONVERTS A CHARACTER STRING OF NUMBERS INTO ACTUAL NUMBERS EITHER
C INTEGERS OR DECIMAL MAY BE READ.
C NUMBER = 1 IF ALL THE REMAINING CHARACTERS ARE BLANK
C        = 2 IF CHARACTERS ARE NOT RECOGNISED AS A NUMBER, LPST IS RESET
C SKK ================================================================== SKK

        DOUBLE PRECISION DEC, TEN
        CHARACTER*(*)   LINE
        CHARACTER       BLANK, COMMA, DOT, MINUS, L
        CHARACTER       CTEN(0:9)
        DATA    CTEN    /'0','1','2','3','4','5','6','7','8','9'/
        PARAMETER       (BLANK = ' ', COMMA = ',')
        PARAMETER       (DOT   = '.', MINUS = '-')
        INTEGER         ITEN
        PARAMETER       (ITEN = 10, TEN = 10.0D0)
        NUM     = 0
        DEC     = 0.0D0
        NP      = 0
        ND      = 0
        MS      = 0
        NUMBER  = 0
        LPEND   = LEN (LINE)
5       IF (LINE(LPST:LPST) .EQ. BLANK) THEN
                LPST    = LPST + 1
                IF (LPST .GT. LPEND) THEN
                        NUMBER  = 1
                        RETURN
                END IF
                GOTO 5
        END IF
        LBEFOR  = LPST

        DO 1 I  = LBEFOR, LPEND
        LPST    = I
        L       = LINE(I:I)
        IF (L .EQ. BLANK .OR. L .EQ. COMMA) THEN
                GOTO 2
        ELSE IF (L .EQ. MINUS) THEN
                MS      = 1
                GOTO 1
        ELSE IF (L .EQ. DOT) THEN
                NP      = 1
                GOTO 1
        END IF
        DO 3 J  = 0, 9
        IF (L .EQ. CTEN(J)) THEN
                N       = J
                GOTO 4
        END IF
3       CONTINUE
        NUMBER  = 2
        LPST    = LBEFOR
        RETURN

4       IF (NP .EQ. 1) THEN
                ND      = ND + 1
                DEC     = DEC + DFLOAT(N)/TEN**ND
        ELSE
                NUM     = NUM*ITEN + N
        END IF
1       CONTINUE

2       DEC     = DFLOAT(NUM) + DEC
        IF (MS .EQ. 0) RETURN
        DEC     = -DEC
        NUM     = -NUM
        RETURN
        END
