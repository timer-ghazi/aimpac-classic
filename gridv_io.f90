SUBROUTINE RDPSI
!
   USE molecular_data
   USE basis_data
   USE string_data
   USE computation_data
   USE io_units
   implicit none
!
!  RDPSI READS AIMPAC WAVEFUNCTION FILE.
!
!   INPUT CONSISTS OF;
!     MODE   - WAVEFUNCTION TYPE (SLATER OR GAUSSIAN)
!     NMO    - NO. OF MOLECULAR ORBITALS
!     NPRIMS - NO. OF (PRIMITIVE) BASIS FUNCTIONS
!     NCENT  - NO. OF NUCLEI
!
!  THEN FOR EACH NUCLEUS;
!     XC/YC/ZC   - COORDINATES
!     CHARG      - ATOMIC NUMBER
!
!  AND FOR EACH BASIS FUNCTION;
!     ICENT  - THE NO. OF THE NUCLEUS ON WHICH IT IS CENTERED
!     ITYPE  - FUNCTION TYPE (SEE ARRAYS LABELS AND LABELG)
!     EX     - EXPONENT
!
!  AND FOR EACH MOLECULAR ORBITAL;
!     PO    - OCCUPATION NUMBER
!     EORB  - ORBITAL ENERGY
!     CO    - COEFFICIENTS OF PRIMITIVE BASIS FUNCTIONS
!             INCLUDING ALL NORMALIZATION AND CONTRACTION COEFFICIENTS
!
!  FOLOWED BY;
!     TOTE  - MOLECULAR ENERGY
!     GAMMA - VIRIAL RATIO (V/T)
!
   ! Local variables
   character(len=80) :: LINE
   integer :: ITYPE(2000), ICENT(2000)
   real(8) :: CO(2000,500), EX(2000)
   integer :: I, J, K, L, M, N, NDUM
   character(len=4) :: MODE
!
!    THE ITYPE ARRAY REPRESENTS THE FOLLOWING GAUSSIAN ORBITAL TYPES:
!     S
!     PX, PY, PZ
!     DXX, DYY, DZZ, DXY, DXZ, DYZ
!     FXXX, FYYY, FZZZ, FXXY, FXXZ, FYYZ, FXYY, FXZZ, FYZZ, FXYZ
!
!    read IN WAVEFUNCTION TITLE
!
   read(IWFN,101) WFNTTL
!
!    READ IN MODE, NUMBER OF MOLECULAR ORBITALS, NUMBER OF PRIMITIVES,
!    AND NUMBER OF ATOMIC CENTERS
!
   read(IWFN,102) MODE,NMO,NPRIMS,NCENT
!
   IF (NMO .GT. MMO) STOP 'Too many MOs'
   IF (NPRIMS .GT. MPRIMS) STOP 'Too many primitives'
   IF (NCENT .GT. MCENT) STOP 'Too many nuclei'
!    READ IN ATOM NAMES, CARTESIAN COORDINATES, AND NUMBER OF
!    ELECTRONS PER ATOM
!
   DO 100 I = 1,NCENT
      read(IWFN,103) ATNAM(I),J,XC(J),YC(J),ZC(J),CHARG(J)
100 CONTINUE
!
!    READ IN FUNCTION CENTERS, FUNCTION TYPES, AND EXPONENTS
!
   read(IWFN,104) (ICENT(I),I=1,NPRIMS)
   read(IWFN,104) (ITYPE(I),I=1,NPRIMS)
   read(IWFN,105) (EX(I),I=1,NPRIMS)
!
!    LOOP OVER MOLECULAR ORBITALS READING ORBITAL OCCUPANICES,
!    ORBITAL ENERGIES, AND COEFFICIENTS
!
   DO 110 I = 1,NMO
      read(IWFN,106) PO(I),EORB(I)
      read(IWFN,107) (CO(J,I),J=1,NPRIMS)
110 CONTINUE
!
   read(IWFN,101) LINE
!
!    READ IN ENERGY AND -(V/T), THE VIRIAL RATIO
!
   read(IWFN,109) TOTE,GAMMA
!
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
180         CONTINUE
         ENDIF
170   CONTINUE
      ITP(K)=N
160 CONTINUE
!
   RETURN
!    FORMATS
!
101 FORMAT (A80)
102 FORMAT (4X,A4,12X,3(I3,17X))
103 FORMAT (A8,11X,I3,2X,3F12.8,10X,F5.1)
104 FORMAT (20X,20I3)
105 FORMAT (10X,5E14.7)
106 FORMAT (35X,F12.8,15X,F12.8)
107 FORMAT (5E16.8)
109 FORMAT (17X,F20.12,18X,F13.8)
END

SUBROUTINE MAKNAME(I,STRING,L,EXT)
   IMPLICIT NONE
   CHARACTER*(*) STRING,EXT
   INTEGER I,J,EXTLEN,DOT_POS,L,N
   CALL GETARG(I,STRING)
   J = LEN(STRING)
   EXTLEN = LEN(EXT)

!     Find the length of the actual argument (before spaces)
   L = 0
   DO 10 N = 1,J
      IF(STRING(N:N) .EQ. ' ') THEN
         L = N - 1
         GOTO 20
      ENDIF
10 CONTINUE
   L = J

20 CONTINUE
   IF (L .EQ. 0) RETURN

!     Check if filename already has the desired extension
   IF (L .GE. EXTLEN) THEN
      IF (STRING(L-EXTLEN+1:L) .EQ. EXT) THEN
!         Already has correct extension, just trim spaces
         STRING = STRING(1:L)
         RETURN
      ENDIF
   ENDIF

!     Check if filename has any extension (contains a dot)
   DOT_POS = 0
   DO 30 N = L,1,-1
      IF (STRING(N:N) .EQ. '.') THEN
         DOT_POS = N
         GOTO 40
      ENDIF
30 CONTINUE

40 CONTINUE
   IF (DOT_POS .GT. 0) THEN
!       Has an extension, keep the full filename as-is
      STRING = STRING(1:L)
   ELSE
!       No extension, add the default one
      STRING = STRING(1:L)//EXT
      L = L + EXTLEN
   ENDIF

   RETURN
END

SUBROUTINE MAKNAME_DIRECT(STRING,L,EXT)
   IMPLICIT NONE
   CHARACTER*(*) STRING,EXT
   INTEGER J,EXTLEN,DOT_POS,L,N
   J = LEN(STRING)
   EXTLEN = LEN(EXT)

!     Find the length of the actual argument (before spaces)
   L = 0
   DO 10 N = 1,J
      IF(STRING(N:N) .EQ. ' ') THEN
         L = N - 1
         GOTO 20
      ENDIF
10 CONTINUE
   L = J

20 CONTINUE
   IF (L .EQ. 0) RETURN

!     Check if filename already has the desired extension
   IF (L .GE. EXTLEN) THEN
      IF (STRING(L-EXTLEN+1:L) .EQ. EXT) THEN
!         Already has correct extension, just trim spaces
         STRING = STRING(1:L)
         RETURN
      ENDIF
   ENDIF

!     Check if filename has any extension (contains a dot)
   DOT_POS = 0
   DO 30 N = L,1,-1
      IF (STRING(N:N) .EQ. '.') THEN
         DOT_POS = N
         GOTO 40
      ENDIF
30 CONTINUE

40 CONTINUE
   IF (DOT_POS .GT. 0) THEN
!       Replace existing extension with new one
      STRING = STRING(1:DOT_POS-1)//EXT
      L = DOT_POS - 1 + EXTLEN
   ELSE
!       No extension, add the default one
      STRING = STRING(1:L)//EXT
      L = L + EXTLEN
   ENDIF

   RETURN
END