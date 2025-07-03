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
   do i = 1, ncent
      read(iwfn,103) atnam(i),j,xc(j),yc(j),zc(j),charg(j)
   end do
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
   do i = 1, nmo
      read(iwfn,106) po(i),eorb(i)
      read(iwfn,107) (co(j,i),j=1,nprims)
   end do
!
   read(IWFN,101) LINE
!
!    READ IN ENERGY AND -(V/T), THE VIRIAL RATIO
!
   read(IWFN,109) TOTE,GAMMA
!
   n = 0
   do k = 1, ntype
      do j = 1, nprims
         if (itype(j) == k) then
            n = n + 1
            exx(n) = ex(j)
            ict(n) = icent(j)
            sum(n) = -2.0*ex(j)
            div(n) = 1./sum(n)
            do l = 1, nmo
               coo(n,l) = co(j,l)
            end do
         end if
      end do
      itp(k) = n
   end do
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
   l = 0
   do n = 1, j
      if (string(n:n) == ' ') then
         l = n - 1
         exit
      end if
   end do
   if (l == 0) l = j
   if (l == 0) return

!     Check if filename already has the desired extension
   if (l >= extlen) then
      if (string(l-extlen+1:l) == ext) then
!         Already has correct extension, just trim spaces
         string = string(1:l)
         return
      end if
   end if

!     Check if filename has any extension (contains a dot)
   DOT_POS = 0
   do n = l, 1, -1
      if (string(n:n) == '.') then
         dot_pos = n
         exit
      end if
   end do
   if (dot_pos > 0) then
!       Has an extension, keep the full filename as-is
      string = string(1:l)
   else
!       No extension, add the default one
      string = string(1:l)//ext
      l = l + extlen
   end if

   RETURN
END

SUBROUTINE MAKNAME_DIRECT(STRING,L,EXT)
   IMPLICIT NONE
   CHARACTER*(*) STRING,EXT
   INTEGER J,EXTLEN,DOT_POS,L,N
   J = LEN(STRING)
   EXTLEN = LEN(EXT)

!     Find the length of the actual argument (before spaces)
   l = 0
   do n = 1, j
      if (string(n:n) == ' ') then
         l = n - 1
         exit
      end if
   end do
   if (l == 0) l = j
   if (l == 0) return

!     Check if filename already has the desired extension
   if (l >= extlen) then
      if (string(l-extlen+1:l) == ext) then
!         Already has correct extension, just trim spaces
         string = string(1:l)
         return
      end if
   end if

!     Check if filename has any extension (contains a dot)
   DOT_POS = 0
   do n = l, 1, -1
      if (string(n:n) == '.') then
         dot_pos = n
         exit
      end if
   end do
   if (dot_pos > 0) then
!       Replace existing extension with new one
      string = string(1:dot_pos-1)//ext
      l = dot_pos - 1 + extlen
   else
!       No extension, add the default one
      string = string(1:l)//ext
      l = l + extlen
   end if

   RETURN
END