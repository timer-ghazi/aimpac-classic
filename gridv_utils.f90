FUNCTION        NUMBER  (LINE, LPST, NUM, DEC)
!
! CONVERTS A CHARACTER STRING OF NUMBERS INTO ACTUAL NUMBERS EITHER
! INTEGERS OR DECIMAL MAY BE READ.
! NUMBER = 1 IF ALL THE REMAINING CHARACTERS ARE BLANK
!        = 2 IF CHARACTERS ARE NOT RECOGNISED AS A NUMBER, LPST IS RESET
! SKK ================================================================== SKK

   IMPLICIT NONE
   DOUBLE PRECISION DEC, TEN
   CHARACTER*(*)   LINE
   CHARACTER       BLANK, COMMA, DOT, MINUS, L
   CHARACTER       CTEN(0:9)
   INTEGER         ITEN, NUM, LPST, NUMBER, NP, ND, MS, LPEND, LBEFOR, I, J, N
   DATA    CTEN    /'0','1','2','3','4','5','6','7','8','9'/
   PARAMETER       (BLANK = ' ', COMMA = ',')
   PARAMETER       (DOT   = '.', MINUS = '-')
   PARAMETER       (ITEN = 10, TEN = 10.0D0)
   NUM     = 0
   DEC     = 0.0D0
   NP      = 0
   ND      = 0
   MS      = 0
   NUMBER  = 0
   LPEND   = LEN (LINE)
5  IF (LINE(LPST:LPST) .EQ. BLANK) THEN
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
3     CONTINUE
      NUMBER  = 2
      LPST    = LBEFOR
      RETURN

4     IF (NP .EQ. 1) THEN
         ND      = ND + 1
         DEC     = DEC + DFLOAT(N)/TEN**ND
      ELSE
         NUM     = NUM*ITEN + N
      END IF
1  CONTINUE

2  DEC     = DFLOAT(NUM) + DEC
   IF (MS .EQ. 0) RETURN
   DEC     = -DEC
   NUM     = -NUM
   RETURN
END