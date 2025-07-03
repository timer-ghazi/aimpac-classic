INTEGER FUNCTION NUMBER  (LINE, LPST, NUM, DEC)
!
! CONVERTS A CHARACTER STRING OF NUMBERS INTO ACTUAL NUMBERS EITHER
! INTEGERS OR DECIMAL MAY BE READ.
! NUMBER = 1 IF ALL THE REMAINING CHARACTERS ARE BLANK
!        = 2 IF CHARACTERS ARE NOT RECOGNISED AS A NUMBER, LPST IS RESET
! SKK ================================================================== SKK

   implicit none
   
   ! Arguments
   character(len=*), intent(in) :: LINE
   integer, intent(inout) :: LPST
   integer, intent(out) :: NUM
   real(8), intent(out) :: DEC
   
   ! Local variables
   real(8) :: TEN
   character :: BLANK, COMMA, DOT, MINUS, L
   character :: CTEN(0:9)
   integer :: ITEN, NP, ND, MS, LPEND, LBEFOR, I, J, N
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
   
   ! Skip leading blanks
   do while (LPST <= LPEND)
      if (LINE(LPST:LPST) /= BLANK) exit
      LPST = LPST + 1
   end do
   
   if (LPST > LPEND) then
      NUMBER = 1
      return
   end if
   
   LBEFOR = LPST

   do i = LBEFOR, LPEND
      LPST = i
      L = LINE(i:i)
      if (L == BLANK .or. L == COMMA) then
         exit
      else if (L == MINUS) then
         MS = 1
         cycle
      else if (L == DOT) then
         NP = 1
         cycle
      end if
      
      ! Search for digit
      do j = 0, 9
         if (L == CTEN(j)) then
            N = j
            exit
         end if
      end do
      
      ! If no digit found, return error
      if (j > 9) then
         NUMBER = 2
         LPST = LBEFOR
         return
      end if
      
      if (NP == 1) then
         ND = ND + 1
         DEC = DEC + DFLOAT(N)/TEN**ND
      else
         NUM = NUM*ITEN + N
      end if
   end do

   DEC = DFLOAT(NUM) + DEC
   if (MS == 0) return
   DEC = -DEC
   NUM = -NUM
   return
END