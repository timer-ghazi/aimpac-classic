! Generated from gridv_calculations.f90
! Contains: GAUS_FULL, GAUS_PSI_ONLY

SUBROUTINE GAUS_FULL(XYZ,PTS,MM)
!
   USE molecular_data
   USE basis_data
   USE work_arrays
   implicit none
!
   ! Arguments
   integer, intent(in) :: PTS, MM
   real(8), intent(in) :: XYZ(2000,3)
   
   ! Local variables
   real(8) :: DX(2000,200), DY(2000,200), DZ(2000,200)
   real(8) :: R2(2000,200), CHI(2000,2000)
   real(8) :: CHIX(2000,2000), CHIY(2000,2000), CHIZ(2000,2000)
   real(8) :: CHID2(2000,2000), CHIMAX(2000)
   real(8), save :: ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0, FIVE = 5.D0
   real(8), save :: SEVEN = 7.D0, THREE = 3.D0, SIX = 6.D0, NINE = 9.D0, CUTOFF = 1.0d-10
   integer :: I, J, K, L, IC, ITY, IK, ICL, IYYY, IZZZ, IXXY, IXXZ, IYYZ, IXYY, IXZZ, IYZZ, IXYZ
   integer :: IS
   real(8) :: EXP, GTMAX, EXPMAX, EXPMIN, R2MIN, R2MAX, RMS, R, R3, R4, R5, R6, R7, R8, R9, R10
   real(8) :: X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, X5, Y5, Z5, X6, Y6, Z6, X7, Y7, Z7
   real(8) :: XX, YY, ZZ, XY, XZ, YZ
   real(8) :: A, B, BXX, BXY, BXZ, BYY, BYZ, BZZ, CHECK, TEMP
!
   DO 110 J = 1,NCENT
      DO 112 I=1,PTS
         DX(I,J) = XYZ(I,1) - XC(J)
         DY(I,J) = XYZ(I,2) - YC(J)
         DZ(I,J) = XYZ(I,3) - ZC(J)
         R2(I,J)= DX(I,J)*DX(I,J)+DY(I,J)*DY(I,J)+DZ(I,J)*DZ(I,J)
112   CONTINUE
110 CONTINUE
!
!       FOR S-TYPE
!
   DO 120 J = 1,ITP(1)
      IS=ict(j)
      DO 122 I=1,PTS
         A=SUM(J)*DEXP(-EXX(J)*R2(I,is))
         CHI(I,J)=A*DIV(J)
         CHIX(I,J)=DX(I,is)*A
         CHIY(I,J)=DY(I,is)*A
         CHIZ(I,J)=DZ(I,is)*A
         CHID2(I,J)=(THREE+SUM(J)*R2(I,is))*A
122   CONTINUE
120 CONTINUE
!
!       FOR Px-TYPE
!
   DO 140 J=ITP(1)+1,ITP(2)
      IS=ict(j)
      DO 142 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         B=DX(I,is)*A*SUM(J)
!
         CHI(I,J)=A*DX(I,is)
         CHIX(I,J)=A+DX(I,is)*B
         CHIY(I,J)=DY(I,is)*B
         CHIZ(I,J)=DZ(I,is)*B
         CHID2(I,J)=(FIVE+SUM(J)*R2(I,is))*B
142   CONTINUE
140 CONTINUE
!
!       FOR Py-TYPE
!
   DO 160 J=ITP(2)+1,ITP(3)
      IS=ict(j)
      DO 162 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         B=DY(I,is)*A*SUM(J)
!
         CHI(I,J)=A*DY(I,is)
         CHIX(I,J)=DX(I,is)*B
         CHIY(I,J)=A+DY(I,is)*B
         CHIZ(I,J)=DZ(I,is)*B
         CHID2(I,J)=(FIVE+SUM(J)*R2(I,is))*B
162   CONTINUE
160 CONTINUE
!
!       FOR Pz-TYPE
!
   DO 180 J=ITP(3)+1,ITP(4)
      IS=ict(j)
      DO 182 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         B=DZ(I,is)*A*SUM(J)
!
         CHI(I,J)=A*DZ(I,is)
         CHIX(I,J)=DX(I,is)*B
         CHIY(I,J)=DY(I,is)*B
         CHIZ(I,J)=A+DZ(I,is)*B
         CHID2(I,J)=(FIVE+SUM(J)*R2(I,is))*B
182   CONTINUE
180 CONTINUE
!
!       FOR Dxx-TYPE
!
   DO 220 J=ITP(4)+1,ITP(5)
      IS=ict(j)
      DO 222 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         B=DX(I,is)*DX(I,is)*A*SUM(J)
!
         CHI(I,J)=B*DIV(J)
         CHIX(I,J)=(TWO*A+B)*DX(I,is)
         CHIY(I,J)=DY(I,is)*B
         CHIZ(I,J)=DZ(I,is)*B
         CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,is))*B
222   CONTINUE
220 CONTINUE
!
!       FOR Dyy-TYPE
!
   DO 240 J=ITP(5)+1,ITP(6)
      IS=ict(j)
      DO 242 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         B=DY(I,is)*DY(I,is)*A*SUM(J)
!
         CHI(I,J)=B*DIV(J)
         CHIX(I,J)=DX(I,is)*B
         CHIY(I,J)=(TWO*A+B)*DY(I,is)
         CHIZ(I,J)=DZ(I,is)*B
         CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,is))*B
242   CONTINUE
240 CONTINUE
!
!       FOR Dzz-TYPE
!
   DO 260 J=ITP(6)+1,ITP(7)
      IS=ict(j)
      DO 262 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         B=DZ(I,is)*DZ(I,is)*A*SUM(J)
!
         CHI(I,J)=B*DIV(J)
         CHIX(I,J)=DX(I,is)*B
         CHIY(I,J)=DY(I,is)*B
         CHIZ(I,J)=(TWO*A+B)*DZ(I,is)
         CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,is))*B
262   CONTINUE
260 CONTINUE
!
!       FOR Dxy-TYPE
!
   DO 280 J=ITP(7)+1,ITP(8)
      IS=ict(j)
      DO 282 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         B=DX(I,is)*DY(I,is)*A*SUM(J)
!
         CHI(I,J)=B*DIV(J)
         CHIX(I,J)=DX(I,is)*B+DY(I,is)*A
         CHIY(I,J)=DY(I,is)*B+DX(I,is)*A
         CHIZ(I,J)=DZ(I,is)*B
         CHID2(I,J)=(SEVEN+SUM(J)*R2(I,is))*B
282   CONTINUE
280 CONTINUE
!
!       FOR Dxz-TYPE
!
   DO 320 J=ITP(8)+1,ITP(9)
      IS=ict(j)
      DO 322 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         B=DX(I,is)*DZ(I,is)*A*SUM(J)
!
         CHI(I,J)=B*DIV(J)
         CHIX(I,J)=DX(I,is)*B+DZ(I,is)*A
         CHIY(I,J)=DY(I,is)*B
         CHIZ(I,J)=DZ(I,is)*B+DX(I,is)*A
         CHID2(I,J)=(SEVEN+SUM(J)*R2(I,is))*B
322   CONTINUE
320 CONTINUE
!
!       FOR Dyz-TYPE
!
   DO 340 J=ITP(9)+1,ITP(10)
      IS=ict(j)
      DO 342 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         B=DY(I,is)*DZ(I,is)*A*SUM(J)
!
         CHI(I,J)=B*DIV(J)
         CHIX(I,J)=DX(I,is)*B
         CHIY(I,J)=DY(I,is)*B+DZ(I,is)*A
         CHIZ(I,J)=DZ(I,is)*B+DY(I,is)*A
         CHID2(I,J)=(SEVEN+SUM(J)*R2(I,is))*B
342   CONTINUE
340 CONTINUE
!
!       FOR Fxxx-TYPE
!
   DO 501 J=ITP(10)+1,ITP(11)
      IS=ict(j)
      DO 502 I=1,PTS

         A=DEXP(-EXX(J)*R2(I,is))
         B=DX(I,is)*DX(I,is)*A
!
         CHI(I,J)=B*DX(I,is)
         CHIX(I,J)=(THREE + SUM(J)*DX(I,is)*DX(I,is))*B
         CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
         CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
         CHID2(I,J)=SIX*A*DX(I,is)+&
         &(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
502   CONTINUE
501 CONTINUE
!
!       FOR Fyyy-TYPE
!
   DO 511 J=ITP(11)+1,ITP(12)
      IS=ict(j)
      DO 512 I=1,PTS
         A=DEXP(-EXX(J)*R2(I,is))
         B=DY(I,is)*DY(I,is)*A
!
         CHI(I,J)=B*DY(I,is)
         CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
         CHIY(I,J)=(THREE + SUM(J)*DY(I,is)*DY(I,is))*B
         CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
         CHID2(I,J)=SIX*A*DY(I,is)+&
         &(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
!
512   CONTINUE
511 CONTINUE
!
!       FOR Fzzz-TYPE
!
   DO 521 J=ITP(12)+1,ITP(13)
      IS=ict(j)
      DO 523 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         B=DZ(I,is)*DZ(I,is)*A
!
         CHI(I,J)=B*DZ(I,is)
         CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
         CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
         CHIZ(I,J)=(THREE + SUM(J)*DZ(I,is)*DZ(I,is))*B
         CHID2(I,J)=SIX*A*DZ(I,is)+&
         &(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
!
523   CONTINUE
521 CONTINUE
!
!       FOR Fxxy-TYPE
!
   DO 531 J=ITP(13)+1,ITP(14)
      IS=ict(j)
      DO 532 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         BXY=DX(I,is)*DY(I,is)*A
         BXX=DX(I,is)*DX(I,is)*A
!
         CHI(I,J)=BXY*DX(I,is)
         CHIX(I,J)=(TWO + SUM(J)*DX(I,is)*DX(I,is))*BXY
         CHIY(I,J)=(ONE + SUM(J)*DY(I,is)*DY(I,is))*BXX
         CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
         CHID2(I,J)=TWO*A*DY(I,is)+&
         &(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
!
532   CONTINUE
531 CONTINUE
!
!       FOR Fxxz-TYPE
!
   DO 541 J=ITP(14)+1,ITP(15)
      IS=ict(j)
      DO 543 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         BXZ=DX(I,is)*DZ(I,is)*A
         BXX=DX(I,is)*DX(I,is)*A
!
         CHI(I,J)=BXZ*DX(I,is)
         CHIX(I,J)=(TWO + SUM(J)*DX(I,is)*DX(I,is))*BXZ
         CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
         CHIZ(I,J)=(ONE + SUM(J)*DZ(I,is)*DZ(I,is))*BXX
         CHID2(I,J)=TWO*A*DZ(I,is)+&
         &(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
!
543   CONTINUE
541 CONTINUE
!
!       FOR Fyyz-TYPE
!
   DO 561 J=ITP(15)+1,ITP(16)
      IS=ict(j)
      DO 563 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         BYZ=DZ(I,is)*DY(I,is)*A
         BYY=DY(I,is)*DY(I,is)*A
!
         CHI(I,J)=BYZ*DY(I,is)
         CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
         CHIY(I,J)=(TWO + SUM(J)*DY(I,is)*DY(I,is))*BYZ
         CHIZ(I,J)=(ONE + SUM(J)*DZ(I,is)*DZ(I,is))*BYY
         CHID2(I,J)=TWO*A*DZ(I,is)+&
         &(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
!
563   CONTINUE
561 CONTINUE
!
!       FOR Fxyy-TYPE
!
   DO 551 J=ITP(16)+1,ITP(17)
      IS=ict(j)
      DO 552 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         BXY=DX(I,is)*DY(I,is)*A
         BYY=DY(I,is)*DY(I,is)*A
!
         CHI(I,J)=BXY*DY(I,is)
         CHIX(I,J)=(ONE + SUM(J)*DX(I,is)*DX(I,is))*BYY
         CHIY(I,J)=(TWO + SUM(J)*DY(I,is)*DY(I,is))*BXY
         CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
         CHID2(I,J)=TWO*A*DX(I,is)+&
         &(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
!
552   CONTINUE
551 CONTINUE
!
!       FOR Fxzz-TYPE
!
   DO 571 J=ITP(17)+1,ITP(18)
      IS=ict(j)
      DO 572 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         BXZ=DZ(I,is)*DX(I,is)*A
         BZZ=DZ(I,is)*DZ(I,is)*A
!
         CHI(I,J)=BXZ*DZ(I,is)
         CHIX(I,J)=(ONE + SUM(J)*DX(I,is)*DX(I,is))*BZZ
         CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
         CHIZ(I,J)=(TWO + SUM(J)*DZ(I,is)*DZ(I,is))*BXZ
         CHID2(I,J)=TWO*A*DX(I,is)+&
         &(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
!
572   CONTINUE
571 CONTINUE
!
!       FOR Fyzz-TYPE
!
   DO 581 J=ITP(18)+1,ITP(19)
      IS=ict(j)
      DO 583 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         BYZ=DZ(I,is)*DY(I,is)*A
         BZZ=DZ(I,is)*DZ(I,is)*A
!
         CHI(I,J)=BYZ*DZ(I,is)
         CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
         CHIY(I,J)=(ONE + SUM(J)*DY(I,is)*DY(I,is))*BZZ
         CHIZ(I,J)=(TWO + SUM(J)*DZ(I,is)*DZ(I,is))*BYZ
         CHID2(I,J)=TWO*A*DY(I,is)+&
         &(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
!
583   CONTINUE
581 CONTINUE
!
!       FOR Fxyz-TYPE
!
   DO 591 J=ITP(19)+1,ITP(20)
      IS=ict(j)
      DO 592 I=1,PTS
!
         A=DEXP(-EXX(J)*R2(I,is))
         BXY=DX(I,is)*DY(I,is)*A
         BYZ=DZ(I,is)*DY(I,is)*A
         BXZ=DX(I,is)*DZ(I,is)*A
!
         CHI(I,J)=DX(I,is)*BYZ
         CHIX(I,J)=(ONE + SUM(J)*DX(I,is)*DX(I,is))*BYZ
         CHIY(I,J)=(ONE + SUM(J)*DY(I,is)*DY(I,is))*BXZ
         CHIZ(I,J)=(ONE + SUM(J)*DZ(I,is)*DZ(I,is))*BXY
         CHID2(I,J)=(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
!
592   CONTINUE
591 CONTINUE
!
   DO 103 J=1,Nprims
      temp=zero
      DO 104 I=1,PTs
         check=DMax1(Dabs(CHi(I,J)),Dabs(CHiX(I,j)),Dabs(CHIY(I,J)),&
         &Dabs(CHIZ(I,j)),Dabs(CHID2(I,J)),Temp)
         If(Check.gt.temp)temp=check
104   Continue
      chimax(j)=temp
103 Continue
!
   DO 105 L = 1,NMO
      DO 107 I=1,PTS
         PSI(I,L) = ZERO
         GX(I,L) = ZERO
         GY(I,L) = ZERO
         GZ(I,L) = ZERO
         D2(I,L) = ZERO
107   CONTINUE
!
      DO 125 J = 1,NPRIMS
         check=dabs(chimax(j)*coo(J,L))
         IF(check.gt.cutoff)THEN
            TEMP=COO(J,L)
            DO 126 I = 1,PTS
               PSI(I,L) = PSI(I,L) + TEMP*CHI(I,J)
               GX(I,L) = GX(I,L) + TEMP*CHIX(I,J)
               GY(I,L) = GY(I,L) + TEMP*CHIY(I,J)
               GZ(I,L) = GZ(I,L) + TEMP*CHIZ(I,J)
               D2(I,L) = D2(I,L) + TEMP*CHID2(I,J)
126         CONTINUE
         ENDIF
125   CONTINUE
105 CONTINUE
!
   RETURN
END

SUBROUTINE GAUS_PSI_ONLY(XYZ,PTS,MM)
!
   USE molecular_data
   USE basis_data
   USE work_arrays
   implicit none
!
   ! Arguments
   integer, intent(in) :: PTS, MM
   real(8), intent(in) :: XYZ(2000,3)
   
   ! Local variables
   real(8) :: DX(2000,200), DY(2000,200), DZ(2000,200)
   real(8) :: R2(2000,200), CHI(2000,2000)
   real(8) :: CHIMAX(2000)
   real(8), save :: ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0, FIVE = 5.D0
   real(8), save :: SEVEN = 7.D0, THREE = 3.D0, SIX = 6.D0, NINE = 9.D0, CUTOFF = 1.0d-10
   integer :: I, J, K, L, IC, ITY, IK, ICL, IYYY, IZZZ, IXXY, IXXZ, IYYZ, IXYY, IXZZ, IYZZ, IXYZ
   integer :: IS
   real(8) :: EXP, GTMAX, EXPMAX, EXPMIN, R2MIN, R2MAX, RMS, R, R3, R4, R5, R6, R7, R8, R9, R10
   real(8) :: X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4, X5, Y5, Z5, X6, Y6, Z6, X7, Y7, Z7
   real(8) :: XX, YY, ZZ, XY, XZ, YZ
   real(8) :: A, B, BXX, BXY, BXZ, BYY, BYZ, BZZ, CHECK, TEMP
!
   DO 410 J = 1,NCENT
      DO 412 I=1,PTS
         DX(I,J) = XYZ(I,1) - XC(J)
         DY(I,J) = XYZ(I,2) - YC(J)
         DZ(I,J) = XYZ(I,3) - ZC(J)
         R2(I,J)= DX(I,J)*DX(I,J)+DY(I,J)*DY(I,J)+DZ(I,J)*DZ(I,J)
412   CONTINUE
410 CONTINUE
!
!       FOR S-TYPE
!
   CALL CALC_S_TYPE_PSI(1, ITP(1), PTS, DX, DY, DZ, R2, CHI)
!
!       FOR Px-TYPE
!
   CALL CALC_P_TYPE_PSI(ITP(1)+1, ITP(2), PTS, 'X', DX, DY, DZ, R2, CHI)
!
!       FOR Py-TYPE
!
   CALL CALC_P_TYPE_PSI(ITP(2)+1, ITP(3), PTS, 'Y', DX, DY, DZ, R2, CHI)
!
!       FOR Pz-TYPE
!
   CALL CALC_P_TYPE_PSI(ITP(3)+1, ITP(4), PTS, 'Z', DX, DY, DZ, R2, CHI)
!
!       FOR Dxx-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(4)+1, ITP(5), PTS, 'XX', DX, DY, DZ, R2, CHI)
!
!       FOR Dyy-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(5)+1, ITP(6), PTS, 'YY', DX, DY, DZ, R2, CHI)
!
!       FOR Dzz-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(6)+1, ITP(7), PTS, 'ZZ', DX, DY, DZ, R2, CHI)
!
!       FOR Dxy-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(7)+1, ITP(8), PTS, 'XY', DX, DY, DZ, R2, CHI)
!
!       FOR Dxz-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(8)+1, ITP(9), PTS, 'XZ', DX, DY, DZ, R2, CHI)
!
!       FOR Dyz-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(9)+1, ITP(10), PTS, 'YZ', DX, DY, DZ, R2, CHI)
!
!       FOR Fxxx-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(10)+1, ITP(11), PTS, 'XXX', DX, DY, DZ, R2, CHI)
!
!       FOR Fyyy-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(11)+1, ITP(12), PTS, 'YYY', DX, DY, DZ, R2, CHI)
!
!       FOR Fzzz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(12)+1, ITP(13), PTS, 'ZZZ', DX, DY, DZ, R2, CHI)
!
!       FOR Fxxy-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(13)+1, ITP(14), PTS, 'XXY', DX, DY, DZ, R2, CHI)
!
!       FOR Fxxz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(14)+1, ITP(15), PTS, 'XXZ', DX, DY, DZ, R2, CHI)
!
!       FOR Fxyy-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(16)+1, ITP(17), PTS, 'XYY', DX, DY, DZ, R2, CHI)
!
!       FOR Fyyz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(15)+1, ITP(16), PTS, 'YYZ', DX, DY, DZ, R2, CHI)
!
!       FOR Fxzz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(17)+1, ITP(18), PTS, 'XZZ', DX, DY, DZ, R2, CHI)
!
!       FOR Fyzz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(18)+1, ITP(19), PTS, 'YZZ', DX, DY, DZ, R2, CHI)
!
!       FOR Fxyz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(19)+1, ITP(20), PTS, 'XYZ', DX, DY, DZ, R2, CHI)
!
   DO 603 J=1,Nprims
      temp=zero
      DO 604 I=1,PTs
         check=DMax1(Dabs(CHi(I,J)),Temp)
         If(Check.gt.temp)temp=check
604   Continue
      chimax(j)=temp
603 Continue
!
   DO 605 L = 1,NMO
      DO 607 I=1,PTS
         PSI(I,L) = ZERO
607   CONTINUE
!
      DO 625 J = 1,NPRIMS
         check=dabs(chimax(j)*coo(J,L))
         IF(check.gt.cutoff)THEN
            TEMP=COO(J,L)
            DO 626 I = 1,PTS
               PSI(I,L) = PSI(I,L) + TEMP*CHI(I,J)
626         CONTINUE
         ENDIF
625   CONTINUE
605 CONTINUE
!
   RETURN
END

!
! Helper subroutines for GAUS_PSI_ONLY orbital calculations
!

SUBROUTINE CALC_S_TYPE_PSI(J_START, J_END, PTS, DX, DY, DZ, R2, CHI)
!
! Calculate S-type orbital values (no angular part)
!
   USE molecular_data
   USE basis_data
   implicit none
   
   ! Arguments
   integer, intent(in) :: J_START, J_END, PTS
   real(8), intent(in) :: DX(2000,200), DY(2000,200), DZ(2000,200)
   real(8), intent(in) :: R2(2000,200)
   real(8), intent(inout) :: CHI(2000,2000)
   
   ! Local variables
   integer :: I, J, IS
   
   DO J = J_START, J_END
      IS = ICT(J)
      DO I = 1, PTS
         CHI(I,J) = DEXP(-EXX(J)*R2(I,IS))
      END DO
   END DO
   
   RETURN
END SUBROUTINE CALC_S_TYPE_PSI

SUBROUTINE CALC_P_TYPE_PSI(J_START, J_END, PTS, COORD_TYPE, DX, DY, DZ, R2, CHI)
!
! Calculate P-type orbital values (linear in coordinates)
!
   USE molecular_data
   USE basis_data
   implicit none
   
   ! Arguments
   integer, intent(in) :: J_START, J_END, PTS
   character(len=1), intent(in) :: COORD_TYPE
   real(8), intent(in) :: DX(2000,200), DY(2000,200), DZ(2000,200)
   real(8), intent(in) :: R2(2000,200)
   real(8), intent(inout) :: CHI(2000,2000)
   
   ! Local variables
   integer :: I, J, IS
   real(8) :: COORD_VAL
   
   DO J = J_START, J_END
      IS = ICT(J)
      DO I = 1, PTS
         SELECT CASE(COORD_TYPE)
            CASE('X')
               COORD_VAL = DX(I,IS)
            CASE('Y')
               COORD_VAL = DY(I,IS)
            CASE('Z')
               COORD_VAL = DZ(I,IS)
         END SELECT
         CHI(I,J) = COORD_VAL * DEXP(-EXX(J)*R2(I,IS))
      END DO
   END DO
   
   RETURN
END SUBROUTINE CALC_P_TYPE_PSI

SUBROUTINE CALC_D_TYPE_PSI(J_START, J_END, PTS, D_TYPE, DX, DY, DZ, R2, CHI)
!
! Calculate D-type orbital values (quadratic in coordinates)
!
   USE molecular_data
   USE basis_data
   implicit none
   
   ! Arguments
   integer, intent(in) :: J_START, J_END, PTS
   character(len=2), intent(in) :: D_TYPE
   real(8), intent(in) :: DX(2000,200), DY(2000,200), DZ(2000,200)
   real(8), intent(in) :: R2(2000,200)
   real(8), intent(inout) :: CHI(2000,2000)
   
   ! Local variables
   integer :: I, J, IS
   real(8) :: COORD_PROD
   
   DO J = J_START, J_END
      IS = ICT(J)
      DO I = 1, PTS
         SELECT CASE(D_TYPE)
            CASE('XX')
               COORD_PROD = DX(I,IS) * DX(I,IS)
            CASE('YY')
               COORD_PROD = DY(I,IS) * DY(I,IS)
            CASE('ZZ')
               COORD_PROD = DZ(I,IS) * DZ(I,IS)
            CASE('XY')
               COORD_PROD = DX(I,IS) * DY(I,IS)
            CASE('XZ')
               COORD_PROD = DX(I,IS) * DZ(I,IS)
            CASE('YZ')
               COORD_PROD = DY(I,IS) * DZ(I,IS)
         END SELECT
         CHI(I,J) = COORD_PROD * DEXP(-EXX(J)*R2(I,IS))
      END DO
   END DO
   
   RETURN
END SUBROUTINE CALC_D_TYPE_PSI

SUBROUTINE CALC_F_TYPE_PSI(J_START, J_END, PTS, F_TYPE, DX, DY, DZ, R2, CHI)
!
! Calculate F-type orbital values (cubic in coordinates)
!
   USE molecular_data
   USE basis_data
   implicit none
   
   ! Arguments
   integer, intent(in) :: J_START, J_END, PTS
   character(len=3), intent(in) :: F_TYPE
   real(8), intent(in) :: DX(2000,200), DY(2000,200), DZ(2000,200)
   real(8), intent(in) :: R2(2000,200)
   real(8), intent(inout) :: CHI(2000,2000)
   
   ! Local variables
   integer :: I, J, IS
   real(8) :: COORD_PROD
   
   DO J = J_START, J_END
      IS = ICT(J)
      DO I = 1, PTS
         SELECT CASE(F_TYPE)
            CASE('XXX')
               COORD_PROD = DX(I,IS)**3
            CASE('YYY')
               COORD_PROD = DY(I,IS)**3
            CASE('ZZZ')
               COORD_PROD = DZ(I,IS)**3
            CASE('XXY')
               COORD_PROD = DX(I,IS)**2 * DY(I,IS)
            CASE('XXZ')
               COORD_PROD = DX(I,IS)**2 * DZ(I,IS)
            CASE('YYZ')
               COORD_PROD = DY(I,IS)**2 * DZ(I,IS)
            CASE('XYY')
               COORD_PROD = DX(I,IS) * DY(I,IS)**2
            CASE('XZZ')
               COORD_PROD = DX(I,IS) * DZ(I,IS)**2
            CASE('YZZ')
               COORD_PROD = DY(I,IS) * DZ(I,IS)**2
            CASE('XYZ')
               COORD_PROD = DX(I,IS) * DY(I,IS) * DZ(I,IS)
         END SELECT
         CHI(I,J) = COORD_PROD * DEXP(-EXX(J)*R2(I,IS))
      END DO
   END DO
   
   RETURN
END SUBROUTINE CALC_F_TYPE_PSI
