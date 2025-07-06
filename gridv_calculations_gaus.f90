! GAUS Calculation Routines - DOCUMENTATION ADDED FOR FUTURE MODERNIZATION
! 
! Generated from gridv_calculations.f90
! Contains: GAUS_FULL, GAUS_PSI_ONLY
!
! NOTE: Converted to dynamic arrays with explicit dimension passing
! Uses local_arrays module for scalable memory allocation

SUBROUTINE GAUS_FULL(XYZ,PTS,MM)
!
   USE molecular_data
   USE basis_data
   USE work_arrays
   USE local_arrays
   implicit none
!
   ! Arguments
   integer, intent(in) :: PTS, MM
   real(8), intent(in) :: XYZ(2000,3)  ! Keep input array as is
   
   ! Local variables - remove all array declarations, keep only scalars
   real(8), save :: ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0, FIVE = 5.D0
   real(8), save :: SEVEN = 7.D0, THREE = 3.D0, SIX = 6.D0, NINE = 9.D0, CUTOFF = 1.0d-10
   integer :: I, J, K, L, IS
   real(8) :: A, B, BXX, BXY, BXZ, BYY, BYZ, BZZ, CHECK, TEMP
!
   ! Allocate local arrays for this calculation
   CALL allocate_local_arrays(PTS, NCENT, NPRIMS)
   DO 110 J = 1,NCENT
      DO 112 I=1,PTS
         DX_LOCAL(I,J) = XYZ(I,1) - XC(J)
         DY_LOCAL(I,J) = XYZ(I,2) - YC(J)
         DZ_LOCAL(I,J) = XYZ(I,3) - ZC(J)
         R2_LOCAL(I,J)= DX_LOCAL(I,J)*DX_LOCAL(I,J)+DY_LOCAL(I,J)*DY_LOCAL(I,J)+DZ_LOCAL(I,J)*DZ_LOCAL(I,J)
112   CONTINUE
110 CONTINUE
!
!       FOR S-TYPE
!
   DO 120 J = 1,ITP(1)
      IS=ict(j)
      DO 122 I=1,PTS
         A=SUM(J)*DEXP(-EXX(J)*R2_LOCAL(I,is))
         CHI_LOCAL(I,J)=A*DIV(J)
         CHIX_LOCAL(I,J)=DX_LOCAL(I,is)*A
         CHIY_LOCAL(I,J)=DY_LOCAL(I,is)*A
         CHIZ_LOCAL(I,J)=DZ_LOCAL(I,is)*A
         CHID2_LOCAL(I,J)=(THREE+SUM(J)*R2_LOCAL(I,is))*A
122   CONTINUE
120 CONTINUE
!
!       FOR Px-TYPE
!
   DO 140 J=ITP(1)+1,ITP(2)
      IS=ict(j)
      DO 142 I=1,PTS
!
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DX_LOCAL(I,is)*A*SUM(J)
!
         CHI_LOCAL(I,J)=A*DX_LOCAL(I,is)
         CHIX_LOCAL(I,J)=A+DX_LOCAL(I,is)*B
         CHIY_LOCAL(I,J)=DY_LOCAL(I,is)*B
         CHIZ_LOCAL(I,J)=DZ_LOCAL(I,is)*B
         CHID2_LOCAL(I,J)=(FIVE+SUM(J)*R2_LOCAL(I,is))*B
142   CONTINUE
140 CONTINUE
!
!       FOR Py-TYPE
!
   DO 160 J=ITP(2)+1,ITP(3)
      IS=ict(j)
      DO 162 I=1,PTS
!
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DY_LOCAL(I,is)*A*SUM(J)
!
         CHI_LOCAL(I,J)=A*DY_LOCAL(I,is)
         CHIX_LOCAL(I,J)=DX_LOCAL(I,is)*B
         CHIY_LOCAL(I,J)=A+DY_LOCAL(I,is)*B
         CHIZ_LOCAL(I,J)=DZ_LOCAL(I,is)*B
         CHID2_LOCAL(I,J)=(FIVE+SUM(J)*R2_LOCAL(I,is))*B
162   CONTINUE
160 CONTINUE
!
!       FOR Pz-TYPE
!
   DO 180 J=ITP(3)+1,ITP(4)
      IS=ict(j)
      DO 182 I=1,PTS
!
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DZ_LOCAL(I,is)*A*SUM(J)
!
         CHI_LOCAL(I,J)=A*DZ_LOCAL(I,is)
         CHIX_LOCAL(I,J)=DX_LOCAL(I,is)*B
         CHIY_LOCAL(I,J)=DY_LOCAL(I,is)*B
         CHIZ_LOCAL(I,J)=A+DZ_LOCAL(I,is)*B
         CHID2_LOCAL(I,J)=(FIVE+SUM(J)*R2_LOCAL(I,is))*B
182   CONTINUE
180 CONTINUE
!
!       FOR Dxx-TYPE
!
   DO 220 J=ITP(4)+1,ITP(5)
      IS=ict(j)
      DO 222 I=1,PTS
!
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DX_LOCAL(I,is)*DX_LOCAL(I,is)*A*SUM(J)
!
         CHI_LOCAL(I,J)=B*DIV(J)
         CHIX_LOCAL(I,J)=(TWO*A+B)*DX_LOCAL(I,is)
         CHIY_LOCAL(I,J)=DY_LOCAL(I,is)*B
         CHIZ_LOCAL(I,J)=DZ_LOCAL(I,is)*B
         CHID2_LOCAL(I,J)=TWO*A+(SEVEN+SUM(J)*R2_LOCAL(I,is))*B
222   CONTINUE
220 CONTINUE
!
!       FOR Dyy-TYPE
!
   DO 240 J=ITP(5)+1,ITP(6)
      IS=ict(j)
      DO 242 I=1,PTS
!
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DY_LOCAL(I,is)*DY_LOCAL(I,is)*A*SUM(J)
!
         CHI_LOCAL(I,J)=B*DIV(J)
         CHIX_LOCAL(I,J)=DX_LOCAL(I,is)*B
         CHIY_LOCAL(I,J)=(TWO*A+B)*DY_LOCAL(I,is)
         CHIZ_LOCAL(I,J)=DZ_LOCAL(I,is)*B
         CHID2_LOCAL(I,J)=TWO*A+(SEVEN+SUM(J)*R2_LOCAL(I,is))*B
242   CONTINUE
240 CONTINUE
!
!       FOR Dzz-TYPE
!
   DO 260 J=ITP(6)+1,ITP(7)
      IS=ict(j)
      DO 262 I=1,PTS
!
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DZ_LOCAL(I,is)*DZ_LOCAL(I,is)*A*SUM(J)
!
         CHI_LOCAL(I,J)=B*DIV(J)
         CHIX_LOCAL(I,J)=DX_LOCAL(I,is)*B
         CHIY_LOCAL(I,J)=DY_LOCAL(I,is)*B
         CHIZ_LOCAL(I,J)=(TWO*A+B)*DZ_LOCAL(I,is)
         CHID2_LOCAL(I,J)=TWO*A+(SEVEN+SUM(J)*R2_LOCAL(I,is))*B
262   CONTINUE
260 CONTINUE
!
!       FOR Dxy-TYPE
!
   DO 280 J=ITP(7)+1,ITP(8)
      IS=ict(j)
      DO 282 I=1,PTS
!
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DX_LOCAL(I,is)*DY_LOCAL(I,is)*A*SUM(J)
!
         CHI_LOCAL(I,J)=B*DIV(J)
         CHIX_LOCAL(I,J)=DX_LOCAL(I,is)*B+DY_LOCAL(I,is)*A
         CHIY_LOCAL(I,J)=DY_LOCAL(I,is)*B+DX_LOCAL(I,is)*A
         CHIZ_LOCAL(I,J)=DZ_LOCAL(I,is)*B
         CHID2_LOCAL(I,J)=(SEVEN+SUM(J)*R2_LOCAL(I,is))*B
282   CONTINUE
280 CONTINUE
!
!       FOR Dxz-TYPE
!
   DO 320 J=ITP(8)+1,ITP(9)
      IS=ict(j)
      DO 322 I=1,PTS
!
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DX_LOCAL(I,is)*DZ_LOCAL(I,is)*A*SUM(J)
!
         CHI_LOCAL(I,J)=B*DIV(J)
         CHIX_LOCAL(I,J)=DX_LOCAL(I,is)*B+DZ_LOCAL(I,is)*A
         CHIY_LOCAL(I,J)=DY_LOCAL(I,is)*B
         CHIZ_LOCAL(I,J)=DZ_LOCAL(I,is)*B+DX_LOCAL(I,is)*A
         CHID2_LOCAL(I,J)=(SEVEN+SUM(J)*R2_LOCAL(I,is))*B
322   CONTINUE
320 CONTINUE
!
!       FOR Dyz-TYPE
!
   DO 340 J=ITP(9)+1,ITP(10)
      IS=ict(j)
      DO 342 I=1,PTS
!
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DY_LOCAL(I,is)*DZ_LOCAL(I,is)*A*SUM(J)
!
         CHI_LOCAL(I,J)=B*DIV(J)
         CHIX_LOCAL(I,J)=DX_LOCAL(I,is)*B
         CHIY_LOCAL(I,J)=DY_LOCAL(I,is)*B+DZ_LOCAL(I,is)*A
         CHIZ_LOCAL(I,J)=DZ_LOCAL(I,is)*B+DY_LOCAL(I,is)*A
         CHID2_LOCAL(I,J)=(SEVEN+SUM(J)*R2_LOCAL(I,is))*B
342   CONTINUE
340 CONTINUE
!
!       FOR Fxxx-TYPE
!
   DO 501 J=ITP(10)+1,ITP(11)
      IS=ict(j)
      DO 502 I=1,PTS

         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DX_LOCAL(I,is)*DX_LOCAL(I,is)*A
!
         CHI_LOCAL(I,J)=B*DX_LOCAL(I,is)
         CHIX_LOCAL(I,J)=(THREE + SUM(J)*DX_LOCAL(I,is)*DX_LOCAL(I,is))*B
         CHIY_LOCAL(I,J)=SUM(J)*DY_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHIZ_LOCAL(I,J)=SUM(J)*DZ_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHID2_LOCAL(I,J)=SIX*A*DX_LOCAL(I,is)+&
         &(NINE+SUM(J)*R2_LOCAL(I,is))*CHI_LOCAL(I,J)*SUM(J)
502   CONTINUE
501 CONTINUE
!
!       FOR Fyyy-TYPE
!
   DO 511 J=ITP(11)+1,ITP(12)
      IS=ict(j)
      DO 512 I=1,PTS
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DY_LOCAL(I,is)*DY_LOCAL(I,is)*A
!
         CHI_LOCAL(I,J)=B*DY_LOCAL(I,is)
         CHIX_LOCAL(I,J)=SUM(J)*DX_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHIY_LOCAL(I,J)=(THREE + SUM(J)*DY_LOCAL(I,is)*DY_LOCAL(I,is))*B
         CHIZ_LOCAL(I,J)=SUM(J)*DZ_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHID2_LOCAL(I,J)=SIX*A*DY_LOCAL(I,is)+&
         &(NINE+SUM(J)*R2_LOCAL(I,is))*CHI_LOCAL(I,J)*SUM(J)
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
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         B=DZ_LOCAL(I,is)*DZ_LOCAL(I,is)*A
!
         CHI_LOCAL(I,J)=B*DZ_LOCAL(I,is)
         CHIX_LOCAL(I,J)=SUM(J)*DX_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHIY_LOCAL(I,J)=SUM(J)*DY_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHIZ_LOCAL(I,J)=(THREE + SUM(J)*DZ_LOCAL(I,is)*DZ_LOCAL(I,is))*B
         CHID2_LOCAL(I,J)=SIX*A*DZ_LOCAL(I,is)+&
         &(NINE+SUM(J)*R2_LOCAL(I,is))*CHI_LOCAL(I,J)*SUM(J)
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
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         BXY=DX_LOCAL(I,is)*DY_LOCAL(I,is)*A
         BXX=DX_LOCAL(I,is)*DX_LOCAL(I,is)*A
!
         CHI_LOCAL(I,J)=BXY*DX_LOCAL(I,is)
         CHIX_LOCAL(I,J)=(TWO + SUM(J)*DX_LOCAL(I,is)*DX_LOCAL(I,is))*BXY
         CHIY_LOCAL(I,J)=(ONE + SUM(J)*DY_LOCAL(I,is)*DY_LOCAL(I,is))*BXX
         CHIZ_LOCAL(I,J)=SUM(J)*DZ_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHID2_LOCAL(I,J)=TWO*A*DY_LOCAL(I,is)+&
         &(NINE+SUM(J)*R2_LOCAL(I,is))*CHI_LOCAL(I,J)*SUM(J)
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
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         BXZ=DX_LOCAL(I,is)*DZ_LOCAL(I,is)*A
         BXX=DX_LOCAL(I,is)*DX_LOCAL(I,is)*A
!
         CHI_LOCAL(I,J)=BXZ*DX_LOCAL(I,is)
         CHIX_LOCAL(I,J)=(TWO + SUM(J)*DX_LOCAL(I,is)*DX_LOCAL(I,is))*BXZ
         CHIY_LOCAL(I,J)=SUM(J)*DY_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHIZ_LOCAL(I,J)=(ONE + SUM(J)*DZ_LOCAL(I,is)*DZ_LOCAL(I,is))*BXX
         CHID2_LOCAL(I,J)=TWO*A*DZ_LOCAL(I,is)+&
         &(NINE+SUM(J)*R2_LOCAL(I,is))*CHI_LOCAL(I,J)*SUM(J)
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
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         BYZ=DZ_LOCAL(I,is)*DY_LOCAL(I,is)*A
         BYY=DY_LOCAL(I,is)*DY_LOCAL(I,is)*A
!
         CHI_LOCAL(I,J)=BYZ*DY_LOCAL(I,is)
         CHIX_LOCAL(I,J)=SUM(J)*DX_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHIY_LOCAL(I,J)=(TWO + SUM(J)*DY_LOCAL(I,is)*DY_LOCAL(I,is))*BYZ
         CHIZ_LOCAL(I,J)=(ONE + SUM(J)*DZ_LOCAL(I,is)*DZ_LOCAL(I,is))*BYY
         CHID2_LOCAL(I,J)=TWO*A*DZ_LOCAL(I,is)+&
         &(NINE+SUM(J)*R2_LOCAL(I,is))*CHI_LOCAL(I,J)*SUM(J)
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
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         BXY=DX_LOCAL(I,is)*DY_LOCAL(I,is)*A
         BYY=DY_LOCAL(I,is)*DY_LOCAL(I,is)*A
!
         CHI_LOCAL(I,J)=BXY*DY_LOCAL(I,is)
         CHIX_LOCAL(I,J)=(ONE + SUM(J)*DX_LOCAL(I,is)*DX_LOCAL(I,is))*BYY
         CHIY_LOCAL(I,J)=(TWO + SUM(J)*DY_LOCAL(I,is)*DY_LOCAL(I,is))*BXY
         CHIZ_LOCAL(I,J)=SUM(J)*DZ_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHID2_LOCAL(I,J)=TWO*A*DX_LOCAL(I,is)+&
         &(NINE+SUM(J)*R2_LOCAL(I,is))*CHI_LOCAL(I,J)*SUM(J)
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
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         BXZ=DZ_LOCAL(I,is)*DX_LOCAL(I,is)*A
         BZZ=DZ_LOCAL(I,is)*DZ_LOCAL(I,is)*A
!
         CHI_LOCAL(I,J)=BXZ*DZ_LOCAL(I,is)
         CHIX_LOCAL(I,J)=(ONE + SUM(J)*DX_LOCAL(I,is)*DX_LOCAL(I,is))*BZZ
         CHIY_LOCAL(I,J)=SUM(J)*DY_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHIZ_LOCAL(I,J)=(TWO + SUM(J)*DZ_LOCAL(I,is)*DZ_LOCAL(I,is))*BXZ
         CHID2_LOCAL(I,J)=TWO*A*DX_LOCAL(I,is)+&
         &(NINE+SUM(J)*R2_LOCAL(I,is))*CHI_LOCAL(I,J)*SUM(J)
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
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         BYZ=DZ_LOCAL(I,is)*DY_LOCAL(I,is)*A
         BZZ=DZ_LOCAL(I,is)*DZ_LOCAL(I,is)*A
!
         CHI_LOCAL(I,J)=BYZ*DZ_LOCAL(I,is)
         CHIX_LOCAL(I,J)=SUM(J)*DX_LOCAL(I,is)*CHI_LOCAL(I,J)
         CHIY_LOCAL(I,J)=(ONE + SUM(J)*DY_LOCAL(I,is)*DY_LOCAL(I,is))*BZZ
         CHIZ_LOCAL(I,J)=(TWO + SUM(J)*DZ_LOCAL(I,is)*DZ_LOCAL(I,is))*BYZ
         CHID2_LOCAL(I,J)=TWO*A*DY_LOCAL(I,is)+&
         &(NINE+SUM(J)*R2_LOCAL(I,is))*CHI_LOCAL(I,J)*SUM(J)
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
         A=DEXP(-EXX(J)*R2_LOCAL(I,is))
         BXY=DX_LOCAL(I,is)*DY_LOCAL(I,is)*A
         BYZ=DZ_LOCAL(I,is)*DY_LOCAL(I,is)*A
         BXZ=DX_LOCAL(I,is)*DZ_LOCAL(I,is)*A
!
         CHI_LOCAL(I,J)=DX_LOCAL(I,is)*BYZ
         CHIX_LOCAL(I,J)=(ONE + SUM(J)*DX_LOCAL(I,is)*DX_LOCAL(I,is))*BYZ
         CHIY_LOCAL(I,J)=(ONE + SUM(J)*DY_LOCAL(I,is)*DY_LOCAL(I,is))*BXZ
         CHIZ_LOCAL(I,J)=(ONE + SUM(J)*DZ_LOCAL(I,is)*DZ_LOCAL(I,is))*BXY
         CHID2_LOCAL(I,J)=(NINE+SUM(J)*R2_LOCAL(I,is))*CHI_LOCAL(I,J)*SUM(J)
!
592   CONTINUE
591 CONTINUE
!
   DO 103 J=1,Nprims
      temp=zero
      DO 104 I=1,PTs
         check=DMax1(Dabs(CHI_LOCAL(I,J)),Dabs(CHIX_LOCAL(I,j)),Dabs(CHIY_LOCAL(I,J)),&
         &Dabs(CHIZ_LOCAL(I,j)),Dabs(CHID2_LOCAL(I,J)),Temp)
         If(Check.gt.temp)temp=check
104   Continue
      chimax_local(j)=temp
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
         check=dabs(chimax_local(j)*coo(J,L))
         IF(check.gt.cutoff)THEN
            TEMP=COO(J,L)
            DO 126 I = 1,PTS
               PSI(I,L) = PSI(I,L) + TEMP*CHI_LOCAL(I,J)
               GX(I,L) = GX(I,L) + TEMP*CHIX_LOCAL(I,J)
               GY(I,L) = GY(I,L) + TEMP*CHIY_LOCAL(I,J)
               GZ(I,L) = GZ(I,L) + TEMP*CHIZ_LOCAL(I,J)
               D2(I,L) = D2(I,L) + TEMP*CHID2_LOCAL(I,J)
126         CONTINUE
         ENDIF
125   CONTINUE
105 CONTINUE
!
   ! Cleanup local arrays
   CALL cleanup_local_arrays()
   
   RETURN
END

SUBROUTINE GAUS_PSI_ONLY(XYZ,PTS,MM)
!
   USE molecular_data
   USE basis_data
   USE work_arrays
   USE local_arrays
   implicit none
!
   ! Arguments
   integer, intent(in) :: PTS, MM
   real(8), intent(in) :: XYZ(2000,3)  ! Keep input array as is
   
   ! Local variables - remove all array declarations, keep only scalars
   real(8), save :: ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0, FIVE = 5.D0
   real(8), save :: SEVEN = 7.D0, THREE = 3.D0, SIX = 6.D0, NINE = 9.D0, CUTOFF = 1.0d-10
   integer :: I, J, K, L, IS
   real(8) :: A, B, BXX, BXY, BXZ, BYY, BYZ, BZZ, CHECK, TEMP
!
   ! Allocate local arrays for this calculation  
   CALL allocate_local_arrays(PTS, NCENT, NPRIMS)  
   DO 410 J = 1,NCENT
      DO 412 I=1,PTS
         DX_LOCAL(I,J) = XYZ(I,1) - XC(J)
         DY_LOCAL(I,J) = XYZ(I,2) - YC(J)
         DZ_LOCAL(I,J) = XYZ(I,3) - ZC(J)
         R2_LOCAL(I,J)= DX_LOCAL(I,J)*DX_LOCAL(I,J)+DY_LOCAL(I,J)*DY_LOCAL(I,J)+DZ_LOCAL(I,J)*DZ_LOCAL(I,J)
412   CONTINUE
410 CONTINUE
!
!       FOR S-TYPE
!
   CALL CALC_S_TYPE_PSI(1, ITP(1), PTS, PTS, NCENT, NPRIMS, &
                         DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Px-TYPE
!
   CALL CALC_P_TYPE_PSI(ITP(1)+1, ITP(2), PTS, PTS, NCENT, NPRIMS, &
                         'X', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Py-TYPE
!
   CALL CALC_P_TYPE_PSI(ITP(2)+1, ITP(3), PTS, PTS, NCENT, NPRIMS, &
                         'Y', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Pz-TYPE
!
   CALL CALC_P_TYPE_PSI(ITP(3)+1, ITP(4), PTS, PTS, NCENT, NPRIMS, &
                         'Z', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Dxx-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(4)+1, ITP(5), PTS, PTS, NCENT, NPRIMS, &
                         'XX', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Dyy-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(5)+1, ITP(6), PTS, PTS, NCENT, NPRIMS, &
                         'YY', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Dzz-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(6)+1, ITP(7), PTS, PTS, NCENT, NPRIMS, &
                         'ZZ', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Dxy-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(7)+1, ITP(8), PTS, PTS, NCENT, NPRIMS, &
                         'XY', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Dxz-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(8)+1, ITP(9), PTS, PTS, NCENT, NPRIMS, &
                         'XZ', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Dyz-TYPE
!
   CALL CALC_D_TYPE_PSI(ITP(9)+1, ITP(10), PTS, PTS, NCENT, NPRIMS, &
                         'YZ', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Fxxx-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(10)+1, ITP(11), PTS, PTS, NCENT, NPRIMS, &
                         'XXX', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Fyyy-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(11)+1, ITP(12), PTS, PTS, NCENT, NPRIMS, &
                         'YYY', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Fzzz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(12)+1, ITP(13), PTS, PTS, NCENT, NPRIMS, &
                         'ZZZ', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Fxxy-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(13)+1, ITP(14), PTS, PTS, NCENT, NPRIMS, &
                         'XXY', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Fxxz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(14)+1, ITP(15), PTS, PTS, NCENT, NPRIMS, &
                         'XXZ', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Fxyy-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(16)+1, ITP(17), PTS, PTS, NCENT, NPRIMS, &
                         'XYY', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Fyyz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(15)+1, ITP(16), PTS, PTS, NCENT, NPRIMS, &
                         'YYZ', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Fxzz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(17)+1, ITP(18), PTS, PTS, NCENT, NPRIMS, &
                         'XZZ', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Fyzz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(18)+1, ITP(19), PTS, PTS, NCENT, NPRIMS, &
                         'YZZ', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
!       FOR Fxyz-TYPE
!
   CALL CALC_F_TYPE_PSI(ITP(19)+1, ITP(20), PTS, PTS, NCENT, NPRIMS, &
                         'XYZ', DX_LOCAL, DY_LOCAL, DZ_LOCAL, R2_LOCAL, CHI_LOCAL)
!
   DO 603 J=1,Nprims
      temp=zero
      DO 604 I=1,PTs
         check=DMax1(Dabs(CHI_LOCAL(I,J)),Temp)
         If(Check.gt.temp)temp=check
604   Continue
      chimax_local(j)=temp
603 Continue
!
   DO 605 L = 1,NMO
      DO 607 I=1,PTS
         PSI(I,L) = ZERO
607   CONTINUE
!
      DO 625 J = 1,NPRIMS
         check=dabs(chimax_local(j)*coo(J,L))
         IF(check.gt.cutoff)THEN
            TEMP=COO(J,L)
            DO 626 I = 1,PTS
               PSI(I,L) = PSI(I,L) + TEMP*CHI_LOCAL(I,J)
626         CONTINUE
         ENDIF
625   CONTINUE
605 CONTINUE
!
   ! Cleanup local arrays
   CALL cleanup_local_arrays()
   
   RETURN
END

!
! Helper subroutines for GAUS_PSI_ONLY orbital calculations
!

SUBROUTINE CALC_S_TYPE_PSI(J_START, J_END, PTS, MAXPTS, MAXCENT, MAXPRIMS, &
                           DX, DY, DZ, R2, CHI)
!
! Calculate S-type orbital values (no angular part)
!
   USE molecular_data
   USE basis_data
   implicit none
   
   ! Arguments
   integer, intent(in) :: J_START, J_END, PTS, MAXPTS, MAXCENT, MAXPRIMS
   real(8), intent(in) :: DX(MAXPTS,MAXCENT), DY(MAXPTS,MAXCENT), DZ(MAXPTS,MAXCENT)
   real(8), intent(in) :: R2(MAXPTS,MAXCENT)
   real(8), intent(inout) :: CHI(MAXPTS,MAXPRIMS)
   
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

SUBROUTINE CALC_P_TYPE_PSI(J_START, J_END, PTS, MAXPTS, MAXCENT, MAXPRIMS, &
                           COORD_TYPE, DX, DY, DZ, R2, CHI)
!
! Calculate P-type orbital values (linear in coordinates)
!
   USE molecular_data
   USE basis_data
   implicit none
   
   ! Arguments
   integer, intent(in) :: J_START, J_END, PTS, MAXPTS, MAXCENT, MAXPRIMS
   character(len=1), intent(in) :: COORD_TYPE
   real(8), intent(in) :: DX(MAXPTS,MAXCENT), DY(MAXPTS,MAXCENT), DZ(MAXPTS,MAXCENT)
   real(8), intent(in) :: R2(MAXPTS,MAXCENT)
   real(8), intent(inout) :: CHI(MAXPTS,MAXPRIMS)
   
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

SUBROUTINE CALC_D_TYPE_PSI(J_START, J_END, PTS, MAXPTS, MAXCENT, MAXPRIMS, &
                           D_TYPE, DX, DY, DZ, R2, CHI)
!
! Calculate D-type orbital values (quadratic in coordinates)
!
   USE molecular_data
   USE basis_data
   implicit none
   
   ! Arguments
   integer, intent(in) :: J_START, J_END, PTS, MAXPTS, MAXCENT, MAXPRIMS
   character(len=2), intent(in) :: D_TYPE
   real(8), intent(in) :: DX(MAXPTS,MAXCENT), DY(MAXPTS,MAXCENT), DZ(MAXPTS,MAXCENT)
   real(8), intent(in) :: R2(MAXPTS,MAXCENT)
   real(8), intent(inout) :: CHI(MAXPTS,MAXPRIMS)
   
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

SUBROUTINE CALC_F_TYPE_PSI(J_START, J_END, PTS, MAXPTS, MAXCENT, MAXPRIMS, &
                           F_TYPE, DX, DY, DZ, R2, CHI)
!
! Calculate F-type orbital values (cubic in coordinates)
!
   USE molecular_data
   USE basis_data
   implicit none
   
   ! Arguments
   integer, intent(in) :: J_START, J_END, PTS, MAXPTS, MAXCENT, MAXPRIMS
   character(len=3), intent(in) :: F_TYPE
   real(8), intent(in) :: DX(MAXPTS,MAXCENT), DY(MAXPTS,MAXCENT), DZ(MAXPTS,MAXCENT)
   real(8), intent(in) :: R2(MAXPTS,MAXCENT)
   real(8), intent(inout) :: CHI(MAXPTS,MAXPRIMS)
   
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

! End of GAUS routines - Successfully converted to dynamic arrays
! Current status: Dynamic allocation with explicit dimension passing
! Implementation complete: Scalable and efficient memory usage
