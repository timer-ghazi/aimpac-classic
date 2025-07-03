! Generated from gridv_calculations.f90
! Contains: GRDD2R, GRDKEG, GRDRHO, GAUS

SUBROUTINE GRDD2R (A,CX,XY)
!
   USE grid_data
   USE molecular_data
   USE work_arrays
   implicit none
!
   ! Arguments
   real(8), intent(in) :: A(3,3), CX(3), XY(4)
   
   ! Local variables
   real(8) :: XYZ(2000,3), SUMR(2000)
   real(8), save :: Zero = 0.0d0, Two = 2.0d0
   integer :: IXSTP, IYSTP, IY, I, ICALC, INDX, J
   real(8) :: XCN, YCN, ZCN, XCNT, YCNT, ZCNT
!
!    CALCULATE NUMBER OF STEPS IN X AND Y FOR THIS GRID AS WELL AS
!    INCREMENT MARKERS FOR CENTERING PLOT
!
   IXSTP = IDINT(XY(1)/XY(2))
   IYSTP = IXSTP
   XCN = -XY(1)/Two
   YCN = XCN
   ZCN = Zero
   YCNT = XCN
   ZCNT=ZEro
!
   DO 1400 IY = 1,IYSTP
!
!    RESET X INCREMENT
!
      XCNT = XCN
!
      DO 1300 I = 1,IXSTP
!
         SUMR(I)=Zero
!     BACK TRANSFORM PLANE INTO ORIGINAL MOLECULAR SYSTEM
!
         XYZ(I,1)=A(1,1)*XCNT+A(1,2)*YCNT+A(1,3)*ZCNT+CX(1)
         XYZ(I,2)=A(2,1)*XCNT+A(2,2)*YCNT+A(2,3)*ZCNT+CX(2)
         XYZ(I,3)=A(3,1)*XCNT+A(3,2)*YCNT+A(3,3)*ZCNT+CX(3)
!
         XCNT=XCNT+XY(2)
1300  CONTINUE
!
      ICALC = 0
      CALL GAUS(XYZ,IXSTP,ICALC)
!
      INDX=(IY-1)*IXSTP
      DO 904 J=1,NMO
         DO 905 I = 1,IXSTP
            SUMR(I)=SUMR(I)+Two*PO(J)*(PSI(I,J)*D2(I,J)&
            &+(GX(I,J)**2+GY(I,J)**2+GZ(I,J)**2))
905      CONTINUE
904   CONTINUE
!
      DO 906 I=1,IXSTP
         GRD(INDX+I)=SUMR(I)
906   CONTINUE
!
!    INCREMENT Y VALUE
!
      YCNT=YCNT+XY(2)
!
1400 CONTINUE
!
   RETURN
END

SUBROUTINE GRDKEG (A,CX,XY)
!
   USE grid_data
   USE molecular_data
   USE work_arrays
   implicit none
!
   ! Arguments
   real(8), intent(in) :: A(3,3), CX(3), XY(4)
   
   ! Local variables
   real(8) :: XYZ(2000,3), SUMR(2000)
   real(8), save :: Zero = 0.0d0, Two = 2.0d0, Pt5 = 0.5d0
   integer :: IXSTP, IYSTP, IY, I, ICALC, INDX, J
   real(8) :: XCN, YCN, ZCN, XCNT, YCNT, ZCNT
!
!    CALCULATE NUMBER OF STEPS IN X AND Y FOR THIS GRID AS WELL AS
!    INCREMENT MARKERS FOR CENTERING PLOT
!
   IXSTP = IDINT(XY(1)/XY(2))
   IYSTP = IXSTP
   XCN = -XY(1)/Two
   YCN = XCN
   ZCN = Zero
   YCNT = XCN
   ZCnt=Zero
!
   DO 1400 IY = 1,IYSTP
!
!    RESET X INCREMENT
!
      XCNT = XCN
!
      DO 1300 I = 1,IXSTP
!
         SUMR(I)=Zero
!     BACK TRANSFORM PLANE INTO ORIGINAL MOLECULAR SYSTEM
!
         XYZ(I,1)=A(1,1)*XCNT+A(1,2)*YCNT+A(1,3)*ZCNT+CX(1)
         XYZ(I,2)=A(2,1)*XCNT+A(2,2)*YCNT+A(2,3)*ZCNT+CX(2)
         XYZ(I,3)=A(3,1)*XCNT+A(3,2)*YCNT+A(3,3)*ZCNT+CX(3)
!
         XCNT=XCNT+XY(2)
1300  CONTINUE
!
      CALL GAUS(XYZ,IXSTP,0)
!
      INDX=(IY-1)*IXSTP
      DO 904 J=1,NMO
         DO 905 I = 1,IXSTP
            SUMR(I)=SUMR(I)+Pt5*PO(J)*&
            &(GX(I,J)**2+GY(I,J)**2+GZ(I,J)**2)
905      CONTINUE
904   CONTINUE
!
      DO 906 I=1,IXSTP
         GRD(INDX+I)=SUMR(I)
906   CONTINUE
!
!    INCREMENT Y VALUE
!
      YCNT=YCNT+XY(2)
!
1400 CONTINUE
!
   RETURN
END

SUBROUTINE GRDRHO (A,CX,XY)
!
   USE grid_data
   USE molecular_data
   USE work_arrays
   implicit none
!
   ! Arguments
   real(8), intent(in) :: A(3,3), CX(3), XY(4)
   
   ! Local variables
   real(8) :: XYZ(2000,3), SUMR(2000)
   real(8), save :: Zero = 0.0d0, Two = 2.0d0
   integer :: IXSTP, IYSTP, IY, I, ICALC, INDX, J
   real(8) :: XCN, YCN, ZCN, XCNT, YCNT, ZCNT
!
!    CALCULATE NUMBER OF STEPS IN X AND Y FOR THIS GRID AS WELL AS
!    INCREMENT MARKERS FOR CENTERING PLOT
!
   IXSTP = IDINT(XY(1)/XY(2))
   IYSTP = IXSTP
   XCN = -XY(1)/Two
   YCN = XCN
   ZCN = Zero
   YCNT = XCN
   ZCnt=Zero
!
   DO 1400 IY = 1,IYSTP
!
!    RESET X INCREMENT
!
      XCNT = XCN
!
      DO 1300 I = 1,IXSTP
!
         SUMR(I)=Zero
!     BACK TRANSFORM PLANE INTO ORIGINAL MOLECULAR SYSTEM
!
         XYZ(I,1)=A(1,1)*XCNT+A(1,2)*YCNT+A(1,3)*ZCNT+CX(1)
         XYZ(I,2)=A(2,1)*XCNT+A(2,2)*YCNT+A(2,3)*ZCNT+CX(2)
         XYZ(I,3)=A(3,1)*XCNT+A(3,2)*YCNT+A(3,3)*ZCNT+CX(3)
!
         XCNT=XCNT+XY(2)
1300  CONTINUE
!
      ICALC = 1
      CALL GAUS(XYZ,IXSTP,ICALC)
!
      INDX=(IY-1)*IXSTP
      DO 904 J=1,NMO
         DO 905 I = 1,IXSTP
            SUMR(I)=SUMR(I)+PO(J)*PSI(I,J)*PSI(I,J)
905      CONTINUE
904   CONTINUE
!
      DO 906 I=1,IXSTP
         GRD(INDX+I)=SUMR(I)
906   CONTINUE
!
!    INCREMENT Y VALUE
!
      YCNT=YCNT+XY(2)
!
1400 CONTINUE
!
   RETURN
END

SUBROUTINE GAUS(XYZ,PTS,MM)
!
   USE molecular_data
   USE basis_data
   USE work_arrays
   implicit none
!
   ! Arguments
   integer, intent(in) :: PTS, MM
   real(8), intent(in) :: XYZ(2000,3)
!
   ! Dispatcher: call appropriate subroutine based on MM
   IF (MM .EQ. 1) THEN
      CALL GAUS_PSI_ONLY(XYZ,PTS,MM)
   ELSE
      CALL GAUS_FULL(XYZ,PTS,MM)
   ENDIF
!
   RETURN
END
