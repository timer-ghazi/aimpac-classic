! Generated from gridv.f
! Contains: GAUS, GRDRHO, GRDD2R, GRDKEG


      SUBROUTINE GAUS(XYZ,PTS,MM)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION NINE
      INTEGER PTS
      PARAMETER(MCENT=200,MMO=500,MPRIMS=2000,MPTS=2000,NTYPE=20)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),NMO
      COMMON /PRIMS/ DIV(MPRIMS),COO(MPRIMS,MMO),EXX(MPRIMS),
     +ICT(MPRIMS),SUM(MPRIMS),ITP(NTYPE),NPRIMS
      COMMON /ZZZZ/ PSI(MPTS,MMO),GX(MPTS,MMO),GY(MPTS,MMO),
     +GZ(MPTS,MMO),D2(MPTS,MMO)
      DIMENSION XYZ(MPTS,3),DX(MPTS,MCENT),DY(MPTS,MCENT),
     +DZ(MPTS,MCENT),R2(MPTS,MCENT),CHI(MPTS,MPRIMS),
     +CHIX(MPTS,MPRIMS),CHIY(MPTS,MPRIMS),CHIZ(MPTS,MPRIMS),
     +CHID2(MPTS,MPRIMS),CHIMAX(MPRIMS)
      Save zero,one,two,four,five,seven,three,six,nine,cutoff
      DATA ZERO /0.D0/,ONE/1.D0/,TWO/2.D0/,FOUR/4.D0/,FIVE/5.D0/,
     +SEVEN/7.D0/,THREE/3.D0/,SIX/6.D0/,NINE/9.D0/,CUTOFF/1.0d-10/
C
      IF (MM .EQ. 1) GOTO 133
C
      DO 110 J = 1,NCENT
       DO 112 I=1,PTS 
        DX(I,J) = XYZ(I,1) - XC(J)
        DY(I,J) = XYZ(I,2) - YC(J)
        DZ(I,J) = XYZ(I,3) - ZC(J)
        R2(I,J)= DX(I,J)*DX(I,J)+DY(I,J)*DY(I,J)+DZ(I,J)*DZ(I,J)
112   CONTINUE
110   CONTINUE
C
C       FOR S-TYPE
C
        DO 120 J = 1,ITP(1) 
        IS=ict(j)
        DO 122 I=1,PTS 
        A=SUM(J)*DEXP(-EXX(J)*R2(I,is))
        CHI(I,J)=A*DIV(J)
        CHIX(I,J)=DX(I,is)*A
        CHIY(I,J)=DY(I,is)*A
        CHIZ(I,J)=DZ(I,is)*A
        CHID2(I,J)=(THREE+SUM(J)*R2(I,is))*A
122     CONTINUE
120     CONTINUE
C
C       FOR Px-TYPE
C
        DO 140 J=ITP(1)+1,ITP(2)
        IS=ict(j)
        DO 142 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DX(I,is)*A*SUM(J)
C
        CHI(I,J)=A*DX(I,is)
        CHIX(I,J)=A+DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=(FIVE+SUM(J)*R2(I,is))*B
142     CONTINUE
140     CONTINUE
C
C       FOR Py-TYPE
C
        DO 160 J=ITP(2)+1,ITP(3)
        IS=ict(j)
        DO 162 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DY(I,is)*A*SUM(J)
C
        CHI(I,J)=A*DY(I,is)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=A+DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=(FIVE+SUM(J)*R2(I,is))*B
162     CONTINUE
160     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 180 J=ITP(3)+1,ITP(4)
        IS=ict(j)
        DO 182 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DZ(I,is)*A*SUM(J)
C
        CHI(I,J)=A*DZ(I,is)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=A+DZ(I,is)*B
        CHID2(I,J)=(FIVE+SUM(J)*R2(I,is))*B
182     CONTINUE
180     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 220 J=ITP(4)+1,ITP(5)
        IS=ict(j)
        DO 222 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DX(I,is)*DX(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=(TWO*A+B)*DX(I,is)
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,is))*B
222     CONTINUE
220     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 240 J=ITP(5)+1,ITP(6)
        IS=ict(j)
        DO 242 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DY(I,is)*DY(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=(TWO*A+B)*DY(I,is)
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,is))*B
242     CONTINUE
240     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 260 J=ITP(6)+1,ITP(7)
        IS=ict(j)
        DO 262 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DZ(I,is)*DZ(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=(TWO*A+B)*DZ(I,is)
        CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,is))*B
262     CONTINUE
260     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 280 J=ITP(7)+1,ITP(8)
        IS=ict(j)
        DO 282 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DX(I,is)*DY(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,is)*B+DY(I,is)*A
        CHIY(I,J)=DY(I,is)*B+DX(I,is)*A
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=(SEVEN+SUM(J)*R2(I,is))*B
282     CONTINUE
280     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 320 J=ITP(8)+1,ITP(9)
        IS=ict(j)
        DO 322 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DX(I,is)*DZ(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,is)*B+DZ(I,is)*A
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B+DX(I,is)*A
        CHID2(I,J)=(SEVEN+SUM(J)*R2(I,is))*B
322     CONTINUE
320     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 340 J=ITP(9)+1,ITP(10)
        IS=ict(j)
        DO 342 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DY(I,is)*DZ(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B+DZ(I,is)*A
        CHIZ(I,J)=DZ(I,is)*B+DY(I,is)*A
        CHID2(I,J)=(SEVEN+SUM(J)*R2(I,is))*B
342     CONTINUE
340     CONTINUE
C
C       FOR Fxxx-TYPE
C
        DO 501 J=ITP(10)+1,ITP(11)
        IS=ict(j)
        DO 502 I=1,PTS
 
        A=DEXP(-EXX(J)*R2(I,is))
        B=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=B*DX(I,is)
        CHIX(I,J)=(THREE + SUM(J)*DX(I,is)*DX(I,is))*B
        CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=SIX*A*DX(I,is)+
     1             (NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
 502    CONTINUE
 501    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 511 J=ITP(11)+1,ITP(12)
        IS=ict(j)
        DO 512 I=1,PTS
        A=DEXP(-EXX(J)*R2(I,is))
        B=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=B*DY(I,is)
        CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(THREE + SUM(J)*DY(I,is)*DY(I,is))*B
        CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=SIX*A*DY(I,is)+
     1             (NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
C
 512    CONTINUE
 511    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 521 J=ITP(12)+1,ITP(13)
        IS=ict(j)
        DO 523 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=B*DZ(I,is)
        CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(THREE + SUM(J)*DZ(I,is)*DZ(I,is))*B
        CHID2(I,J)=SIX*A*DZ(I,is)+
     1             (NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
C
 523    CONTINUE
 521    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 531 J=ITP(13)+1,ITP(14)
        IS=ict(j)
        DO 532 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BXX=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=BXY*DX(I,is)
        CHIX(I,J)=(TWO + SUM(J)*DX(I,is)*DX(I,is))*BXY
        CHIY(I,J)=(ONE + SUM(J)*DY(I,is)*DY(I,is))*BXX
        CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=TWO*A*DY(I,is)+
     1             (NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
C
 532    CONTINUE
 531    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 541 J=ITP(14)+1,ITP(15)
        IS=ict(j)
        DO 543 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BXZ=DX(I,is)*DZ(I,is)*A
        BXX=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=BXZ*DX(I,is)
        CHIX(I,J)=(TWO + SUM(J)*DX(I,is)*DX(I,is))*BXZ
        CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,is)*DZ(I,is))*BXX
        CHID2(I,J)=TWO*A*DZ(I,is)+
     1             (NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
C
 543    CONTINUE
 541    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 561 J=ITP(15)+1,ITP(16)
        IS=ict(j)
        DO 563 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BYZ=DZ(I,is)*DY(I,is)*A
        BYY=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=BYZ*DY(I,is)
        CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(TWO + SUM(J)*DY(I,is)*DY(I,is))*BYZ
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,is)*DZ(I,is))*BYY
        CHID2(I,J)=TWO*A*DZ(I,is)+
     1             (NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
C
 563    CONTINUE
 561    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 551 J=ITP(16)+1,ITP(17)
        IS=ict(j)
        DO 552 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BYY=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=BXY*DY(I,is)
        CHIX(I,J)=(ONE + SUM(J)*DX(I,is)*DX(I,is))*BYY
        CHIY(I,J)=(TWO + SUM(J)*DY(I,is)*DY(I,is))*BXY
        CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=TWO*A*DX(I,is)+
     1             (NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
C
 552    CONTINUE
 551    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 571 J=ITP(17)+1,ITP(18)
        IS=ict(j)
        DO 572 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BXZ=DZ(I,is)*DX(I,is)*A
        BZZ=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=BXZ*DZ(I,is)
        CHIX(I,J)=(ONE + SUM(J)*DX(I,is)*DX(I,is))*BZZ
        CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(TWO + SUM(J)*DZ(I,is)*DZ(I,is))*BXZ
        CHID2(I,J)=TWO*A*DX(I,is)+
     1             (NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
C
 572    CONTINUE
 571    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 581 J=ITP(18)+1,ITP(19)
        IS=ict(j)
        DO 583 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BYZ=DZ(I,is)*DY(I,is)*A
        BZZ=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=BYZ*DZ(I,is)
        CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(ONE + SUM(J)*DY(I,is)*DY(I,is))*BZZ
        CHIZ(I,J)=(TWO + SUM(J)*DZ(I,is)*DZ(I,is))*BYZ
        CHID2(I,J)=TWO*A*DY(I,is)+
     1             (NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
C
 583    CONTINUE
 581    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 591 J=ITP(19)+1,ITP(20)
        IS=ict(j)
        DO 592 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BYZ=DZ(I,is)*DY(I,is)*A
        BXZ=DX(I,is)*DZ(I,is)*A
C
        CHI(I,J)=DX(I,is)*BYZ
        CHIX(I,J)=(ONE + SUM(J)*DX(I,is)*DX(I,is))*BYZ
        CHIY(I,J)=(ONE + SUM(J)*DY(I,is)*DY(I,is))*BXZ
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,is)*DZ(I,is))*BXY
        CHID2(I,J)=(NINE+SUM(J)*R2(I,is))*CHI(I,J)*SUM(J)
C
 592    CONTINUE
 591    CONTINUE
C
       DO 103 J=1,Nprims
       temp=zero
       DO 104 I=1,PTs
       check=DMax1(Dabs(CHi(I,J)),Dabs(CHiX(I,j)),Dabs(CHIY(I,J)),
     + Dabs(CHIZ(I,j)),Dabs(CHID2(I,J)),Temp)
       If(Check.gt.temp)temp=check
104    Continue
       chimax(j)=temp
103    Continue
C
       DO 105 L = 1,NMO
       DO 107 I=1,PTS 
       PSI(I,L) = ZERO
       GX(I,L) = ZERO
       GY(I,L) = ZERO
       GZ(I,L) = ZERO
       D2(I,L) = ZERO
107    CONTINUE 
C
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
126    CONTINUE
       ENDIF
125    CONTINUE
105    CONTINUE
C
       GOTO 999
C
133     CONTINUE
C
      DO 410 J = 1,NCENT
       DO 412 I=1,PTS 
        DX(I,J) = XYZ(I,1) - XC(J)
        DY(I,J) = XYZ(I,2) - YC(J)
        DZ(I,J) = XYZ(I,3) - ZC(J)
        R2(I,J)= DX(I,J)*DX(I,J)+DY(I,J)*DY(I,J)+DZ(I,J)*DZ(I,J)
412   CONTINUE
410   CONTINUE
C
C       FOR S-TYPE
C
        DO 420 J = 1,ITP(1) 
        is=ict(j)
        DO 422 I=1,PTS 
        CHI(I,J)=DEXP(-EXX(J)*R2(I,is))
422     CONTINUE
420     CONTINUE
C
C       FOR Px-TYPE
C
        DO 440 J=ITP(1)+1,ITP(2)
        is=ict(j)
        DO 442 I=1,PTS
        CHI(I,J)=DX(I,is)*DEXP(-EXX(J)*R2(I,is))
442     CONTINUE
440     CONTINUE
C
C       FOR Py-TYPE
C
        DO 460 J=ITP(2)+1,ITP(3)
        is=ict(j)
        DO 462 I=1,PTS
C
        CHI(I,J)=DY(I,is)*DEXP(-EXX(J)*R2(I,is))
462     CONTINUE
460     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 480 J=ITP(3)+1,ITP(4)
        is=ict(j)
        DO 482 I=1,PTS
        CHI(I,J)=DZ(I,is)*DEXP(-EXX(J)*R2(I,is))
482     CONTINUE
480     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 520 J=ITP(4)+1,ITP(5)
        is=ict(j)
        DO 522 I=1,PTS
        CHI(I,J)=DX(I,is)*DX(I,is)*DEXP(-EXX(J)*R2(I,is))
522     CONTINUE
520     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 540 J=ITP(5)+1,ITP(6)
        is=ict(j)
        DO 542 I=1,PTS
        CHI(I,J)=DY(I,is)*DY(I,is)*DEXP(-EXX(J)*R2(I,is))
542     CONTINUE
540     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 560 J=ITP(6)+1,ITP(7)
        is=ict(j)
        DO 562 I=1,PTS
        CHI(I,J)=DZ(I,is)*DZ(I,is)*DEXP(-EXX(J)*R2(I,is))
562     CONTINUE
560     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 580 J=ITP(7)+1,ITP(8)
        is=ict(j)
        DO 582 I=1,PTS
        CHI(I,J)=DX(I,is)*DY(I,is)*DEXP(-EXX(J)*R2(I,is))
582     CONTINUE
580     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 620 J=ITP(8)+1,ITP(9)
        is=ict(j)
        DO 622 I=1,PTS
        CHI(I,J)=DX(I,is)*DZ(I,is)*DEXP(-EXX(J)*R2(I,is))
622     CONTINUE
620     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 640 J=ITP(9)+1,ITP(10)
        is=ict(j)
        DO 642 I=1,PTS
        CHI(I,J)=DY(I,is)*DZ(I,is)*DEXP(-EXX(J)*R2(I,is))
642     CONTINUE
640     CONTINUE
C
C       FOR Fxxx-TYPE
C
        DO 503 J=ITP(10)+1,ITP(11)
        is=ict(j)
        DO 504 I=1,PTS
        CHI(I,J)=DEXP(-EXX(J)*R2(I,is))*DX(I,is)**3
 504    CONTINUE
 503    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 513 J=ITP(11)+1,ITP(12)
        is=ict(j)
        DO 514 I=1,PTS
        CHI(I,J)=DEXP(-EXX(J)*R2(I,is))*DY(I,is)**3
 514    CONTINUE
 513    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 524 J=ITP(12)+1,ITP(13)
        is=ict(j)
        DO 525 I=1,PTS
        CHI(I,J)=DEXP(-EXX(J)*R2(I,is))*DZ(I,is)**3
 525    CONTINUE
 524    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 533 J=ITP(13)+1,ITP(14)
        is=ict(j)
        DO 534 I=1,PTS
        CHI(I,J)=Dexp(-Exx(j)*r2(I,is))*DY(I,is)*DX(I,is)**2
 534    CONTINUE
 533    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 544 J=ITP(14)+1,ITP(15)
        is=ict(j)
        DO 545 I=1,PTS
        CHI(I,J)=Dexp(-Exx(j)*r2(I,is))*DZ(I,is)*DX(I,is)**2
 545    CONTINUE
 544    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 553 J=ITP(16)+1,ITP(17)
        is=ict(j)
        DO 554 I=1,PTS
        CHI(I,J)=Dexp(-Exx(j)*r2(I,is))*DX(I,is)*DY(I,is)**2
 554    CONTINUE
 553    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 564 J=ITP(15)+1,ITP(16)
        is=ict(j)
        DO 565 I=1,PTS
        CHI(I,J)=Dexp(-Exx(j)*r2(I,is))*DZ(I,is)*DY(I,is)**2
 565    CONTINUE
 564    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 573 J=ITP(17)+1,ITP(18)
        is=ict(j)
        DO 574 I=1,PTS
        CHI(I,J)=Dexp(-Exx(j)*r2(I,is))*DX(I,is)*DZ(I,is)**2
 574    CONTINUE
 573    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 584 J=ITP(18)+1,ITP(19)
        is=ict(j)
        DO 585 I=1,PTS
        CHI(I,J)=Dexp(-Exx(j)*r2(I,is))*DY(I,is)*DZ(I,is)**2
 585    CONTINUE
 584    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 593 J=ITP(19)+1,ITP(20)
        is=ict(j)
        DO 594 I=1,PTS
        CHI(I,J)=Dexp(-Exx(j)*r2(I,is))*DX(I,is)*DY(I,is)*DZ(I,is)
 594    CONTINUE
 593    CONTINUE
C
       DO 603 J=1,Nprims
       temp=zero
       DO 604 I=1,PTs
       check=DMax1(Dabs(CHi(I,J)),Temp)
       If(Check.gt.temp)temp=check
604    Continue
       chimax(j)=temp
603    Continue
C
       DO 605 L = 1,NMO
       DO 607 I=1,PTS 
       PSI(I,L) = ZERO
607    CONTINUE 
C
       DO 625 J = 1,NPRIMS
       check=dabs(chimax(j)*coo(J,L))
       IF(check.gt.cutoff)THEN
       TEMP=COO(J,L)
       DO 626 I = 1,PTS
       PSI(I,L) = PSI(I,L) + TEMP*CHI(I,J)
626    CONTINUE
       ENDIF
625    CONTINUE
605    CONTINUE
C
       GOTO 999
C
999    RETURN
       END

      SUBROUTINE GRDRHO (A,CX,XY)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL GRD
      PARAMETER(MCENT=200,MMO=500,MPRIMS=2000,MPTS=2000)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),NMO
      COMMON /ZZZZ/ PSI(MPTS,MMO),GX(MPTS,MMO),GY(MPTS,MMO),
     +GZ(MPTS,MMO),D2(MPTS,MMO)
      DIMENSION A(3,3),CX(3),XY(4),XYZ(MPTS,3),SUMR(MPTS)
      COMMON /GRID/ GRD(MPTS*MPTS)
      Save Zero,Two
      Data Zero/0.0d0/,Two/2.0d0/
C
C    CALCULATE NUMBER OF STEPS IN X AND Y FOR THIS GRID AS WELL AS
C    INCREMENT MARKERS FOR CENTERING PLOT
C
      IXSTP = IDINT(XY(1)/XY(2))
      IYSTP = IXSTP
      XCN = -XY(1)/Two
      YCN = XCN
      ZCN = Zero
      YCNT = XCN
      ZCnt=Zero
C
      DO 1400 IY = 1,IYSTP
C
C    RESET X INCREMENT
C
        XCNT = XCN
C
        DO 1300 I = 1,IXSTP
C
          SUMR(I)=Zero
C     BACK TRANSFORM PLANE INTO ORIGINAL MOLECULAR SYSTEM
C
          XYZ(I,1)=A(1,1)*XCNT+A(1,2)*YCNT+A(1,3)*ZCNT+CX(1)
          XYZ(I,2)=A(2,1)*XCNT+A(2,2)*YCNT+A(2,3)*ZCNT+CX(2)
          XYZ(I,3)=A(3,1)*XCNT+A(3,2)*YCNT+A(3,3)*ZCNT+CX(3)
C
        XCNT=XCNT+XY(2)
1300    CONTINUE
C
       ICALC = 1
       CALL GAUS(XYZ,IXSTP,ICALC)
C
        INDX=(IY-1)*IXSTP
        DO 904 J=1,NMO
        DO 905 I = 1,IXSTP
          SUMR(I)=SUMR(I)+PO(J)*PSI(I,J)*PSI(I,J)
905     CONTINUE
904     CONTINUE
C
        DO 906 I=1,IXSTP
        GRD(INDX+I)=SUMR(I)
906     CONTINUE 
C
C    INCREMENT Y VALUE
C
        YCNT=YCNT+XY(2)
C
1400  CONTINUE
C
      RETURN 
      END

      SUBROUTINE GRDD2R (A,CX,XY)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL GRD
      PARAMETER(MCENT=200,MMO=500,MPRIMS=2000,MPTS=2000)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),NMO
      COMMON /ZZZZ/ PSI(MPTS,MMO),GX(MPTS,MMO),GY(MPTS,MMO),
     +GZ(MPTS,MMO),D2(MPTS,MMO)
      DIMENSION A(3,3),CX(3),XY(4),XYZ(MPTS,3),SUMR(MPTS)
      COMMON /GRID/ GRD(MPTS*MPTS)
      Save Zero,Two
      Data Zero/0.0d0/,Two/2.0d0/
C
C    CALCULATE NUMBER OF STEPS IN X AND Y FOR THIS GRID AS WELL AS
C    INCREMENT MARKERS FOR CENTERING PLOT
C
      IXSTP = IDINT(XY(1)/XY(2))
      IYSTP = IXSTP
      XCN = -XY(1)/Two
      YCN = XCN
      ZCN = Zero
      YCNT = XCN
      ZCNT=ZEro
C
      DO 1400 IY = 1,IYSTP
C
C    RESET X INCREMENT
C
        XCNT = XCN
C
        DO 1300 I = 1,IXSTP
C
          SUMR(I)=Zero
C     BACK TRANSFORM PLANE INTO ORIGINAL MOLECULAR SYSTEM
C
          XYZ(I,1)=A(1,1)*XCNT+A(1,2)*YCNT+A(1,3)*ZCNT+CX(1)
          XYZ(I,2)=A(2,1)*XCNT+A(2,2)*YCNT+A(2,3)*ZCNT+CX(2)
          XYZ(I,3)=A(3,1)*XCNT+A(3,2)*YCNT+A(3,3)*ZCNT+CX(3)
C
        XCNT=XCNT+XY(2)
1300    CONTINUE
C
       ICALC = 0
       CALL GAUS(XYZ,IXSTP,ICALC)
C
        INDX=(IY-1)*IXSTP
        DO 904 J=1,NMO
        DO 905 I = 1,IXSTP
        SUMR(I)=SUMR(I)+Two*PO(J)*(PSI(I,J)*D2(I,J)
     +          +(GX(I,J)**2+GY(I,J)**2+GZ(I,J)**2))
905     CONTINUE
904     CONTINUE
C
        DO 906 I=1,IXSTP
        GRD(INDX+I)=SUMR(I)
906     CONTINUE
C 
C    INCREMENT Y VALUE
C
        YCNT=YCNT+XY(2)
C
 1400 CONTINUE
C
      RETURN 
      END

      SUBROUTINE GRDKEG (A,CX,XY)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL GRD
      PARAMETER (MCENT = 200, MMO = 500, MPRIMS = 2000, MPTS=2000)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),NMO
      COMMON /ZZZZ/ PSI(MPTS,MMO),GX(MPTS,MMO),GY(MPTS,MMO),
     +GZ(MPTS,MMO),D2(MPTS,MMO)
      DIMENSION A(3,3),CX(3),XY(4),XYZ(MPTS,3),SUMR(MPTS)
      COMMON /GRID/ GRD(MPTS*MPTS)
      Save Zero,Two,Pt5
      Data Zero/0.0d0/,Two/2.0d0/,Pt5/0.5d0/
C
C    CALCULATE NUMBER OF STEPS IN X AND Y FOR THIS GRID AS WELL AS
C    INCREMENT MARKERS FOR CENTERING PLOT
C
      IXSTP = IDINT(XY(1)/XY(2))
      IYSTP = IXSTP
      XCN = -XY(1)/Two
      YCN = XCN
      ZCN = Zero
      YCNT = XCN
      ZCnt=Zero
C
      DO 1400 IY = 1,IYSTP
C
C    RESET X INCREMENT
C
        XCNT = XCN
C
        DO 1300 I = 1,IXSTP
C
          SUMR(I)=Zero
C     BACK TRANSFORM PLANE INTO ORIGINAL MOLECULAR SYSTEM
C
          XYZ(I,1)=A(1,1)*XCNT+A(1,2)*YCNT+A(1,3)*ZCNT+CX(1)
          XYZ(I,2)=A(2,1)*XCNT+A(2,2)*YCNT+A(2,3)*ZCNT+CX(2)
          XYZ(I,3)=A(3,1)*XCNT+A(3,2)*YCNT+A(3,3)*ZCNT+CX(3)
C
        XCNT=XCNT+XY(2)
1300    CONTINUE
C
       CALL GAUS(XYZ,IXSTP,0)
C
        INDX=(IY-1)*IXSTP
        DO 904 J=1,NMO
        DO 905 I = 1,IXSTP
        SUMR(I)=SUMR(I)+Pt5*PO(J)*
     +  (GX(I,J)**2+GY(I,J)**2+GZ(I,J)**2)
905     CONTINUE
904     CONTINUE
C
        DO 906 I=1,IXSTP
        GRD(INDX+I)=SUMR(I)
906     CONTINUE
C
C    INCREMENT Y VALUE
C
        YCNT=YCNT+XY(2)
C
 1400 CONTINUE
C
      RETURN 
      END
