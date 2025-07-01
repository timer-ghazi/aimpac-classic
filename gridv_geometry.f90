SUBROUTINE EULER (E,A)
!
   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!   CALCULATE THE ROTATION MATRICES USING
!   EULER ANGLES
!
   DIMENSION E(3), A(9)
!
   Save One,Two
   DATA One/1.0d0/,Two/2.0d0/
!
   PI=DACOS(-ONE)
!
   RADIAN=TWO*PI/360
!
   E1 = E(1)*RADIAN
   E2 = E(2)*RADIAN
   E3 = E(3)*RADIAN
!
   SA = DSIN(E1)
   SB = DSIN(E2)
   SC = DSIN(E3)
   CA = DCOS(E1)
   CB = DCOS(E2)
   CC = DCOS(E3)
!
   A(1) =  CC*CA - CB*SA*SC
   A(2) =  CC*SA + CB*CA*SC
   A(3) =  SC*SB
   A(4) = -SC*CA - CB*SA*CC
   A(5) = -SC*SA + CB*CA*CC
   A(6) =  CC*SB
   A(7) =  SB*SA
   A(8) = -SB*CA
   A(9) =  CB
!
   RETURN
END

SUBROUTINE GROCKLE (N, X, IR, S, E)
   IMPLICIT DOUBLE PRECISION (A-H, O-Z)
   DIMENSION X(3,N), S(3,N), C(3), E(3,3), EV(3), R(3,3)
   INTEGER IR(N)
   save Zero,One
   Data Zero/0.0d0/,One/1.0d0/
!
!    ZERO OUT CENTROID AND EIGENVECTORS
!
   C(1) = ZERO
   C(2) = ZERO
   C(3) = ZERO
   E(1,1) = ZERO
   E(1,2) = ZERO
   E(1,3) = ZERO
   E(2,2) = ZERO
   E(2,3) = ZERO
   E(3,3) = ZERO
!
!    CENTROID OF FRAGMENT
!
   M = ZERO
   DO 100 I = 1, N
      IF (IR(I) .GT. 0) THEN
         C(1) = C(1) + X(1,IR(I))
         C(2) = C(2) + X(2,IR(I))
         C(3) = C(3) + X(3,IR(I))
         M = M + 1
      ELSE IF (IR(I) .LT. 0) THEN
         C(1) = C(1) + S(1,-IR(I))
         C(2) = C(2) + S(2,-IR(I))
         C(3) = C(3) + S(3,-IR(I))
         M = M + 1
      END IF
100 CONTINUE
!
   DD = One/DFLOAT(M)
   C(1) = DD*C(1)
   C(2) = DD*C(2)
   C(3) = DD*C(3)
!
!    CALCULATE INERTIAL MATRIX.
!
   DO 200 I = 1, N
      IF (IR(I) .GT. 0) THEN
         X1 = X(1,IR(I)) - C(1)
         X2 = X(2,IR(I)) - C(2)
         X3 = X(3,IR(I)) - C(3)
      ELSE IF (IR(I) .LT. 0) THEN
         X1 = S(1,-IR(I)) - C(1)
         X2 = S(2,-IR(I)) - C(2)
         X3 = S(3,-IR(I)) - C(3)
      END IF
      E(1,1) = E(1,1) + X2**2 + X3**2
      E(2,2) = E(2,2) + X1**2 + X3**2
      E(3,3) = E(3,3) + X1**2 + X2**2
      E(1,2) = E(1,2) - X1*X2
      E(1,3) = E(1,3) - X1*X3
      E(2,3) = E(2,3) - X2*X3
      X1 = ZERO
      X2 = ZERO
      X3 = ZERO
200 CONTINUE
!
   E(2,1) = E(1,2)
   E(3,1) = E(1,3)
   E(3,2) = E(2,3)
!
!    GENERATES EIGENVALUES AND EIGENVECTORS OF THE INERTIAL MATRIX.
!
   CALL TRACE(E, EV, C, 3, IFAIL)
!
!    CHECK FOR RIGHT HAND CONVENTION FOR EIGEN-AXES
!
   DET = E(1,1)*(E(2,2)*E(3,3) - E(3,2)*E(2,3)) +&
   &E(1,2)*(E(3,1)*E(2,3) - E(2,1)*E(3,3)) +&
   &E(1,3)*(E(2,1)*E(3,2) - E(3,1)*E(2,2))
!
   IF (DET .LT. ZERO) THEN
      E(1,2) = -E(1,2)
      E(2,2) = -E(2,2)
      E(3,2) = -E(3,2)
   END IF
!
   RETURN
END