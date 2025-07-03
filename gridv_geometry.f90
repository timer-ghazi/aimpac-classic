SUBROUTINE EULER (E,A)
!
   implicit none
!
!   CALCULATE THE ROTATION MATRICES USING
!   EULER ANGLES
!
   ! Arguments
   real(8), intent(in) :: E(3)
   real(8), intent(out) :: A(9)
   
   ! Local variables
   real(8), save :: One = 1.0d0, Two = 2.0d0
   real(8) :: PI, RADIAN, E1, E2, E3, SA, SB, SC, CA, CB, CC
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
   implicit none
   
   ! Arguments  
   integer, intent(in) :: N
   real(8), intent(in) :: X(3,N), S(3,N)
   integer, intent(in) :: IR(N)
   real(8), intent(out) :: E(3,3)
   
   ! Local variables
   real(8) :: C(3), EV(3), R(3,3)
   real(8), save :: Zero = 0.0d0, One = 1.0d0
   integer :: I, M, IFAIL
   real(8) :: DD, X1, X2, X3, DET
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
   M = 0
   do i = 1, n
      if (ir(i) > 0) then
         c(1) = c(1) + x(1,ir(i))
         c(2) = c(2) + x(2,ir(i))
         c(3) = c(3) + x(3,ir(i))
         m = m + 1
      else if (ir(i) < 0) then
         c(1) = c(1) + s(1,-ir(i))
         c(2) = c(2) + s(2,-ir(i))
         c(3) = c(3) + s(3,-ir(i))
         m = m + 1
      end if
   end do
!
   dd = one/dfloat(m)
   c(1) = dd*c(1)
   c(2) = dd*c(2)
   c(3) = dd*c(3)
!
!    CALCULATE INERTIAL MATRIX.
!
   do i = 1, n
      if (ir(i) > 0) then
         x1 = x(1,ir(i)) - c(1)
         x2 = x(2,ir(i)) - c(2)
         x3 = x(3,ir(i)) - c(3)
      else if (ir(i) < 0) then
         x1 = s(1,-ir(i)) - c(1)
         x2 = s(2,-ir(i)) - c(2)
         x3 = s(3,-ir(i)) - c(3)
      end if
      e(1,1) = e(1,1) + x2**2 + x3**2
      e(2,2) = e(2,2) + x1**2 + x3**2
      e(3,3) = e(3,3) + x1**2 + x2**2
      e(1,2) = e(1,2) - x1*x2
      e(1,3) = e(1,3) - x1*x3
      e(2,3) = e(2,3) - x2*x3
      x1 = zero
      x2 = zero
      x3 = zero
   end do
!
   e(2,1) = e(1,2)
   e(3,1) = e(1,3)
   e(3,2) = e(2,3)
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
   if (det < zero) then
      e(1,2) = -e(1,2)
      e(2,2) = -e(2,2)
      e(3,2) = -e(3,2)
   end if
!
   RETURN
END