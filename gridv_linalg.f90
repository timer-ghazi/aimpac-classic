SUBROUTINE      TQLGRM  (N, D, E, Z, IERR)
   implicit none
   
   ! Arguments
   integer, intent(in) :: N
   real(8), intent(inout) :: D(N), E(N), Z(N,N)
   integer, intent(out) :: IERR
   
   ! Parameters
   real(8), parameter :: AMACH = 16.0E-13
   real(8), parameter :: ZERO = 0.0D0, ONE = 1.0D0
   
   ! Local variables
   integer :: I, J, L, M, L1, II, MML, K
   real(8) :: F, B, H, G, P, R, C, S
!
   IERR    = 0
   IF (N .EQ. 1) RETURN
!
   do i = 2, n
      e(i-1) = e(i)
   end do

   F       = ZERO
   B       = ZERO
   E(N)    = ZERO

   do l = 1, n
      j = 0
      h = amach*(dabs(d(l)) + dabs(e(l)))
      if (b < h) b = h

      ! Find first small subdiagonal element
      do m = l, n
         if (dabs(e(m)) <= b) exit
      end do

      if (m /= l) then
         ! Iterate until convergence
         do
            if (j == 30) then
               ierr = l
               return
            end if

            j = j + 1
            l1 = l + 1
            g = d(l)
            p = (d(l1) - g)/(2*e(l))
            if (dabs(p*amach) > one) then
               r = p
            else
               r = dsqrt(p*p + 1)
            end if
            d(l) = e(l)/(p + dsign(r,p))
            h = g - d(l)

            do i = l1, n
               d(i) = d(i) - h
            end do

            f = f + h
            p = d(m)
            c = one
            s = zero
            mml = m - l

            do ii = 1, mml
               i = m - ii
               g = c*e(i)
               h = c*p
               if (dabs(p) >= dabs(e(i))) then
                  c = e(i)/p
                  r = dsqrt(c*c + 1)
                  e(i+1) = s*p*r
                  s = c/r
                  c = one/r
               else
                  c = p/e(i)
                  r = dsqrt(c*c + 1)
                  e(i+1) = s*e(i)*r
                  s = 1.d0/r
                  c = c*s
               end if
               p = c*d(i) - s*g
               d(i+1) = h + s*(c*g + s*d(i))

               do k = 1, n
                  h = z(k,i+1)
                  z(k,i+1) = s*z(k,i) + c*h
                  z(k,i) = c*z(k,i) - s*h
               end do
            end do

            e(l) = s*p
            d(l) = c*p
            if (dabs(e(l)) <= b) exit
         end do
      end if

      d(l) = d(l) + f
   end do

   ! Sort eigenvalues and eigenvectors
   do ii = 2, n
      i = ii - 1
      k = i
      p = d(i)

      do j = ii, n
         if (d(j) < p) then
            k = j
            p = d(j)
         end if
      end do

      if (k /= i) then
         d(k) = d(i)
         d(i) = p

         do j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
         end do
      end if
   end do
   RETURN
END

SUBROUTINE      TRACE   (H, E, W, N, IERR)
!
! TRACE CALLS TREDIG AND TLQGRM TO DIAGONALIZE A SYMMETRIC REAL MATRIX.
! THE MATRIX IS PASSED DOWN IN H AND IS REPLACED BY THE EIGENVECTORS.
! THE EIGENVALUES IN E ARE STORED SMALLEST FIRST.
! THE WORK STORE W SHOULD BE AT LEAST OF DIMENSION N.
! SKK ==================================================================
!
   implicit none
   
   ! Arguments
   integer, intent(in) :: N
   real(8), intent(inout) :: H(N,N), E(N), W(N)
   integer, intent(out) :: IERR
!
   CALL TREDIG     (N, E, W, H)
   CALL TQLGRM     (N, E, W, H, IERR)
!
   RETURN
END

SUBROUTINE      TREDIG  (N, D, E, Z)
   implicit none
   
   ! Arguments
   integer, intent(in) :: N
   real(8), intent(out) :: D(N), E(N)
   real(8), intent(inout) :: Z(N,N)
   
   ! Parameters
   real(8), parameter :: ZERO = 0.0D0, ONE = 1.0D0
   
   ! Local variables
   integer :: II, I, L, K, J, JP1
   real(8) :: H, SCALE, RSCALE, F, G, RH, RHSCALE, HH

   if (n == 1) then
      d(1) = zero
   else
      do ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = zero
         scale = zero

         if (l >= 2) then
            ! Calculate scale factor
            do k = 1, l
               scale = scale + dabs(z(i,k))
            end do

            if (scale /= zero) then
               rscale = one/scale
               do k = 1, l
                  z(i,k) = z(i,k)*rscale
                  h = h + z(i,k)*z(i,k)
               end do
               f = z(i,l)
               g = -dsign(dsqrt(h),f)
               e(i) = scale*g
               h = h - f*g
               z(i,l) = f - g
               f = zero
               rh = one/h
               rhscale = rh*rscale

               do j = 1, l
                  z(j,i) = z(i,j)*rhscale
                  g = zero

                  do k = 1, j
                     g = g + z(j,k)*z(i,k)
                  end do

                  jp1 = j + 1
                  if (l >= jp1) then
                     do k = jp1, l
                        g = g + z(k,j)*z(i,k)
                     end do
                  end if

                  e(j) = g*rh
                  f = f + e(j)*z(i,j)
               end do

               hh = f/(h + h)

               do j = 1, l
                  f = z(i,j)
                  g = e(j) - hh*f
                  e(j) = g
                  do k = 1, j
                     z(j,k) = z(j,k) - f*e(k) - g*z(i,k)
                  end do
               end do

               do k = 1, l
                  z(i,k) = scale*z(i,k)
               end do
            else
               e(i) = z(i,l)
            end if
         else
            e(i) = z(i,l)
         end if

         d(i) = h
      end do

      d(1) = zero
   end if
   e(1) = zero

   do i = 1, n
      l = i - 1
      if (d(i) /= zero) then
         do j = 1, l
            g = zero

            do k = 1, l
               g = g + z(i,k)*z(k,j)
            end do

            do k = 1, l
               z(k,j) = z(k,j) - g*z(k,i)
            end do
         end do
      end if

      d(i) = z(i,i)
      z(i,i) = one
      if (l >= 1) then
         do j = 1, l
            z(j,i) = zero
            z(i,j) = zero
         end do
      end if
   end do
   RETURN
END