! modules.f90 - Local arrays module for dynamic allocation in GAUS routines
! This module provides allocatable arrays for coordinate transformations and GAUS calculations

MODULE local_arrays
   IMPLICIT NONE
   
   ! Coordinate arrays for distance calculations (grid_points, atoms)
   DOUBLE PRECISION, ALLOCATABLE :: DX_LOCAL(:,:), DY_LOCAL(:,:), DZ_LOCAL(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: R2_LOCAL(:,:)
   
   ! Basis function value arrays (grid_points, basis_functions)
   DOUBLE PRECISION, ALLOCATABLE :: CHI_LOCAL(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: CHIX_LOCAL(:,:), CHIY_LOCAL(:,:), CHIZ_LOCAL(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: CHID2_LOCAL(:,:)
   
   ! Maximum basis function values (basis_functions)
   DOUBLE PRECISION, ALLOCATABLE :: CHIMAX_LOCAL(:)
   
CONTAINS

   SUBROUTINE allocate_local_arrays(npts, ncent, nprims)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: npts, ncent, nprims
      INTEGER :: stat
      
      ! Deallocate if already allocated
      CALL cleanup_local_arrays()
      
      ! Allocate coordinate arrays
      ALLOCATE(DX_LOCAL(npts, ncent), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating DX_LOCAL'
      
      ALLOCATE(DY_LOCAL(npts, ncent), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating DY_LOCAL'
      
      ALLOCATE(DZ_LOCAL(npts, ncent), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating DZ_LOCAL'
      
      ALLOCATE(R2_LOCAL(npts, ncent), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating R2_LOCAL'
      
      ! Allocate basis function arrays
      ALLOCATE(CHI_LOCAL(npts, nprims), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating CHI_LOCAL'
      
      ALLOCATE(CHIX_LOCAL(npts, nprims), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating CHIX_LOCAL'
      
      ALLOCATE(CHIY_LOCAL(npts, nprims), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating CHIY_LOCAL'
      
      ALLOCATE(CHIZ_LOCAL(npts, nprims), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating CHIZ_LOCAL'
      
      ALLOCATE(CHID2_LOCAL(npts, nprims), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating CHID2_LOCAL'
      
      ALLOCATE(CHIMAX_LOCAL(nprims), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating CHIMAX_LOCAL'
      
   END SUBROUTINE allocate_local_arrays
   
   SUBROUTINE cleanup_local_arrays()
      IMPLICIT NONE
      
      IF (ALLOCATED(DX_LOCAL)) DEALLOCATE(DX_LOCAL)
      IF (ALLOCATED(DY_LOCAL)) DEALLOCATE(DY_LOCAL)
      IF (ALLOCATED(DZ_LOCAL)) DEALLOCATE(DZ_LOCAL)
      IF (ALLOCATED(R2_LOCAL)) DEALLOCATE(R2_LOCAL)
      IF (ALLOCATED(CHI_LOCAL)) DEALLOCATE(CHI_LOCAL)
      IF (ALLOCATED(CHIX_LOCAL)) DEALLOCATE(CHIX_LOCAL)
      IF (ALLOCATED(CHIY_LOCAL)) DEALLOCATE(CHIY_LOCAL)
      IF (ALLOCATED(CHIZ_LOCAL)) DEALLOCATE(CHIZ_LOCAL)
      IF (ALLOCATED(CHID2_LOCAL)) DEALLOCATE(CHID2_LOCAL)
      IF (ALLOCATED(CHIMAX_LOCAL)) DEALLOCATE(CHIMAX_LOCAL)
      
   END SUBROUTINE cleanup_local_arrays

END MODULE local_arrays