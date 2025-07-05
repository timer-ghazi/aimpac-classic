MODULE grid_data
   IMPLICIT NONE
   INTEGER :: MPTS  ! Now a variable, not parameter
   INTEGER :: CURRENT_GRID_SIZE = 0  ! Track actual allocated grid size
   REAL, ALLOCATABLE :: GRD(:)
END MODULE grid_data

MODULE molecular_data
   IMPLICIT NONE
   INTEGER, PARAMETER :: MCENT = 200
   INTEGER, PARAMETER :: MMO = 500
   DOUBLE PRECISION :: XC(MCENT), YC(MCENT), ZC(MCENT), CHARG(MCENT)
   DOUBLE PRECISION :: EORB(MMO), PO(MMO)
   INTEGER :: NCENT, NMO
END MODULE molecular_data

MODULE basis_data
   IMPLICIT NONE
   INTEGER, PARAMETER :: MPRIMS = 2000
   INTEGER, PARAMETER :: NTYPE = 20
   DOUBLE PRECISION :: DIV(MPRIMS), COO(MPRIMS,500), EXX(MPRIMS), SUM(MPRIMS)
   INTEGER :: ICT(MPRIMS), ITP(NTYPE), NPRIMS
END MODULE basis_data

MODULE string_data
   IMPLICIT NONE
   CHARACTER*80 :: WFNTTL, JOBTTL
   CHARACTER*8 :: ATNAM(200)
   INTEGER :: NAT
END MODULE string_data

MODULE io_units
   IMPLICIT NONE
   INTEGER :: INPT, IOUT, IWFN, IDBG
END MODULE io_units

MODULE computation_data
   IMPLICIT NONE
   DOUBLE PRECISION :: THRESH1, THRESH2, GAMMA, TOTE
END MODULE computation_data

MODULE work_arrays
   IMPLICIT NONE
   DOUBLE PRECISION, ALLOCATABLE :: PSI(:,:), GX(:,:), GY(:,:), GZ(:,:), D2(:,:)

CONTAINS

   SUBROUTINE allocate_work_arrays(npts, nmo)
      INTEGER, INTENT(IN) :: npts, nmo
      INTEGER :: stat
      
      ! Deallocate arrays if already allocated
      IF (ALLOCATED(PSI)) DEALLOCATE(PSI)
      IF (ALLOCATED(GX)) DEALLOCATE(GX)
      IF (ALLOCATED(GY)) DEALLOCATE(GY)
      IF (ALLOCATED(GZ)) DEALLOCATE(GZ)
      IF (ALLOCATED(D2)) DEALLOCATE(D2)
      
      ! Allocate arrays with error checking
      ALLOCATE(PSI(npts,nmo), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating PSI array'
      
      ALLOCATE(GX(npts,nmo), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating GX array'
      
      ALLOCATE(GY(npts,nmo), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating GY array'
      
      ALLOCATE(GZ(npts,nmo), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating GZ array'
      
      ALLOCATE(D2(npts,nmo), STAT=stat)
      IF (stat /= 0) STOP 'Error allocating D2 array'
   END SUBROUTINE allocate_work_arrays

   SUBROUTINE cleanup_work_arrays()
      ! Deallocate all work arrays
      IF (ALLOCATED(PSI)) DEALLOCATE(PSI)
      IF (ALLOCATED(GX)) DEALLOCATE(GX)
      IF (ALLOCATED(GY)) DEALLOCATE(GY)
      IF (ALLOCATED(GZ)) DEALLOCATE(GZ)
      IF (ALLOCATED(D2)) DEALLOCATE(D2)
   END SUBROUTINE cleanup_work_arrays

END MODULE work_arrays