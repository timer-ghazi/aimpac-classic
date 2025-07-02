! Memory management routines for GRIDV

SUBROUTINE allocate_grid(grid_size)
   USE grid_data
   implicit none
   INTEGER, INTENT(IN) :: grid_size
   INTEGER :: ierr
   
   ! Deallocate if already allocated
   IF (ALLOCATED(GRD)) DEALLOCATE(GRD)
   
   ! Allocate grid array
   ALLOCATE(GRD(grid_size * grid_size), STAT=ierr)
   IF (ierr /= 0) STOP 'Memory allocation failed for GRD'
   
   ! Update size variable
   MPTS = grid_size
   
END SUBROUTINE allocate_grid

SUBROUTINE cleanup_grid()
   USE grid_data
   implicit none
   
   IF (ALLOCATED(GRD)) DEALLOCATE(GRD)
   
END SUBROUTINE cleanup_grid