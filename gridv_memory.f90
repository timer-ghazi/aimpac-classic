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
   
   ! Update size variables
   MPTS = grid_size
   CURRENT_GRID_SIZE = grid_size
   
   ! Note: work arrays need to be allocated separately after NMO is known
   
END SUBROUTINE allocate_grid

SUBROUTINE cleanup_grid()
   USE grid_data
   implicit none
   
   IF (ALLOCATED(GRD)) DEALLOCATE(GRD)
   
END SUBROUTINE cleanup_grid

SUBROUTINE allocate_all_arrays(grid_size, nmo)
   USE grid_data
   USE work_arrays
   implicit none
   INTEGER, INTENT(IN) :: grid_size, nmo
   
   ! Allocate grid array first
   CALL allocate_grid(grid_size)
   
   ! Allocate work arrays with proper error checking
   CALL allocate_work_arrays(grid_size, nmo)
   
END SUBROUTINE allocate_all_arrays

SUBROUTINE cleanup_all_arrays()
   USE grid_data
   USE work_arrays
   implicit none
   
   ! Clean up grid arrays
   CALL cleanup_grid()
   
   ! Clean up work arrays
   CALL cleanup_work_arrays()
   
   ! Reset grid size tracking
   CURRENT_GRID_SIZE = 0
   
END SUBROUTINE cleanup_all_arrays

FUNCTION get_current_grid_size() RESULT(grid_size)
   USE grid_data
   implicit none
   INTEGER :: grid_size
   
   grid_size = CURRENT_GRID_SIZE
   
END FUNCTION get_current_grid_size