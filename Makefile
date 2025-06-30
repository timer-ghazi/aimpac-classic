# Makefile for the gridv refactoring project

# --- Variables ---

# Compiler
FC = gfortran

# Compiler flags
# Basic flags for debugging and warnings
FFLAGS = -g -Wall -fcheck=all
# OpenMP flag for parallelization phase
FFLAGS_OMP = -fopenmp
# Optimization flags for final version
FFLAGS_OPT = -O3

# Linker flags (e.g., for LAPACK)
# This will be changed to -llapack during the LAPACK refactoring step
LDFLAGS =

# Source files
# Note: Legacy source is treated separately
SRC_LEGACY = gridv.f77
SRCS_MODERN = modules.f90 gridv.f90

# Object files are derived from modern sources
OBJS_MODERN = $(SRCS_MODERN:.f90=.o)

# Executable names
EXEC_LEGACY = gridv_legacy
EXEC_MODERN = gridv_modern

# Test case files
TEST_INF = hco2.inf
TEST_WFN = hco2.wfn
REF_GRD = hco2.ref.grd
TEST_GRD = hco2.grd
VALIDATE_SCRIPT = ./validate.py

# --- Targets ---

# Phony targets do not correspond to actual files
.PHONY: all clean test validate baseline

# The default target, executed when you just run 'make'
all: $(EXEC_MODERN)

# Rule to build the modern executable
$(EXEC_MODERN): $(OBJS_MODERN)
	$(FC) $(FFLAGS) $(FFLAGS_OMP) -o $@ $^ $(LDFLAGS)
	@echo "Modern executable '$(EXEC_MODERN)' is ready."

# Rule to build the legacy executable (for baseline creation)
$(EXEC_LEGACY): $(SRC_LEGACY)
	$(FC) -o $@ $^
	@echo "Legacy executable '$(EXEC_LEGACY)' is ready."

# Pattern rule to compile .f90 source files into .o object files
# The -c flag means 'compile only, do not link'
# This automatically handles dependencies (e.g., gridv.o depends on modules.mod)
%.o: %.f90
	$(FC) $(FFLAGS) $(FFLAGS_OMP) -c $< -o $@

# --- Workflow Targets ---

# Target to generate the reference/baseline grid file
baseline: $(EXEC_LEGACY)
	@echo "--- Generating baseline grid file... ---"
	./$(EXEC_LEGACY) $(TEST_INF) $(TEST_WFN)
	mv $(TEST_GRD) $(REF_GRD)
	@echo "Baseline file '$(REF_GRD)' created successfully."

# Target to run the full validation workflow
test: $(EXEC_MODERN)
	@echo "--- Running modern executable to generate test output... ---"
	./$(EXEC_MODERN) $(TEST_INF) $(TEST_WFN)
	@echo "--- Validating new grid against baseline... ---"
	@$(VALIDATE_SCRIPT) $(REF_GRD) $(TEST_GRD)

# 'validate' is an alias for 'test'
validate: test

# Target to clean up all generated files
clean:
	@echo "--- Cleaning up project directory... ---"
	rm -f $(EXEC_LEGACY) $(EXEC_MODERN) *.o *.mod $(TEST_GRD)
