The gridv program is a legacy computational chemistry utility written in FORTRAN 77. Its primary function is to calculate the value of scalar fields derived from a molecular wavefunction on a two-dimensional grid of points. Specifically, it computes the electron density (ρ(r)), its Laplacian (∇2 ρ(r)), or the kinetic energy density (G(r)).

The program is a component of the AIMPAC (Atoms in Molecules PACkage) suite, a set of tools based on the Quantum Theory of Atoms in Molecules (QTAIM) developed by R.F.W. Bader and his group. The output .grd file is intended for visualization programs (like contor or relief) to generate 2D contour or 3D surface plots, which are used to analyze the topological features of the electron density, such as bond critical points.

The code's structure is characteristic of its era (early 1990s), designed for performance on vector supercomputers. This is explicitly mentioned in the manual and is evident in the code's looping structures. Modernization efforts should focus on replacing outdated language features, improving data structures, and adapting the parallelism strategy to modern multi-core CPUs.

Our overall goal is to modernize the code without breaking any functionality.

Guiding Principle: Every single change, no matter how small, will be followed by a compilation and a rigorous validation against a baseline result. We will use version control (git) to checkpoint our progress and allow for easy rollbacks.

Already compiled and working legacy code is `./gridv_legacy` and the test files are `hco2.wfn` and `hco2.inf`, and the reference output produced by `./gridv_legacy` is `hco2.ref.grd`. The validation tool `./validate.py` compares the .grd file with the reference. Familiarize yourself with this tool. 

### **For Each Refactoring Step:**

1. **Perform the code modification** on `gridv.f90` and/or `modules.f90`.
2. **Run the test:** `make test`
  This single command will:
    - Recompile only the necessary files (`.f90` -> `.o`).
    - Link the objects to create the `./gridv_modern` executable.
    - Run the modern executable to create `hco2.grd`.
    - Run the validation script to compare `hco2.grd` with `hco2.ref.grd`.
3. **Analyze the output:**
    - If `make test` completes and the validation script reports `VALIDATION PASSED`, the change was successful. Commit it:
               git commit -am "Refactor: [Description of change]"
    - If `make test` fails (either during compilation or validation), the change was bad. Revert it:
                git reset --hard
      Then, re-attempt the refactoring step.
