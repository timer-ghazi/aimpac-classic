#!/usr/bin/env python3
import sys
import os
import argparse
import numpy as np

# ANSI color codes for prettier output
C_RED = '\033[91m'
C_GREEN = '\033[92m'
C_YELLOW = '\033[93m'
C_BLUE = '\033[94m'
C_END = '\033[0m'

def print_fail(message):
    """Prints a failure message and exits."""
    print(f"{C_RED}❌ VALIDATION FAILED: {message}{C_END}")
    sys.exit(1)

def print_pass(message):
    """Prints a success message."""
    print(f"{C_GREEN}✅ {message}{C_END}")

def print_info(message):
    """Prints an informational message."""
    print(f"{C_BLUE}ℹ️  {message}{C_END}")

def parse_header_and_load(filepath):
    """
    Parses the .grd file header to get dimensions and number of centers,
    then loads the grid data.
    """
    with open(filepath, 'r') as f:
        # Line 1: Title (discard)
        f.readline()
        # Line 2: Blank (discard)
        f.readline()
        # Line 3: Grid dimensions
        dim_line = f.readline()
        nx = int(dim_line.strip().split()[0])
        # Line 4: Function type and number of centers
        ncent_line = f.readline()
        ncent = int(ncent_line.strip().split()[1])

    # The grid data starts after the header.
    # Header lines = 1(title) + 1(blank) + 1(dims) + 1(nfunc) + 1(center) + ncent(atoms)
    skip_rows = 5 + ncent
    grid_data = np.loadtxt(filepath, skiprows=skip_rows)
    return nx, grid_data


def main():
    """Main function to run the validation."""
    parser = argparse.ArgumentParser(
        description="Compares two AIMPAC .grd files for numerical consistency.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("reference_file", help="Path to the reference .grd file (the 'golden' standard).")
    parser.add_argument("new_file", help="Path to the new .grd file to be validated.")
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-8,
        help="Absolute tolerance for floating-point comparisons (default: 1e-8)."
    )
    args = parser.parse_args()

    print_info(f"Starting validation of '{args.new_file}' against '{args.reference_file}'")

    # --- 1. Pre-computation Sanity Checks ---
    for f in [args.reference_file, args.new_file]:
        if not os.path.exists(f):
            print_fail(f"File not found: '{f}'")

    if os.path.getsize(args.new_file) == 0:
        print_fail(f"New output file '{args.new_file}' is empty. The program likely crashed.")

    # --- 2. Header Validation and Data Loading ---
    try:
        nx_ref, ref_data = parse_header_and_load(args.reference_file)
        print_info(f"Reference grid dimensions: {nx_ref}x{nx_ref}")
    except Exception as e:
        print_fail(f"Could not parse reference file '{args.reference_file}': {e}")

    try:
        nx_new, new_data = parse_header_and_load(args.new_file)
        print_info(f"New grid dimensions: {nx_new}x{nx_new}")
    except Exception as e:
        print_fail(f"Could not parse new file '{args.new_file}': {e}")

    if nx_ref != nx_new:
        print_fail(f"Grid dimensions do not match! Reference is {nx_ref}x{nx_ref}, New is {nx_new}x{nx_new}.")

    if ref_data.shape != new_data.shape:
        print_fail(f"Data array shapes do not match! Reference is {ref_data.shape}, New is {new_data.shape}.")

    print_pass("File parsing and dimension checks passed.")

    # --- 3. Numerical Comparison ---
    print_info(f"Comparing numerical data with tolerance={args.tol}...")
    are_close = np.allclose(ref_data, new_data, atol=args.tol, rtol=0)

    if are_close:
        print_pass("Numerical comparison passed.")
        print(f"\n{C_GREEN}✅✅✅ VALIDATION PASSED ✅✅✅{C_END}")
        sys.exit(0)
    else:
        # Provide a bit more detail on failure
        diff = np.abs(ref_data - new_data)
        max_abs_diff = np.max(diff)
        mean_abs_diff = np.mean(diff)
        print_fail(
            "Numerical comparison failed.\n"
            f"       {C_YELLOW}Max absolute difference: {max_abs_diff:.2e}\n"
            f"       {C_YELLOW}Mean absolute difference: {mean_abs_diff:.2e}{C_END}"
        )


if __name__ == "__main__":
    main()