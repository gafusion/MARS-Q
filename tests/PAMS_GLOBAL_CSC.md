# Optional global steady-response solve

The default `PAMS_CSC=0` build retains the original PAMS block solver and has no
fortsparse dependency. Selecting an unavailable or unknown backend fails clearly.
The new path is limited to steady `NCASE=2`, `ALAM=1`, `INORMSOL=0`.

Build using the existing native makefile:

```sh
make clean
make PAMS_CSC=1 FORTSPARSE_PREFIX=/path/to/fortsparse F95=mpif90 \
  F95FLAGS='-O2 -fopenmp -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-none -ffree-line-length-none -fallow-argument-mismatch'
```

Use the native CMake build prefix (containing include/, lib/ and bin/) from
fortsparse commit5a8d93c2676260bc562d6233f2ae6e8ce8dcbccb or later with the
same UMFPACK IPC option/budget API, and its matching helper.
Supply `FORTSPARSE_LIBS` when a static library requires extra platform libraries.
The adapter uses the existing public concrete UMFPACK backend API and does not
link UMFPACK into MARS. At runtime:

```sh
export MARS_PAMS_BACKEND=UMFPACK_IPC
export FORTSPARSE_UMFPACK_HELPER=/path/to/fortsparse_umfpack_helper
export MARS_PAMS_UMFPACK_BUDGET_BYTES=... # positive symbolic peak estimate limit
export MARS_PAMS_ROW_ERRORS=original_row_errors.txt
```

The helper budget is a pre-numeric estimate gate, not a hard process memory cap.
Reserve separate memory for original dense blocks, caller CSC, IPC and vectors.
Do not set the factor budget equal to the whole job allocation. A budget failure
returns nonzero; there is no automatic fallback to the original block solver.
Pivot tolerance1 permits global numerical pivoting without an explicit inverse.

Native `UBXRHS` forms the identical forcing before elimination. The adapter reads
CALPAM's original block views: A/H have already been shifted by their caller.
It uses ordering X1,Y1,X2,Y2,...,Xlast, retains every exact nonzero coefficient,
and counts dimensions/nonzeros in int64 before casting to the public int32 CSC.
Original block arrays are not mutated. ALAM, ALNORM and complex phase are unchanged;
NIT=1 and NONCON=0 indicate a certified sparse solve. The default retains its
historical NCASE2 NONCON=3 semantics.

Every original row is checked using |b-Ax|/(|b|+sum|a_ij||x_j|), without a floor.
A zero denominator and zero residual give zero. Extended exponent-range
accumulation handles products of finite doubles, including tiny rows. Any
nonfinite result or row error above1e-11 aborts before a result is accepted.
A finite quality failure preserves the maximum row/value and optional full row
file before exiting nonzero. Failures without a valid certificate do not dump it.
This backward-error certificate alone does not establish small forward error for
an ill-conditioned physical matrix.

Run one focused native test with an existing library build:

```sh
python tests/test_pams_global_csc.py \
  --fortsparse-prefix /path/to/fortsparse --build-parent /path/to/scratch \
  --baseline-revision 8824bb18e1514b4a27f6357d5fe859d4fe690542 \
  --report /path/to/native-caller-before-after.json
```

The cases cover zero D, the nearly-zero-D failure family, unequal X/Y dimensions,
complex/scaled rows, both radial boundaries, a genuinely singular matrix and tiny
row/zero-row certificates. CSC action is compared with original native MUAMAT.
The optional baseline comparison compiles the same caller fixture against the
unmodified original pams.f with the adapter/library still present: zero and
near-zero local pivots fail the original-system oracle before the hook repair
and pass afterward. A Fortran STOP with status0 cannot pass without the required
completion marker. The full-matrix native LU oracle passes on both variants.

The CALPAM/PAMERA/UBXRHS lifecycle test uses an explicitly manufactured force
callback; it checks one RHS call, shifted guard views, complex phase, unchanged
ALAM/ALNORM, NIT/NONCON and dependency-free default behavior. Unrelated callbacks
abort if invoked. This fixture does not qualify the physical MUBMAT force producer.

Before production use, a held native qualification must run the actual physical
MUBMAT/LINEAR/CALPAM path with unchanged source equilibrium, profiles, coil input
and Fourier conventions. Record binary/source/library/helper pins, all original
row errors, symbolic/factor memory, actual peak RSS, and output readback closure.
No revised field or torque curve is admitted solely by this small test.
