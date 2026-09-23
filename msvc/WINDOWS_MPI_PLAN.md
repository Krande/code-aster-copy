# MPI support for the Windows (clang-cl + ifx) build — findings and plan

Branch `cf-win-mpi` (based on `cf-win-openmp` = 18.1.7 + cf-base + Windows patch + OpenMP).
Status: **Phase 1 (investigation only)**, 2026-09-23. Nothing below has been built yet;
the only experiments were tiny ifx compile/link checks and `mpiexec` smoke tests.

## 1. Package availability on conda-forge win-64 (anaconda.org API, 2026-09-23)

| Package | win-64 MPI builds? | Latest win-64 build(s) | Notes |
|---|---|---|---|
| `impi_rt` / `impi-devel` | yes (Intel MPI) | 2021.18.0 `h57928b3_749` / `hac47afa_749` | MPI 4.0, `mpi 1.0 impi`; `impi-devel` run-exports `impi_rt >=2021.18.0,<2021.18.1` (x.x.x pin) |
| `msmpi` | yes (MS-MPI) | 10.1.1 `h571195b_8` | MPI **2.2** headers (`MPI_VERSION 2`); depends on libgfortran |
| `mpich`, `openmpi` | **no** | – | 0 win-64 files |
| `mpi4py` | yes, both | 4.1.2 `py3xx…_101` for impi (`impi_rt >=2021.10`) and msmpi | loose impi pin, good |
| `hdf5` | **impi only** | 2.2.0 `mpi_impi_hbad3c58_0` (impi_rt 2021.17.2.*), 1.14.6 `mpi_impi_h8156f85_5` (impi_rt 2021.16.0.*) | hdf5-feedstock win: `mpi: [nompi, impi]`, built with ifx |
| `libmed` | **no** (nompi only) | 4.2 `nompi_*_29` | libmed-feedstock: `skip: mpi != "nompi" and (osx or win)`; nompi `medC.dll` has **no `MEDparFileOpen`** (checked exports) |
| `medcoupling` | **no** (nompi only) | 9.16.0 `py3xx_nompi_*_0` | depends `libmed * nompi_*` → cannot coexist with an MPI libmed |
| `libptscotch` / `ptscotch` | yes, impi and msmpi, int32/int64 | 7.0.13 `int64_ha3fd15e_0` (impi_rt **2021.16.0.***), 7.0.13 `int64_hd61e444_0` (msmpi) | libs: `ptscotch.lib ptscotcherr.lib …` |
| `mumps-mpi`, `mumps-mpi-fortran-devel` | **no** | – | mumps-feedstock win: `mpi: [nompi]` only (both flang and ifx `module_abi`) |
| `mumps-seq-fortran-devel` | (seq) | 5.8.2 `ifx_h4868554_3` → `mumps-seq 5.8.2 ifx_h94ad166_3` (MKL) | what the sequential Windows build uses |
| `scalapack` | **no** | – | but **MKL provides it**: `mkl_scalapack_lp64.2.dll`, `mkl_blacs_lp64.2.dll` + `mkl_blacs_intelmpi_lp64.2.dll` / `mkl_blacs_msmpi_lp64.2.dll` are in `mkl`/`mkl-devel` 2026.1 (already in our h_env); exports have `pdgesv_`/`blacs_gridinit_` lowercase+underscore |
| `parmetis` | **no** | – | optional for code_aster (`--disable-parmetis`, `--metis-libs=metis`) and for MUMPS |
| `petsc`, `slepc`, `petsc4py`, `slepc4py` | **no** | – | 0 win-64 files; PETSc-on-Windows is out of scope |

Solver check: `libptscotch=*=int64*` + `hdf5=*=mpi_impi*` + `mpi=*=impi` resolves to
`hdf5 1.14.6 mpi_impi_h8156f85_5` + `impi_rt 2021.16.0` (hdf5 2.x mpi_impi needs impi_rt 2021.17.2 →
conflict with current ptscotch). A ptscotch rebuild against impi 2021.17/18 would remove that constraint.

## 2. Recommended MPI implementation: Intel MPI (`impi-devel` / `impi_rt`)

Reasons:
- Only MPI with a **parallel HDF5** on win-64 (needed for parallel MED, see §5).
- MPI 4.0 complete. MS-MPI 10.1 lacks things code_aster uses:
  `MPI_Comm_create_group` (bibcxx/ParallelUtilities/MPIGroup.cxx) and `MPI_CXX_BOOL`
  (bibcxx/ParallelUtilities/AsterMPI.h) — checked against `msmpi.dll` exports / `mpi.h`.
  With impi all 42 MPI C functions used by bibc/bibcxx are exported (only the f2c/c2f macros "missing").
- Same vendor as ifx + MKL; MKL's BLACS dispatcher (`mkl_blacs_lp64_dll`) defaults to Intel MPI
  (`MKL_BLACS_MPI=INTELMPI`), MS-MPI would need `MKL_BLACS_MPI=MSMPI`.
- mpi4py and libptscotch have impi builds.

Downsides: `impi_rt` run-export pins x.x.x, so all MPI deps must be built against the same impi
patch release (currently ptscotch=2021.16, hdf5 2.x=2021.17.2, impi-devel latest=2021.18); impi.dll
exports COMMON blocks only in UPPERCASE (see §3).

## 3. Fortran ↔ MPI interface with ifx (`/names:lowercase /assume:underscore /integer-size:64`)

**code_aster itself never calls the MPI Fortran bindings.** All MPI calls go through C
(`bibc/supervis/aster_mpi.c`, `asmpi_*` wrappers declared in `bibfor/include/asterc/asmpi_*.h`,
communicators passed as `MPI_Fint` = `mpi_int` = `integer(4)`, converted with `MPI_Comm_f2c`).
`#include "mpif.h"` appears in 14 bibfor files, only for PARAMETER constants
(`MPI_SUM/MAX/MIN`, `MPI_INTEGER`, `MPI_INTEGER8`, `MPI_DOUBLE_PRECISION`, `MPI_DOUBLE_COMPLEX`)
converted via `to_mpi_int`. `ASTER_MPI_INT_SIZE` comes from `sizeof(MPI_Fint)` = 4 on impi.

Verified with ifx 2026.1 + impi 2021.18 `mpif.h` (tiny DLL link test in the scratchpad):
- constants-only use under `/integer-size:64 /names:lowercase /assume:underscore`: links, and the
  DLL has **no** import from impi.dll (unused DLLIMPORT COMMONs are not referenced). OK for code_aster.
- using `MPI_IN_PLACE` (COMMON `/MPIPRIV1/`): `LNK2019: unresolved external symbol __imp_mpipriv1_`.
  impi.dll exports Fortran *procedures* in all spellings (`mpi_allreduce`, `mpi_allreduce_`,
  `mpi_allreduce__`, `MPI_ALLREDUCE`) but COMMON blocks only as `MPIPRIV1`, `MPIPRIV2`, `MPIPRIVC`,
  `MPIFCMB5`, `MPIFCMB9`, `MPIFCMBA`. (msmpi.dll exports lowercase `mpipriv1_` too.)
- fix: append `!DEC$ ATTRIBUTES ALIAS:'MPIPRIV1' :: /MPIPRIV1/` (same for the 5 others) to a copy of
  `mpif.h` → links, imports `mpi_allreduce_` + `MPIPRIV1` from impi.dll.
  **This only matters for MUMPS** (5.8.2 uses `MPI_IN_PLACE` 112 times, no `MPI_STATUS_IGNORE`/`MPI_BOTTOM`).

C/C++: impi `mpi.h` includes the C++ bindings (`mpicxx.h`, needs `impicxx.lib`) unless
`MPICH_SKIP_MPICXX` is defined; code_aster does not use `MPI::` → define `MPICH_SKIP_MPICXX`.
`AsterMPI::mpi_type<>` is type-based (`long`→`MPI_LONG`, `long long`→`MPI_LONG_LONG`), so LLP64 is fine;
bibc uses `MPI_INTEGER8` for `ASTERINTEGER`.

MUMPS-MPI: code_aster calls MUMPS through its Fortran interface (`dmumps_struc.h`, `mumps%comm`),
so it needs an **ifx-ABI, `/names:lowercase /assume:underscore`, MKL, LP64** MUMPS built against impi,
with ScaLAPACK/BLACS from MKL (`mkl_scalapack_lp64_dll.lib mkl_blacs_lp64_dll.lib`).

## 4. Runtime (`run_aster` / `run_ctest` / `mpiexec`)

Smoke tests with conda `impi-devel 2021.18` on this machine (no hydra_service installed):
- `mpiexec -n 2 python -m mpi4py.bench helloworld` works (local launch needs no service/smpd).
- `PMI_RANK` / `PMI_SIZE` are set for each rank.
- `mpiexec -n 2 path\to\x.bat args` launches a .bat directly; exit code propagates (`exit /b 5` → 5).
- flags: `-prepend-rank`, `-ordered-output`, `--bind-to none`, `-localonly` accepted; `--tag-output` rejected.
- `mpiexec -n 1 echo x` **fails** (`echo` is a cmd builtin → "unable to create process");
  `hostname` works.
- singleton `MPI_Init` without mpiexec works (impi and msmpi, via mpi4py) → `require_mpiexec = 0`.

POSIX assumptions to fix:
- `data/wscript:check_mpiexec` probes with `echo` → all flag checks fail on Windows (would yield
  `mpiexec {program}` without `-n`). Use `hostname` on Windows; skip ssh/rsh and `--allow-run-as-root`.
- `data/wscript:check_mpi_get_rank` runs `mpiexec env` and stores `echo ${PMI_RANK}` (sh syntax).
  Windows: `echo %PMI_RANK%` (unset → literal text → `int()` fails → -1, which is the intended
  "not under mpiexec" value), or better make `run_aster/run.py:get_procid()` read
  `PMI_RANK`/`OMPI_COMM_WORLD_RANK`/`SLURM_PROCID` from `os.environ` directly.
- `waftools/parallel.py:check_vmsize` compiles a `/proc` + `unistd.h` fragment when MPI is on
  (called from `check_memory_stats`) → must be skipped on Windows (and set `require_mpiexec` 0).
- `run_aster/run_aster_main.py` re-runs `osp.join(RUNASTER_ROOT, "bin", "run_aster")` under mpiexec →
  must be `run_aster.bat` on Windows (the .bat does `chcp 65001` + profile.bat per rank; fine).
- `run_aster/run.py` `_gen_cmd_file("cmd{idx}.sh")`/ddt paths are only for dry-run/DDT output — leave.
- `run_ctest` handles `mpi_nbcpu` via ctest `PROCESSORS`; no change expected.

## 5. Scope and blockers

Parallel testcases (`P testlist … parallel`, 340 of 4617; heuristic grep of their comm files):
151 need neither PETSc nor ParallelMesh, 89 need ParallelMesh (parallel MED) but not PETSc,
100 need PETSc.

Blockers, in order of impact:
1. **No `mumps-mpi` on win-64** (hard blocker). Must extend mumps-feedstock (local checkout
   `C:\Work\code\mumps-feedstock`, branch `win/ifx-mkl`): its `recipe/CMakeLists.txt` has a
   `WITH_MPI` option but no ScaLAPACK/BLACS linking and no impi handling.
2. **No MPI libmed on win-64**: `bibcxx/IOManager/MedFilePointer.cxx:openParallel` calls
   `MEDparFileOpen`, absent from nompi `medC.dll` → an MPI build does not even link against nompi libmed.
   Needs libmed-feedstock win `mpi: impi` (hdf5 `mpi_impi` exists) — or, as a stopgap, a configure
   check + `#ifdef` so `openParallel` throws when libmed is sequential.
3. **medcoupling** has no MPI win build and its nompi build requires `libmed * nompi_*` → once libmed
   is MPI, medcoupling must also be an impi build (ParaMEDMEM on Windows = extra port), or be dropped
   (it is optional at configure: `check_medcoupling` failure is reverted; only Coupling/some tests use it).
4. **No PETSc/SLEPc/petsc4py on win-64** → `--disable-petsc` (loses ~100 parallel tests; out of scope).
5. **impi_rt x.x.x pin skew** between ptscotch (2021.16), hdf5 2.x (2021.17.2) and impi-devel (2021.18).
6. **No ParMETIS** → `--disable-parmetis --metis-libs=metis`; MUMPS without `-Dparmetis`.

Minimal viable scope (Phase A): impi + MUMPS-MPI(ifx, MKL ScaLAPACK) + ptscotch, **sequential**
hdf5/libmed/medcoupling with a guarded `openParallel`, no PETSc, no ParMETIS. That covers distributed
MUMPS / parallel elementary computations on replicated meshes (~151 parallel tests + all sequential
tests run with 1 rank). Phase B adds libmed (+hdf5) `mpi_impi` → ParallelMesh (+~89 tests).

## 6. Implementation plan (ordered)

### Step 1 — MUMPS-MPI for Windows (mumps-feedstock, local package first)
- `recipe/conda_build_config.yaml`: add `- impi  # [win]` to `mpi`; skip `win and mpi != 'nompi' and
  module_abi == 'flang'` (no netlib ScaLAPACK on win).
- `recipe/recipe.yaml` (`mumps-mpi`, `mumps-mpi-fortran-devel`): on win host use `impi-devel` (not
  `${{ mpi }}`), `mkl-devel` (ScaLAPACK/BLACS), `libptscotch * int${{ scotch_intsize }}_*`, `metis`;
  no `parmetis`, no `scalapack` on win; build string keeps `ifx_` prefix; test with `mpiexec -n 2`.
- `recipe/CMakeLists.txt` `WITH_MPI` on WIN32+ifx: don't build `libseq`; include `%LIBRARY_PREFIX%/include`,
  link `impi.lib`, `mkl_scalapack_lp64_dll.lib mkl_blacs_lp64_dll.lib` (+ existing MKL BLAS/LAPACK);
  `-Dscotch -Dptscotch -Dmetis`; keep `-DAdd_`, `/names:lowercase /assume:underscore`, `/MD`.
- Generate `build/mpi_include/mpif.h` = impi `mpif.h` + 6 `!DEC$ ATTRIBUTES ALIAS:'<UPPER>' :: /<common>/`
  lines, first in the Fortran include path (fixed-form safe: directives start in column 1).
- `build-mumps.bat`: `-DWITH_MPI=ON` when `%mpi%` != nompi; `make_integers_explicit.py` as for seq;
  run simpletests via `mpiexec -n 2 … < input` (stdin goes to rank 0).
- Check `mumps_c.c` compiles with `MPI_Comm_f2c` (macro in impi `mpi.h`), and that `esmumps.lib`
  still comes from somewhere (seq build ships it; ptscotch package does not).
- Build locally (`rattler-build` → local channel), then PR upstream.

### Step 2 — code_aster source (this branch, becomes part of `windows_msvc_ifx.patch`)
1. `waftools/parallel.py:load_compilers_mpi`: `if Utils.is_win32:` skip `check_cfg(--showme…)`;
   set `INCLUDES_MPI` (`%I_MPI_ROOT%\include` or `LIBRARY_PREFIX/include`), `LIBPATH_MPI`,
   `LIB_MPI = ["impi"]`, `DEFINES_MPI = ["MPICH_SKIP_MPICXX", "OMPI_SKIP_MPICXX"]`. Keep
   `check_mpi_fortran_interface` (mpi_f08 check is optional).
2. `waftools/parallel.py:check_vmsize`: return early on Windows (no /proc); set `require_mpiexec` 0
   (or compile a Windows fragment and still call `check_require_mpiexec`).
3. `data/wscript:check_mpiexec` / `check_mpi_get_rank`: Windows probe program `hostname`; rank command
   `echo %PMI_RANK%`; skip ssh/rsh/root checks. Expected config:
   `mpiexec -n {mpi_nbcpu} --bind-to none -prepend-rank {program}` (decide on `-ordered-output`).
   Optionally `run_aster/run.py:get_procid()` → env-var lookup (portable, no subprocess).
4. `run_aster/run_aster_main.py`: `run_aster` → `run_aster.bat` when `RUNASTER_PLATFORM == "win"`.
5. `waftools/mathematics.py:detect_mkl`: with MPI on msvc, don't append plain `scalapack`
   (would fail and fall back to `detect_math_lib` → `_scalapack()` fails). Either nothing (MUMPS DLL
   carries ScaLAPACK) or `mkl_scalapack_lp64_dll` + `mkl_blacs_lp64_dll`.
6. `waftools/med_cfg.py` + `bibcxx/IOManager/MedFilePointer.cxx`: define `ASTER_HAVE_MED_PARALLEL`
   when `MEDparFileOpen` links; otherwise `openParallel` throws (Phase A only; upstreamable).
7. `msvc/def_gen_*.py`: re-check that no MPI COMMON/`__imp_` symbols leak into `.def` exports.
8. Nothing expected in `bibc/supervis/aster_mpi.c` (only `<dlfcn.h>` under `OPEN_MPI`).

### Step 3 — feedstock (`recipe/`, new branch off upstream/main, after Step 1 packages exist)
- `conda_build_config.yaml`: `mpi: [nompi, openmpi # [not win], impi # [win]]`.
- `recipe.yaml`: `skip: osx or (python_impl == 'pypy')` (drop the win MPI skip);
  win+impi host: `impi-devel`, `mpi4py`, `mumps-mpi-fortran-devel * ifx_*`, `libptscotch * int64_*`,
  hdf5/libmed/medcoupling `*nompi*` (Phase A) → `*mpi_impi*` (Phase B); run: `mpi4py` (impi_rt comes
  from `impi-devel` run_exports). `build: ${{ mpi }}` must not add a package named `impi` on win.
- `build.bat`: `--enable-mpi` when `%mpi%` != nompi (`--disable-petsc --disable-parmetis
  --metis-libs=metis`), `INCLUDES_MUMPS` without `include/mumps_seq`, keep CC/FC = clang-cl/ifx
  (not the impi `mpicc.bat` wrappers).
- package test: a few `mpi_nbcpu 2` testcases through `run_ctest` (+ sequential smoke set).

### Step 4 — dependency feedstocks (Phase B, upstream PRs)
- libmed-feedstock: win `mpi: impi` (currently lists `msmpi # [win]` but skips win MPI), host
  `hdf5 * mpi_impi_*` + `impi-devel`.
- scotch-feedstock: rebuild ptscotch against a current impi (align impi_rt with hdf5).
- medcoupling-feedstock: win impi variant (ParaMEDMEM) — larger effort; optional.

## 7. Open questions
- msmpi as a second win MPI variant? Needs code changes (MPI-2.2: `MPI_Comm_create_group`,
  `MPI_CXX_BOOL`) and has no parallel hdf5 — not recommended now.
- `-ordered-output` buffers each rank's output until exit — acceptable for run_aster logs?
- Does conda-forge accept `mumps-mpi` win builds that depend on `mkl` (ScaLAPACK) only for ifx?
- Should Phase A (stubbed `openParallel`) ever ship to conda-forge, or only be used locally until
  libmed `mpi_impi` exists?
- impi at runtime in an activated env sets `I_MPI_ROOT` (activate.d); check that `python -m
  run_aster…` works from a non-activated env (DLL search for `impi.dll`/`libfabric.dll` in `Library\bin`).
