# Windows (MSVC/ifx) patch — code review

Review of the minimal Windows support patch (`cf-win-patch` vs base `7825a46e08`, i.e.
release 18.1.7) as shipped to the conda-forge `code-aster-feedstock` as
`recipe/patches/windows_msvc_ifx.patch`.

Context at time of review: first full Windows run of the `submit` testcase suite gave
**2165 / 2313 passing, 148 failures**. The largest failure cluster (71 tests) is a native
access violation on the nonlinear solve path. An initial hypothesis that this was a MUMPS
struct ABI mismatch was **falsified** (5 of the 71 crashers use `LDLT`, never MUMPS; and a
rebuild with the MUMPS fix did not change the crash).

---

## 0. The nonlinear access-violation crashes: root-cause analysis

Two candidates both fit "71 nonlinear tests, identical crash right after
`VARI_ELGA initialisé a zéro`, no native frames, solver-independent". Discriminating them
takes minutes, not a rebuild.

### Candidate A (test first, ~5 min, no rebuild): stack overflow from ifx automatic arrays

- ifx on Windows puts automatic/temporary arrays on the **stack** by default; there is no
  `/heap-arrays` anywhere in the build (only `/integer-size:64 /real-size:64
  /names:lowercase /assume:underscore /MD`).
- The host process is `python.exe` with a ~2 MB stack reserve. On Linux `run_aster` runs
  under `ulimit -s unlimited`; there is no Windows equivalent in the patch, and you cannot
  relink `python.exe`.
- Release `bibfor` is compiled **`/Od`** (`bibfor/wscript:600-606`, `flags = ["/Od"]` in the
  *release* branch of `check_optimization_fcflags_msvc`), which massively amplifies stack
  usage — everything lives in memory slots, no reuse.
- Fit: the crash point is the first behaviour-integration/assembly (the deepest Fortran call
  chain in the code, `nmcomp`/te-routine stack), it is independent of MUMPS/LDLT, and a blown
  stack (`0xC00000FD`) is unwalkable — hence "no native frames". Linear controls
  (`sslp114a`) pass because their chains are shallower.

**Cheap verification, cheapest first:**

1. **Get the exception code** — this alone discriminates:
   `cdb -g -G -o python.exe -m run_aster ...` (or attach WinDbg, or check the Windows Event
   Log / WER reports for the crashed `python.exe`).
   `0xC00000FD` = stack overflow → Candidate A confirmed. `0xC0000005` = true AV →
   Candidate B.
2. **Patch the exe header, no rebuild**: back up `python.exe`, then
   `editbin /STACK:134217728 %PREFIX%\python.exe`, re-run `hsnv100n`. If it passes, done.
3. Permanent fix: add `/heap-arrays:0` (or `/heap-arrays:64`) to FCFLAGS in the recipe,
   and/or ship a `/STACK`-bumped launcher; also fix the `/Od` release bug (below), which
   reduces stack pressure for free.

### Candidate B: `aster_logical(kind=1)` ↔ C++ `bool` representation mismatch (ifx `.TRUE.` = 0xFF)

Evidence chain (all verified):

- `bibfor/include/asterf_types.h:29` — `#define aster_logical logical(kind=ASTER_LOGICAL_SIZE)`,
  with `ASTER_LOGICAL_SIZE=1`; `asterf_types.h:32` even bit-puns via `transfer`
  (`int_to_logical`).
- C++ side: `bibcxx/Supervis/astercxx.h:40` `using ASTERBOOL = bool;`
  `bibcxx/MemoryManager/JeveuxVector.h:890` `JeveuxVectorLogical = JeveuxVector<bool>` —
  **the same bytes** Fortran writes as `logical(kind=1)` are read as C++ `bool`, and vice
  versa.
- Concrete consumer on a hot path: `bibcxx/DataFields/SimpleFieldOnCells.h:442` `hasValue()`
  returns `(*_allocated)[position]` — the Fortran-written `.CESL` mask read as `bool`. C++
  writing `true` (0x01) into a mask later tested by ifx code using `.eqv. .TRUE.`
  (bit-compare against 0xFF) evaluates **false**.
- No mitigation flag exists: zero hits for `fpscomp` / `standard-semantics` in the tree or
  the build flags.
- ifx default: `.TRUE.` = 0xFF, truth test = low-bit/odd; clang-cl `bool` contract: value
  must be 0/1 (0xFF is UB; `b == true` compares against 1 and fails).

**Confidence:** HIGH (~85%) that this is the root cause of the *other* silent-failure family
— the 16 tests where `getValuesWithDescription()` / SimpleField extraction returns **empty**
(`ssnv177a`, `mtlp104a`, `ssls12a`, `zzzz505d/506j/509d/509r`, `sslx100a`, `ssna125a-d`,
`ssls142a/b`, `zzzz268a/b`; `<F> CALCULEL2_12` "CHAM_NO_S est vide" from
`bibfor/calculel/cnscno.F90:545` shows even Fortran sees an all-false mask).
**MODERATE (~40-50%)** for the AV crashers: the nonlinear path crosses this boundary heavily
and a wrong logical steering a branch/size could AV, but no specific corrupting call was
pinned — whereas Candidate A fits the crash signature (no frames, depth-correlated) better.
If the `cdb` exception code comes back `0xC0000005`, this becomes the prime suspect.

**Exact fix:** add **`/fpscomp:logicals`** to FCFLAGS (recipe-side, next to
`/integer-size:64`). It makes ifx emit `.TRUE.`=1 and treat any nonzero as true — fixing both
directions of the boundary. Do **not** use `/standard-semantics` (it drags in a dozen
unrelated semantic changes).

Side effects: `.eqv./.neqv.` and logical I/O become value-based (safe here); old saved bases
written with 0xFF logicals still read as true (nonzero); one residual risk is third-party
Fortran libs (MUMPS) receiving 1 instead of -1 — harmless, MUMPS tests nonzero. A source-side
fix (normalizing in `JeveuxVector<bool>` / `hasValue`) is strictly worse — dozens of crossing
points.

**Cheap verification without a full rebuild:**

1. In-situ probe on the **existing** install (no compile): in `python -m run_aster` console or
   a 5-line `.py` test, build any `FieldOnCells`, convert `toSimpleFieldOnCells()`, and dump
   the mask buffer: `np.frombuffer(sfield.toNumpy()[1], dtype=np.uint8)` (the mask array from
   `SimpleFieldOnCells.h:564` is `NPY_BOOL` over the raw `.CESL` bytes). Bytes = `0xFF` →
   hypothesis confirmed for the mask family.
2. Standalone 2-file repro (~2 min):
   `f.f90`: `subroutine setl(l) bind(c); use iso_c_binding; logical(c_bool)::l; l=.true._c_bool; end`
   — but compile it as plain `logical(kind=1)` (no `c_bool`) to mirror the code;
   `main.cpp`: `extern "C" void setl(bool*); bool b=false; setl(&b); printf("%d %d\n", (int)(unsigned char)b, b==true);`
   → `ifx /c f.f90 & clang-cl main.cpp f.obj` → prints `255 0` without `/fpscomp:logicals`,
   `1 1` with.
3. Targeted recompile is *not* possible for the real fix (the flag must apply to all of
   bibfor), so confirm via 1+2 before committing to the rebuild — and fold `/heap-arrays`,
   `/O2`-release, and the `iodr.c` fix into the same rebuild.

### Does it explain HHO / te0146 FPEs?

Mostly **separate**, with one connection. The 18 HHO failures + `ssls134d` are floating-point
traps (real divides by zero in `bibfor/hho/HHO_inertia_module.F90` and `bibfor/te/te0146.F90:538`),
enabled by `inisig.c:97-102` unmasking `_EM_ZERODIVIDE|_EM_OVERFLOW` on MSVC64. Two
Windows-specific aggravators in the patch:

- `bibc/include/aster_depend.h:116` defines `ASTER_HAVE_SUPPORT_FPE` only for LINUX/MINGW →
  **`matfpe()` is a silent no-op on MSVC64**, so all 112 existing FPE shields around
  BLAS/LAPACK are dead — this is also the likely cause of `hrom001a/b` crashing inside MKL
  `eigh`. Fix: add `ASTER_PLATFORM_WINDOWS` to that guard.
- If the logical-mask bug feeds zeroed data into these elements, the divides may vanish once
  Candidate B is fixed — retest HHO after, before excluding anything.

---

## 1. Ranked findings from the whole-patch review

### CRITICAL (silent data corruption)

1. **`bibc/utilitai/iodr.c`** — Jeveux base I/O broken twice on LLP64:
   (i) `:26/:29/:115/:211` — 8-byte `fread/fwrite` (`OFF_INIT`=8) into 4-byte `long nenr[]`
   elements: clobbers the neighbour slot and writes a garbage record-length header;
   (ii) `:145-148/:176` — `ASTER_HAVE_LONG_LONG` is only defined for MINGW
   (`aster_depend.h:74`; verified absent from the real build's `asterc_config.h`), so MSVC
   uses 32-bit `long offset` + `fseek` → any base >2 GB silently corrupts.
   Fix: `ASTERINTEGER`/`int64_t` + `_fseeki64` under `ASTER_PLATFORM_WINDOWS`.
2. **The `aster_logical`/`bool` mismatch** (above) — 16 tests currently return
   silently-empty/wrong data.
3. **`matfpe` no-op on Windows** (`aster_depend.h:116`) — FPE shields dead (above).

### HIGH

4. `bibfor/wscript:604` — **release Fortran built `/Od`**: the entire solver unoptimized.
   (Release C/C++ correctly `/O2`, `bibc/wscript:317`, `bibcxx/wscript`.) Also release lacks
   `/traceback`.
5. **`PyLong_FromLong((long)…)` truncation family** (Fortran→Python):
   `bibc/supervis/aster_module.c:1045,1168`, `aster_utils.c:394,416`,
   `aster_core_module.c:384,536-537` — jeveux integer vectors ("entiers codés" use high
   bits!) and integer `TEST_RESU` values truncated to 32 bits, silently.
   Fix: `PyLong_FromLongLong`.
6. `run_aster/utils.py:203-206` — on win, exit status is shifted `>>8` although `os.system`
   on Windows returns the code directly (the patch's own comment at `:218` says so, then
   keeps the shift): a raw exit code < 256 with no exitcode-file becomes 0 → failures
   reported as success.

### MEDIUM

7. **MINGW regressions from `MINGW→WINDOWS` guard swaps**: `bibc/utilitai/debugging.c:21`
   (MinGW now pulls `<execinfo.h>` — and its body at `:40` still tests MINGW, so MSVC selects
   the `backtrace()` branch, dormant only because `ASTER_HAVE_BACKTRACE` is undefined),
   `envima.c:41` (MinGW `ISMAX` falls back to 32-bit `LONG_MAX`), `libinfos.c:35`.
   Use the `mempid.c:38` pattern (`WINDOWS || MINGW`) consistently.
8. **Unconditional Linux behaviour changes in a "Windows" patch** — upstream will push back:
   `bibfor/echange/as_med_module.F90:76-81` (closes any Fortran unit holding the MED file, on
   all platforms, silently leaving the unit closed), `bibfor/op/op0039.F90:185-196` (same),
   `bibfor/supervis/lxinit.F90:114` (CR reclassified as blank everywhere),
   `code_aster/Cata/Commands/variable.py:24-32` (**byte-level `\r\n→\n` rewrite of a binary
   pickle stream on all platforms — can corrupt any pickle whose payload contains 0x0D0A**),
   `code_aster/MedUtils/MedConverter/field_converter.py:208,218` (unconditional int32 cast),
   `run_aster/run_aster_main.py:509` (`tee=False` for *all* Windows runs, not just ctest),
   top-level `wscript:291` (`add_os_flags("LDFLAGS")` now honoured on Linux too),
   `run_aster/export.py:231` (quote-stripping everywhere).
   Guard these `sys.platform == "win32"` or justify each.
9. `data/run_aster.bat:9-14` / `run_ctest.bat` — classic cmd.exe delayed-expansion bug:
   `%RUNASTER_ROOT%` is expanded when the `if not defined CONDA_PREFIX (...)` block is
   *parsed*, before the `set` inside it runs → the non-conda fallback path sets
   `ASTER_ROOT=\..` etc. Broken for anyone outside conda.
10. `data/profile.bat.tmpl` / `run_ctest.bat` — PYTHONPATH lacks the trailing `;.` that
    `profile.sh.tmpl:18` has → workdir imports fail (`supv001e`).
11. **Cross-DLL Fortran module variables** are not covered by the jeveux COMMON export scheme
    (`bibfor/include/jeveux.h:28-34` covers exactly its 6 COMMONs — complete for COMMONs;
    audit confirmed no other COMMON crosses). But three state-carrying modules sit in
    `bibfor.dll` with users in `bibfor_ext.dll`: `lmp_data_module` (live ref today at
    `third_party_interf/appcpr.F90:661`, harmless dead read), and
    `ldlt_xp_data_module.ap2foi_called` / `saddle_point_data_module` which **silently diverge
    the day PETSc is enabled** (writer `ap2foi.F90:72` in bibfor_ext, reader
    `nonlinear/nmresd.F90:139-141` in bibfor). Document or relocate.

### `msvc/` machinery (works, but fragile — flag to upstream reviewers)

12. **Export generation**: bibcxx uses a **positive allow-list**
    (`msvc/def_gen_cpp.py:66-89`) — exactly the mechanism that already silently dropped
    `BP*_tr6_Fortran` once; any new `extern "C"` in bibcxx fails only at consumer link time.
    bibfor/bibc use deny-lists (safer). `def_gen_c.py:78` has a real Python bug: missing comma
    → `"vsnprintf" "_scanf_l"` concatenates, so neither is excluded. `def_gen_c.py`'s
    `is_data` heuristic is effectively inverted (data lines contain "notype" → classified as
    functions), so **no data symbol except the 6 hardcoded jeveux COMMONs
    (`def_gen_fc.py:168-170`) is ever exported `DATA`** — the correct long-term approach is
    `__declspec(dllexport)` via the existing DEF*/CALL* macro layer, not dumpbin scraping.
    `run_dumpbin_for_file` (`def_gen_cache.py:209`) returns `(False, [])` on failure → an
    object whose dump fails silently contributes nothing. `msvc/msvc_lib.py:380-384`: bibfor's
    defgen subtracts symbols read from `msvc/bibfor_ext.def`/`asterGC.def` **with no ordering
    between the defgen tasks** → race/stale-file dependency; generated `.def`s and caches are
    written into the *source* tree (`create_defgen_task`, `msvc_lib.py:263-275`).
13. `msvc/mfront.def` — 46 committed **MSVC-STL-mangled** exports (`??$_Getvals@…`, RTTI
    descriptors): tied to one exact STL version, breaks on toolset bump; and as far as can be
    traced it is **dead** (only referenced via `clean_name_map`, no task ever produces a
    "mfront" libgen output). Same concern in miniature for `allowed_msvc_mangled` `toLower` in
    `def_gen_cpp.py:93-96`.
14. `config/ifort.py` / `config/msvc.py`: they take effect only via `@conf` name re-binding
    over stock waf 2.1.5 — a waf upgrade silently reverts individual functions; configure runs
    the whole detection **twice**, discarding the first pass (`wscript:274-277` loads before
    `setenv("default")`); `ifort.py:150` registry typo (`1AFortran`); `msvc.py:112-125` sorts
    VS versions by `float()` → 17.9 beats 17.14; `CC` env ignored (`msvc.py:903-913` forces
    `CC=CXX`); ~half of `msvc.py` is dead (WinCE/WinPhone/icl/libtool). Trim before
    upstreaming.
15. `msvc/c_entrypoints/entry_helpers.cxx:166-181` — computes the proxy-relative DLL path then
    loads by bare name anyway (relies on CONDA_PREFIX/PATH fallbacks); leaks `_dupenv_s`
    buffers; `SetDllDirectory(NULL)` side effect; `entry_helpers.h:1-3` personal attribution
    comment instead of the project header.
16. **Stale artifacts**: `data/update_pyd_links.bat` (hardlinks all `.pyd`s to `aster.dll` —
    the *old* single-DLL design, contradicts the proxy design, still installed by
    `data/wscript:324-328`); orphaned `bibfor/include/asterf_mumps.h` +
    `bibfor/include/mumps/*.h` (zero includes remain); `wscript:430` platform detection
    requires an external `DEFINES=ASTER_PLATFORM_MSVC64` env token — three inconsistent
    MSVC-detection mechanisms coexist (`wscript:360-361`, `:430`, `CC_NAME=='msvc'`).

### Genuine latent-bug fixes ifx exposed

All verified correct and platform-independent; **submit upstream as plain bugfixes, don't
bury them in the Windows patch**:

- `libs/gc/LibAsterGC.cxx` `long`→`int64_t` (correct; real stack corruption on LLP64)
- `bibfor/algeline/mstget.F90:119` `nb_mode_appui` init (correct)
- `bibfor/hho/HHO_matrix_module.F90:267` off-by-one slice (correct — the old slice wrote
  n+1 rows/cols)
- `bibfor/fracture/cgComputeGtheta.F90:170-171` (correct — callee declares
  `character(len=24)` per `exicp.h`, caller passed a len-0 literal → OOB read)
- `recursive` on `dismoi.F90` / `dismdy.F90` (correct)
- `bibfor/include/asterfort/chckVari.h` — the **base tree at 7825a46e08 contains literal
  merge-conflict markers** (verified via `git show`) — an upstream packaging bug worth
  reporting on its own
- `bibcxx` `stol`→`stoll` / `NPY_LONG`→`NPY_INT64` (correct; keep)
- `uthk.F90` `max()`-splitting and `mfront/VISC_ISOT_PLAS.mfront` `pascal`→`stress_unit`
  (windows.h macro collision) are legitimate but deserve a one-line comment saying why

---

## 2. `known_failures` split

### Must FIX before shipping (do not exclude)

- the 16 empty-mask tests (`ssnv177a`, `mtlp104a`, `ssls12a`, `zzzz505d`, `zzzz506j`,
  `zzzz509d`, `zzzz509r`, `sslx100a`, `ssna125a-d`, `ssls142a/b`, `zzzz268a/b`) → logical fix
- the 71 AV crashers → stack / logical fix
- `supv001e` → PYTHONPATH `.;`
- `distr01a` → Fortran `SHARE='DENYNONE'` on user units (`ulopen.F90`) or release user units
  in `ReservedUnitUsed`
- `zzzz264b/c` + `sdlv135c` → just apply the existing feedstock `numpy_py313_compat.patch` on
  Windows too
- `hrom001a/b` → the `matfpe`/`ASTER_HAVE_SUPPORT_FPE` fix (or `disable_fpe()` around
  `pod_analysis_base.py:650`)
- HHO 18 + `ssls134d` → re-triage after the logical + matfpe fixes

### Legitimately exclude (Windows `known_failures` list)

- HOMARD group (`forma01b/c`, `forma11a`, `tpll01j`, `tplp305d`, `ssnv173k`, `wtnl100c`,
  `zzzz121c/d/e`, `zzzz175b`, `zzzz319a/b`, `zzzz356a`) — HOMARD not packaged
- `sdls118a` (MISS3D)
- `func01e` (`/proc`)
- `vocab01b` (spawns `python3`)
- `asrun01b` (POSIX path asserts; also fix `run_aster/timer.py:58` `os.times()` →
  `time.monotonic()` for elapsed=0)
- `erreu09e` (upstream #34305 save-db semantics, not Windows)
- `ssnv187x` (MFront on-the-fly compile timeout)
- numerical-drift NOOKs `sdll159a`, `sdlv126e/f`, `sdnd107e`, `ssls134a` (all within a hair of
  tolerance — but re-check `ssls134a` after the logical fix)
- `mesh001m` (octree FP-branch platform sensitivity, 46/50 pass)

---

## 3. Recommendation on `e797886662` (MUMPS module fix)

**Keep it, but rewrite its commit message.**

Measured facts: with today's conda-forge MUMPS 5.8.2 headers (all kinds explicit), the struct
layout under the global promotion flags was already byte-identical to the library's — so it
never was an ABI mismatch *against these headers*, which is consistent with the rebuild not
fixing anything.

It is not, however, dead weight:

1. it future-proofs against any MUMPS header set with bare `INTEGER`/`REAL` (the
   `win-support` history suggests the original crash *was* against patched headers of that
   kind);
2. it removed the last COMMON blocks from `bibfor_ext` (alignment/exportability class of
   risk);
3. it is verified safe on both platforms (single state copy in `bibfor_ext.dll`; Linux
   `use`-list change is required by the module and net-equivalent).

Reverting buys ~80 lines of minimality at the cost of another rebuild and re-validation, and
reintroduces a known-fragile pattern. Just stop claiming it fixes the crashes, and note the
two genuine leftovers: delete or mark orphaned `bibfor/include/asterf_mumps.h` +
`bibfor/include/mumps/*.h`, and add `_mp_` module-data filtering to
`def_gen_bibfor_ext.py:95` before anyone regenerates the defs.

---

## Immediate action order

1. Get the crash exception code via `cdb` (minutes).
2. `editbin /STACK` test.
3. Mask-byte probe + 2-file `/fpscomp:logicals` repro.
4. One rebuild bundling: `/fpscomp:logicals`, `/heap-arrays:0`, release `/O2` + `/traceback`
   (`bibfor/wscript:604`), `iodr.c` LLP64 fixes, `ASTER_HAVE_SUPPORT_FPE` for WINDOWS,
   `PyLong_FromLongLong` family, PYTHONPATH `.`, numpy patch.
