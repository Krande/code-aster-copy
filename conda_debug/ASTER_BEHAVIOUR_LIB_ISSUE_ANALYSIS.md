# ASTER_BEHAVIOUR_LIB Environment Variable Issue - Analysis & Solution

**Date:** 2025-11-07  
**Issue:** Missing `ASTER_BEHAVIOUR_LIB` environment variable causing 20 MFront/MGIS test failures

---

## Problem Summary

20 tests are failing with the error:
```
MGIS_CONFIG_ERROR ('MFront/MGIS configuration error - ASTER_BEHAVIOUR_LIB not set')
```

**Affected Tests:**
```
ssnv518a, hsna106a, comp012f, mfron02j, ssnv163e, hsnv135a, ssnv205c, comp012g, 
ssnv232c, ssnv205b, hsnv136d, comp012h, ssnv520a, comp012i, zzzz509f, wtnv122e, 
ssnv519a, ssnv163d, ssnv181b, wtnv122d
```

---

## Root Cause Analysis

### 1. How ASTER_BEHAVIOUR_LIB is Set (Build Time)

In the main `wscript` file, the `check_variant_vars()` function sets the library name:

```python
@Configure.conf
def check_variant_vars(self):
    self.setenv("debug")
    self.env["ASTER_BEHAVIOUR_LIB"] = "AsterMFrOfficialDebug"
    self.define("ASTER_BEHAVIOUR_LIB", self.env["ASTER_BEHAVIOUR_LIB"])

    self.setenv("release")
    self.env["ASTER_BEHAVIOUR_LIB"] = "AsterMFrOfficial"
    self.define("ASTER_BEHAVIOUR_LIB", self.env["ASTER_BEHAVIOUR_LIB"])
```

This value is written to `aster_config.py` during the build process via `write_config_h("Python", variant)`.

### 2. How ASTER_BEHAVIOUR_LIB is Used (Runtime)

In `code_aster/Helpers/MGISHelper.py`, line 45:

```python
@classmethod
def from_embedded(cls, behaviour_name):
    """Create MGISBehaviour for an embedded behaviour.

    Arguments:
        behaviour_name (str): Name of the behaviour.
    """
    libpath = osp.join(
        os.environ["ASTER_LIBDIR"], "lib" + config["ASTER_BEHAVIOUR_LIB"] + ".so"
    )
    return MGISBuilder.from_library(behaviour_name, libpath)
```

The code directly accesses `config["ASTER_BEHAVIOUR_LIB"]` which will raise a `KeyError` if the key doesn't exist.

### 3. The Conda/Portability Issue

The current approach relies on:
1. `profile.sh` being sourced to set environment variables
2. The `config_env()` function in `data/wscript` which exports variables prefixed with `CONFIG_ENV_`

This is problematic for conda builds because:
- Conda packages should not rely on shell scripts (profile.sh) for configuration
- Environment variables set via profile.sh are not portable across different environments
- The `config` dictionary in Python should contain the value but may not in certain build scenarios

---

## Proposed Solutions

### Solution 1: Hardcoded Fallback in MGISHelper.py (Recommended)

**Advantages:**
- Simple, localized change
- Maintains backward compatibility
- Works in all conda/conda-forge environments
- No build system changes needed

**Implementation:**

Modify `code_aster/Helpers/MGISHelper.py`, line 45:

```python
@classmethod
def from_embedded(cls, behaviour_name):
    """Create MGISBehaviour for an embedded behaviour.

    Arguments:
        behaviour_name (str): Name of the behaviour.
    """
    # Use hardcoded fallback for conda builds where config might not have this key
    # The library name follows the pattern: AsterMFrOfficial (release) or AsterMFrOfficialDebug (debug)
    behaviour_lib = config.get("ASTER_BEHAVIOUR_LIB", "AsterMFrOfficial")
    libpath = osp.join(
        os.environ["ASTER_LIBDIR"], "lib" + behaviour_lib + ".so"
    )
    return MGISBuilder.from_library(behaviour_name, libpath)
```

**Note:** This assumes release builds. For debug builds, additional logic might be needed to detect the debug mode.

### Solution 2: Enhanced Fallback with Debug Detection

**Implementation:**

```python
@classmethod
def from_embedded(cls, behaviour_name):
    """Create MGISBehaviour for an embedded behaviour.

    Arguments:
        behaviour_name (str): Name of the behaviour.
    """
    # Determine library name with fallback
    if "ASTER_BEHAVIOUR_LIB" in config:
        behaviour_lib = config["ASTER_BEHAVIOUR_LIB"]
    else:
        # Fallback: detect debug mode by checking if debug libraries exist
        aster_libdir = os.environ.get("ASTER_LIBDIR", "")
        debug_lib = osp.join(aster_libdir, "libAsterMFrOfficialDebug.so")
        release_lib = osp.join(aster_libdir, "libAsterMFrOfficial.so")
        
        if osp.exists(debug_lib):
            behaviour_lib = "AsterMFrOfficialDebug"
        elif osp.exists(release_lib):
            behaviour_lib = "AsterMFrOfficial"
        else:
            # Final fallback to release
            behaviour_lib = "AsterMFrOfficial"
    
    libpath = osp.join(
        os.environ["ASTER_LIBDIR"], "lib" + behaviour_lib + ".so"
    )
    return MGISBuilder.from_library(behaviour_name, libpath)
```

### Solution 3: Ensure ASTER_BEHAVIOUR_LIB is Always in aster_config.py

This would require verifying that the build system always populates this value correctly during the waf configure/build process. However, this is less portable for conda builds.

---

## Recommended Action

**Implement Solution 1** - the simple hardcoded fallback:

1. **File to modify:** `code_aster/Helpers/MGISHelper.py`
2. **Line:** ~45 (in the `from_embedded` classmethod)
3. **Change:** Replace `config["ASTER_BEHAVIOUR_LIB"]` with `config.get("ASTER_BEHAVIOUR_LIB", "AsterMFrOfficial")`

This change:
- Is minimal and focused
- Provides a sensible default (release library)
- Maintains backward compatibility
- Works in both traditional and conda builds
- Follows Python best practices (using `.get()` with default values)

---

## Testing Plan

After implementing the fix:

1. Build the conda package
2. Run the 20 failing tests:
   ```bash
   ssnv518a, hsna106a, comp012f, mfron02j, ssnv163e, hsnv135a, ssnv205c, 
   comp012g, ssnv232c, ssnv205b, hsnv136d, comp012h, ssnv520a, comp012i, 
   zzzz509f, wtnv122e, ssnv519a, ssnv163d, ssnv181b, wtnv122d
   ```
3. Verify that the library is loaded correctly
4. Confirm no regression in existing tests

---

## Additional Notes

- The library naming convention is:
  - Release: `libAsterMFrOfficial.so`
  - Debug: `libAsterMFrOfficialDebug.so`
  
- These libraries are built by `mfront/wscript` when `BUILD_MFRONT` is enabled

- The `ASTER_LIBDIR` environment variable must still be set correctly (usually by conda activation scripts)

---

## Files Involved

- **Source of issue:** `code_aster/Helpers/MGISHelper.py` (line ~45)
- **Build configuration:** `wscript` (lines 609-616, `check_variant_vars` function)
- **Library build:** `mfront/wscript` (uses `env["ASTER_BEHAVIOUR_LIB"]` as target name)
- **Config generation:** `wscript` (lines 674-738, config file generation)

