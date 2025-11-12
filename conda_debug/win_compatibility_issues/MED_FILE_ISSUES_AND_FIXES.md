# MED File Issues and Fixes Summary

## Issues Identified

### 1. **HDF5 Floating-Point Initialization Failure** (Critical)
**Error:**
```
H5T__init_native_float_types(): can't save floating-point environment, 
errno = 0, error message = 'No error', Win32 GetLastError() = 0
```

**Root Cause:**
- HDF5 built with `HDF5_ENABLE_THREADSAFE:BOOL=ON` tries to save/restore the floating-point control word
- On Windows with MSVC, especially in debug mode, this operation can fail
- The debug CRT has stricter validation that rejects certain FPU states

**Fix Applied:**
- Added FPU initialization in `PyInit_aster()` in `bibc/supervis/aster_module.c`
- Sets floating-point control word to `_CW_DEFAULT` before HDF5 is first used
- Forces early HDF5 initialization to catch errors at module load time

### 2. **Missing Error Checking in MED File Opening** (Critical)
**Error:**
```
D:\a\artifacts\bld\rattler-build_libmed_1739718231\work\src\ci\MEDfileVersionOpen.c [72] : 
Erreur à l'ouverture du fichier fort.20
```

**Root Cause:**
- `MEDfileOpen()` and `MEDparFileOpen()` can fail (return negative file ID)
- Code_Aster didn't check the return value
- Continued execution with invalid file ID (-1) set `_isOpen = true`
- Led to access violations when trying to use the invalid file handle

**Fix Applied:**
- Added error checking in `bibcxx/IOManager/MedFilePointer.cxx`
- Throws exception with descriptive message if file opening fails
- Sets `_isOpen = false` on failure

### 3. **Missing Bounds Checking in Mesh Access** (High)
**Error:**
```
Exception Code: 0xC0000005 (Access Violation)
std::_Ptr_base<MedMesh>::_Incref at memory:1370
```

**Root Cause:**
- `MedFileReader::getMesh(index)` didn't validate the index
- When file opening fails, `_meshes` vector is empty
- Accessing `_meshes[0]` causes invalid pointer dereference
- In debug mode, this triggers heap corruption detection

**Fix Applied:**
- Added bounds checking in `bibcxx/IOManager/MedFileReader.h`
- Validates index range before accessing `_meshes` vector
- Throws exception with helpful error message including available mesh count

### 4. **Debug Heap Corruption** (High)
**Error:**
```
Debug Assertion Failed!
Expression: _CrtIsValidHeapPointer(block)
File: debug_heap.cpp
```

**Root Cause:**
- Memory allocated in one CRT being freed in another
- Common issue when debug and release builds are mixed
- HDF5 with thread safety allocates memory that may cross DLL boundaries
- MED library allocates/frees memory that Code_Aster tries to manage

**Mitigation:**
- Ensure all dependencies (HDF5, MED, Code_Aster) use same build type
- Conda recipe already specifies build type matching: `hdf5 * *${{ build_type }}*`
- The FPU initialization fix should prevent HDF5 from entering error states
- Proper error checking prevents code from continuing after failures

### 5. **Missing Line Info in Stack Traces** (Low)
**Issue:**
- Stack traces show "(no line info)" unless build directory is kept
- PDB files reference source files by absolute path

**Explanation:**
- This is expected MSVC behavior, not a bug
- PDB files contain absolute paths to source files
- If build directory is deleted, debugger can't find sources
- This is normal for distributed packages

**Options:**
1. Keep build directory for development debugging
2. Use Windows symbol server for distributed builds
3. Accept limitation for release builds

## Files Modified

### 1. `bibcxx/IOManager/MedFilePointer.cxx`
```cpp
// Added error checking after MEDfileOpen
_fileId = MEDfileOpen( filename.string().c_str(), medAccessMode );
if ( _fileId < 0 ) {
    _isOpen = false;
    throw std::runtime_error( "Failed to open MED file: " + filename.string() );
}

// Added error checking after MEDparFileOpen
_fileId = MEDparFileOpen( filename.string().c_str(), medAccessMode, comm, MPI_INFO_NULL );
if ( _fileId < 0 ) {
    _isOpen = false;
    throw std::runtime_error( "Failed to open MED file in parallel mode: " + filename.string() );
}
```

### 2. `bibcxx/IOManager/MedFileReader.h`
```cpp
// Added bounds checking in getMesh
MedMeshPtr getMesh( int index ) const {
    if ( index < 0 || index >= static_cast< int >( _meshes.size() ) ) {
        throw std::runtime_error( "Mesh index out of bounds: " + std::to_string( index ) +
                                  " (available meshes: " + std::to_string( _meshes.size() ) +
                                  ")" );
    }
    return _meshes[index];
};
```

### 3. `bibc/supervis/aster_module.c`
```cpp
// Added Windows-specific includes
#ifdef ASTER_PLATFORM_MSVC64
#include <float.h>
#ifdef ASTER_HAVE_HDF5
#include "hdf5.h"
#endif
#endif

// Added FPU and HDF5 initialization in PyInit_aster
#ifdef ASTER_PLATFORM_MSVC64
    unsigned int current_word = 0;
    _controlfp_s(&current_word, 0, 0);
    _controlfp_s(&current_word, _CW_DEFAULT, _MCW_EM | _MCW_RC | _MCW_PC);
    
#ifdef ASTER_HAVE_HDF5
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl >= 0) {
        H5Pclose(fapl);
    }
#endif
#endif
```

## Recommendations for HDF5 Build

### Issue: Thread-Safe HDF5 on Windows
The current HDF5 build uses `HDF5_ENABLE_THREADSAFE:BOOL=ON`, which causes issues on Windows:

**Option 1: Disable Thread Safety (Simplest)**
```bat
-D HDF5_ENABLE_THREADSAFE:BOOL=OFF ^
```
- Eliminates FPU state saving issues
- May impact multi-threaded performance
- Code_Aster primarily uses sequential HDF5 access

**Option 2: Use Non-Thread-Safe HDF5 for Debug Builds**
```bat
set HDF5_THREADSAFE=ON
if "%build_type%" == "debug" (
    set HDF5_THREADSAFE=OFF
)
-D HDF5_ENABLE_THREADSAFE:BOOL=%HDF5_THREADSAFE% ^
```

**Option 3: Keep Thread Safety but Add Workarounds** (Current Approach)
- Keep HDF5 thread-safe
- Add FPU initialization in Code_Aster (already done)
- May still have edge cases

## Testing

After applying these fixes, test with:
```cmd
pixi run buildd
pixi run installd
pixi run ctestd adlv100a
```

Expected behavior:
1. HDF5 should initialize without floating-point errors
2. MED file opening failures should produce clear error messages
3. No access violations or heap corruption
4. Helpful error messages if file doesn't exist or is invalid

## Additional Notes

### Why the file might not be opening:
1. **File doesn't exist** - Check if fort.20 is actually created before opening
2. **HDF5 initialization failed** - Fixed by FPU initialization
3. **File permissions** - Unlikely in test environment
4. **Invalid file format** - HDF5 version mismatch or corrupted file

### Next Steps:
1. Rebuild Code_Aster with the fixes
2. Run the test suite
3. Check if HDF5 error is gone
4. If file still doesn't open, the error message will now be clear
5. Consider disabling HDF5 thread safety if issues persist

