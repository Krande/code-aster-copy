// GNU __atomic_* builtin stubs for MSVC linkers
// These stubs provide the minimal set used by MUMPS when compiled with
// compilers that emit GNU-style atomic builtins, while linking with MSVC.
//
// Notes:
// - Implementations below cover 4- and 8-byte loads/stores and a 32-bit
//   compare-exchange which appears sufficient for the current MUMPS usage
//   triggering unresolved externals like __atomic_load.
// - On 64-bit Windows, Interlocked*64 are used. On 32-bit, 64-bit ops are
//   not truly atomic here but are provided for completeness.
// - Memory order parameters are accepted for signature compatibility but
//   are not used to model precise ordering semantics.

#include <windows.h>
#include <stdint.h>
#include <stdbool.h>

// __atomic_compare_exchange(ptr, expected, desired, weak, success_memorder, failure_memorder)
// Returns true if *ptr was equal to *expected and the swap was performed
__declspec(dllexport) bool __atomic_compare_exchange(void* ptr, void* expected, const void* desired,
                               bool weak, int success_memorder, int failure_memorder) {
    (void)weak; (void)success_memorder; (void)failure_memorder;
    volatile LONG* p = (volatile LONG*)ptr;
    LONG* exp = (LONG*)expected;
    LONG des = *(const LONG*)desired;

    LONG old_val = InterlockedCompareExchange(p, des, *exp);

    if (old_val == *exp) {
        return true;
    } else {
        *exp = old_val;
        return false;
    }
}

// __atomic_load_8 and __atomic_load_4 return the value at *ptr
__declspec(dllexport) int64_t __atomic_load_8(const void* ptr, int memorder) {
    (void)memorder;
    volatile LONG64* p = (volatile LONG64*)ptr;
#ifdef _WIN64
    return (int64_t)InterlockedCompareExchange64(p, 0, 0);
#else
    return *(const volatile int64_t*)p; // Not truly atomic on 32-bit
#endif
}

__declspec(dllexport) int32_t __atomic_load_4(const void* ptr, int memorder) {
    (void)memorder;
    volatile LONG* p = (volatile LONG*)ptr;
    return (int32_t)InterlockedCompareExchange(p, 0, 0);
}

// Generic __atomic_load that dispatches based on size
__declspec(dllexport) void __atomic_load(int size, const void* ptr, void* ret, int memorder) {
    if (size == 4) {
        *(int32_t*)ret = __atomic_load_4(ptr, memorder);
    } else if (size == 8) {
        *(int64_t*)ret = __atomic_load_8(ptr, memorder);
    } else {
        // Fallback: plain memcpy-like load for other sizes
        const volatile char* src = (const volatile char*)ptr;
        char* dst = (char*)ret;
        for (int i = 0; i < size; ++i) dst[i] = src[i];
    }
}

// __atomic_store_8 and __atomic_store_4 store val into *ptr
__declspec(dllexport) void __atomic_store_8(void* ptr, int64_t val, int memorder) {
    (void)memorder;
    volatile LONG64* p = (volatile LONG64*)ptr;
#ifdef _WIN64
    InterlockedExchange64(p, val);
#else
    *(volatile int64_t*)p = val; // Not truly atomic on 32-bit
#endif
}

__declspec(dllexport) void __atomic_store_4(void* ptr, int32_t val, int memorder) {
    (void)memorder;
    volatile LONG* p = (volatile LONG*)ptr;
    InterlockedExchange(p, val);
}

// Generic __atomic_store that dispatches based on size
__declspec(dllexport) void __atomic_store(int size, void* ptr, const void* val, int memorder) {
    if (size == 4) {
        __atomic_store_4(ptr, *(const int32_t*)val, memorder);
    } else if (size == 8) {
        __atomic_store_8(ptr, *(const int64_t*)val, memorder);
    } else {
        // Fallback: plain memcpy-like store for other sizes
        volatile char* dst = (volatile char*)ptr;
        const char* src = (const char*)val;
        for (int i = 0; i < size; ++i) dst[i] = src[i];
    }
}
