#include "entry_helpers.h"

PyInit_func_t original_PyInit_libaster = NULL;

extern "C" PyObject* PyInit_libaster()
{
    if (!original_PyInit_libaster)
    {
        LoadDllAndGetFunction("bibcxx.dll", "PyInit_libaster", original_PyInit_libaster);
    }

    if (original_PyInit_libaster)
    {
        std::cout << "Calling original PyInit_libaster" << std::endl;
        return original_PyInit_libaster();
    }
    std::cout << "Returning NULL" << std::endl;
    return NULL;
}
