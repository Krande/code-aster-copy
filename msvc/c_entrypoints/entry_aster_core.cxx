#include "entry_helpers.h"

PyInit_func_t original_PyInit_aster_core = NULL;

extern "C" PyObject* PyInit_aster_core()
{
    if (!original_PyInit_aster_core)
    {
        LoadDllAndGetFunction("bibc.dll", "PyInit_aster_core", original_PyInit_aster_core);
    }

    if (original_PyInit_aster_core)
    {
        std::cout << "Calling original PyInit_aster_core" << std::endl;
        return original_PyInit_aster_core();
    }
    std::cout << "Returning NULL" << std::endl;
    return NULL;
}
