#include "entry_helpers.h"

PyInit_func_t original_PyInit_aster = NULL;

extern "C" PyObject* PyInit_aster()
{
    if (!original_PyInit_aster)
    {
        LoadDllAndGetFunction("bibc.dll", "PyInit_aster", original_PyInit_aster);
    }

    if (original_PyInit_aster)
    {
    #ifdef _DEBUG
        std::cout << "Calling original PyInit_aster" << std::endl;
    #endif
        return original_PyInit_aster();
    }
    std::cout << "Returning NULL" << std::endl;
    return NULL;
}
