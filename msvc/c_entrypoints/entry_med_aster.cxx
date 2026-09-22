#include "entry_helpers.h"

PyInit_func_t original_PyInit_med_aster = NULL;

extern "C" PyObject* PyInit_med_aster()
{
    if (!original_PyInit_med_aster)
    {
        LoadDllAndGetFunction("bibc.dll", "PyInit_med_aster", original_PyInit_med_aster);
    }

    if (original_PyInit_med_aster)
    {
        return original_PyInit_med_aster();
    }
    std::cout << "Returning NULL" << std::endl;
    return NULL;
}
