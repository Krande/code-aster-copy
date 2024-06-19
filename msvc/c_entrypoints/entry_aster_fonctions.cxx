#include "entry_helpers.h"

PyInit_func_t original_PyInit_aster_fonctions = NULL;

extern "C" PyObject* PyInit_aster_fonctions()
{
    if (!original_PyInit_aster_fonctions)
    {
        LoadDllAndGetFunction("bibc.dll", "PyInit_aster_fonctions", original_PyInit_aster_fonctions);
    }

    if (original_PyInit_aster_fonctions)
    {
        std::cout << "Calling original PyInit_aster_fonctions" << std::endl;
        return original_PyInit_aster_fonctions();
    }
    std::cout << "Returning NULL" << std::endl;
    return NULL;
}
