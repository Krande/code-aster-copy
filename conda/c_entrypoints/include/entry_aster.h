#ifndef ASTER_ENTRY_ASTER_H_
#define ASTER_ENTRY_ASTER_H_

#include "Windows.h"
#include "astercxx.h"

#include "aster_init.h"
#include "aster_numpy.h"
#include "aster_pybind.h"


#ifdef __cplusplus
extern "C" {
#endif
// export to dll extern PyObject *PyInit_aster( void );
extern PyObject *PyInit_aster( void );

} // extern "C"

#endif  // ASTER_ENTRY_ASTER_H_