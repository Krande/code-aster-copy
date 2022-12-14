/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2017 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

/* person_in_charge: nicolas.sellenet at edf.fr */

#include "Python.h"
#include "aster.h"
#include "aster_fort_utils.h"
#include "aster_utils.h"
#include "definition_pt.h"

#include "dll_register.h"
#include "dll_mfront.h"

#ifdef ASTER_HAVE_MFRONT
#include "MFrontBehaviour.h"
#endif

#ifdef ASTER_PLATFORM_POSIX
#include <dlfcn.h>
#endif
PyObject* get_dll_register_dict();


/* *********************************************************************
 *
 *                          MFRONT interface
 *
 * *********************************************************************/

void DEFSSSSP(MFRONT_SET_DOUBLE_PARAMETER, mfront_set_double_parameter,
    char* nomlib, STRING_SIZE lnomlib, char* nomsub, STRING_SIZE lnomsub,
    char* nommod, STRING_SIZE lnommod,
    char* nomparam, STRING_SIZE lnomparam, ASTERDOUBLE* value)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper : wrapper to the MFRONT set function through the function pointer
     * Load the library if necessary (at the first call).
    */
    char *libname, *symbol, *model, *symbname=NULL, *nom_param;
    FUNC_MFRONT_SET_DOUBLE(f_mfront) = NULL;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    libname = MakeCStrFromFStr(nomlib, lnomlib);
    symbol = MakeCStrFromFStr(nomsub, lnomsub);
    model = MakeCStrFromFStr(nommod, lnommod);
    nom_param = MakeCStrFromFStr(nomparam, lnomparam);

    mfront_name(libname, symbol, model, "_setParameter", &symbname);
    if ( symbname == NULL ) return;

    f_mfront = (FUNC_MFRONT_SET_DOUBLE())libsymb_get_symbol(DLL_DICT, libname, symbname);
    CALLMFRONTSETDOUBLE(f_mfront, nom_param, *value);
    FreeStr(libname);
    FreeStr(symbol);
    FreeStr(model);
    FreeStr(nom_param);
    FreeStr(symbname);
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

void DEFSSSSP(MFRONT_SET_INTEGER_PARAMETER, mfront_set_integer_parameter,
    char* nomlib, STRING_SIZE lnomlib, char* nomsub, STRING_SIZE lnomsub,
    char* nommod, STRING_SIZE lnommod,
    char* nomparam, STRING_SIZE lnomparam, ASTERINTEGER* value)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper : wrapper to the MFRONT set function through the function pointer
     * Load the library if necessary (at the first call).
    */
    char *libname, *symbol, *model, *symbname=NULL, *nom_param;
    FUNC_MFRONT_SET_INTEGER(f_mfront) = NULL;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    libname = MakeCStrFromFStr(nomlib, lnomlib);
    symbol = MakeCStrFromFStr(nomsub, lnomsub);
    model = MakeCStrFromFStr(nommod, lnommod);
    nom_param = MakeCStrFromFStr(nomparam, lnomparam);

    mfront_name(libname, symbol, model, "_setUnsignedShortParameter", &symbname);
    if ( symbname == NULL ) return;

    f_mfront = (FUNC_MFRONT_SET_INTEGER())libsymb_get_symbol(DLL_DICT, libname, symbname);
    CALLMFRONTSETINTEGER(f_mfront, nom_param, (unsigned short)(*value));
    FreeStr(libname);
    FreeStr(symbol);
    FreeStr(model);
    FreeStr(nom_param);
    FreeStr(symbname);
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

void DEFPPSP(MFRONT_GET_EXTERNAL_STATE_VARIABLE,
             mfront_get_external_state_variable,
             ASTERINTEGER* pliesv, ASTERINTEGER* pnbesv,
             char* txval, STRING_SIZE ltx, ASTERINTEGER* nbvarc)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper
    */

    char** ext_var = (char**)*pliesv;

    unsigned short* nb_ext_var = (unsigned short*)*pnbesv;
    AS_ASSERT(*nb_ext_var <= ltx);
    *nbvarc = *nb_ext_var;

    unsigned short i;
    for ( i = 0; i < *nb_ext_var; ++i )
    {
        SetTabFStr( txval, i, ext_var[i], 8 );
    }
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

int mfront_get_number_of_internal_state_variables(char* nomlib, STRING_SIZE lnomlib,
                                                  char* nomsub, STRING_SIZE lnomsub,
                                                  char* nommod, STRING_SIZE lnommod)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper
    */
    char *libname, *symbol, *model, *symbname=NULL;
    char * name1;
    int retour = 0;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    libname = MakeCStrFromFStr(nomlib, lnomlib);
    symbol = MakeCStrFromFStr(nomsub, lnomsub);
    model = MakeCStrFromFStr(nommod, lnommod);

    mfront_name(libname, symbol, model, "_nInternalStateVariables", &symbname);
    if ( symbname == NULL )
    {
        name1 = (char *)malloc(strlen(symbol)\
                               + strlen("_nInternalStateVariables") + 1);
        strcpy(name1, symbol);
        strcat(name1, "_nInternalStateVariables");
        error_symbol_not_found(libname, name1);
    }

    unsigned short* nbvari2 = (unsigned short*)libsymb_get_symbol(DLL_DICT, libname, symbname);

    FreeStr(libname);
    FreeStr(symbol);
    FreeStr(symbname);
    FreeStr(model);
    return *nbvari2;
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}


int mfront_get_strain_model(char* nomlib, STRING_SIZE lnomlib,
                            char* nomsub, STRING_SIZE lnomsub,
                            char* nommod, STRING_SIZE lnommod)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper
    */
    char *libname, *symbol, *model, *symbname=NULL;
    char * name1;
    int retour = 0;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    libname = MakeCStrFromFStr(nomlib, lnomlib);
    symbol = MakeCStrFromFStr(nomsub, lnomsub);
    model = MakeCStrFromFStr(nommod, lnommod);

    mfront_name(libname, symbol, model, "_FiniteStrainFormulation", &symbname);
    if ( symbname == NULL )
    {
        FreeStr(libname);
        FreeStr(symbol);
        FreeStr(symbname);
        FreeStr(model);
        return 0;
    }
    else
    {
        unsigned short* type = (unsigned short*)libsymb_get_symbol(DLL_DICT, libname, symbname);
        FreeStr(libname);
        FreeStr(symbol);
        FreeStr(symbname);
        FreeStr(model);
        return *type;
    }
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

void DEFSSSP(MFRONT_GET_NUMBER_OF_INTERNAL_STATE_VARIABLES,
             mfront_get_number_of_internal_state_variables,
             char* nomlib, STRING_SIZE lnomlib,
             char* nomsub, STRING_SIZE lnomsub,
             char* nommod, STRING_SIZE lnommod,
             ASTERINTEGER* nbvari)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper
    */
    int nbvari2 = mfront_get_number_of_internal_state_variables(nomlib, lnomlib,
                                                                nomsub, lnomsub,
                                                                nommod, lnommod);
    *nbvari = nbvari2;
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

void DEFSSSP(MFRONT_GET_STRAIN_MODEL,
             mfront_get_strain_model,
             char* nomlib, STRING_SIZE lnomlib,
             char* nomsub, STRING_SIZE lnomsub,
             char* nommod, STRING_SIZE lnommod,
             ASTERINTEGER* strain_model)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper
    */
    int strain_model2 = mfront_get_strain_model(nomlib, lnomlib,
                                                nomsub, lnomsub,
                                                nommod, lnommod);
    *strain_model = strain_model2;
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

void DEFSSSS(MFRONT_GET_INTERNAL_STATE_VARIABLES_TYPES,
              mfront_get_internal_state_variables_types,
              char* nomlib, STRING_SIZE lnomlib,
              char* nomsub, STRING_SIZE lnomsub,
              char* nommod, STRING_SIZE lnommod,
              char* txval, STRING_SIZE ltx)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper
    */
    char *libname, *symbol, *model, *symbname=NULL;
    char * name1;
    int retour = 0;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    libname = MakeCStrFromFStr(nomlib, lnomlib);
    symbol = MakeCStrFromFStr(nomsub, lnomsub);
    model = MakeCStrFromFStr(nommod, lnommod);

    mfront_name(libname, symbol, model, "_nInternalStateVariables", &symbname);
    if ( symbname == NULL ) {
        name1 = (char *)malloc(strlen(symbol) + strlen(model)
                               + strlen("_nInternalStateVariables") + 1);
        strcpy(name1, symbol);
        strcat(name1, model);
        strcat(name1, "_nInternalStateVariables");
        error_symbol_not_found(libname, name1);
    }
    unsigned short* nb_int_var = (unsigned short*)libsymb_get_symbol(DLL_DICT, libname, symbname);

    mfront_name(libname, symbol, model, "_InternalStateVariablesTypes", &symbname);
    if ( symbname == NULL ) {
        name1 = (char *)malloc(strlen(symbol) + strlen(model)
                               + strlen("_InternalStateVariablesTypes") + 1);
        strcpy(name1, symbol);
        strcat(name1, model);
        strcat(name1, "_InternalStateVariablesTypes");
        error_symbol_not_found(libname, name1);
    }

    int* int_var = (int*)libsymb_get_symbol(DLL_DICT, libname, symbname);

    unsigned short i;
    for ( i = 0; i < *nb_int_var; ++i )
    {
        if ( int_var[i] == 0 )
        {
            SetTabFStr( txval, i, "scalar", ltx );
        }
        else if ( int_var[i] == 1 )
        {
            SetTabFStr( txval, i, "vector", ltx );
        }
        else if ( int_var[i] == 3 )
        {
            SetTabFStr( txval, i, "tensor", ltx );
        }
        else
        {
            AS_ASSERT( int_var[i] == 0 || int_var[i] == 1 || int_var[i] == 3);
        }
    }

    FreeStr(libname);
    FreeStr(symbol);
    FreeStr(symbname);
    FreeStr(model);
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

void DEFSSSSP(MFRONT_GET_INTERNAL_STATE_VARIABLES,
              mfront_get_internal_state_variables,
              char* nomlib, STRING_SIZE lnomlib,
              char* nomsub, STRING_SIZE lnomsub,
              char* nommod, STRING_SIZE lnommod,
              char* txval, STRING_SIZE ltx, ASTERINTEGER* nbintvar)
{
#ifdef ASTER_PLATFORM_POSIX
#define MIN(A,B)  ((A) < (B) ? (A) : (B))

    /* MFRONT Wrapper
    */
    char *libname, *symbol, *model, *symbname=NULL;
    char * name1;
    int retour = 0;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    int nbcheck = mfront_get_number_of_internal_state_variables(nomlib, lnomlib,
                                                                nomsub, lnomsub,
                                                                nommod, lnommod);
    nbcheck = MIN(nbcheck, *nbintvar);

    libname = MakeCStrFromFStr(nomlib, lnomlib);
    symbol = MakeCStrFromFStr(nomsub, lnomsub);
    model = MakeCStrFromFStr(nommod, lnommod);

    mfront_name(libname, symbol, model, "_InternalStateVariables", &symbname);
    if ( symbname == NULL )
    {
        name1 = (char *)malloc(strlen(symbol)\
                               + strlen("_InternalStateVariables") + 1);
        strcpy(name1, symbol);
        strcat(name1, "_InternalStateVariables");
        error_symbol_not_found(libname, name1);
    }

    char** int_var = (char**)libsymb_get_symbol(DLL_DICT, libname, symbname);

    unsigned short i;
    for ( i = 0; i < nbcheck; ++i )
    {
        AS_ASSERT(strlen(int_var[i]) <= ltx);
        SetTabFStr( txval, i, int_var[i], ltx );
    }

    FreeStr(libname);
    FreeStr(symbol);
    FreeStr(symbname);
    FreeStr(model);
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

void DEFSSSPPPPP(MFRONT_GET_POINTERS,
                 mfront_get_pointers,
                 char* nomlib, STRING_SIZE lnomlib,
                 char* nomsub, STRING_SIZE lnomsub,
                 char* nommod, STRING_SIZE lnommod,
                 ASTERINTEGER* pliesv, ASTERINTEGER* pnbesv, ASTERINTEGER* pfcmfr,
                 ASTERINTEGER* pmatprop, ASTERINTEGER* pnbprop)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper
    */
    char *libname, *symbol, *model, *symbname=NULL;
    char * name1;
    int retour = 0;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    libname = MakeCStrFromFStr(nomlib, lnomlib);
    symbol = MakeCStrFromFStr(nomsub, lnomsub);
    model = MakeCStrFromFStr(nommod, lnommod);

    if ( ! libsymb_is_known(DLL_DICT, libname, symbol) ) {
        retour = load_mfront_lib(libname, symbol);
        if (retour == 1)
        {
            error_symbol_not_found(libname, symbname);
        }
    }
    *pfcmfr = (ASTERINTEGER)libsymb_get_symbol(DLL_DICT, libname, symbol);

    mfront_name(libname, symbol, model, "_ExternalStateVariables", &symbname);
    if ( symbname == NULL )
    {
        name1 = (char *)malloc(strlen(symbol) + strlen(model)\
                               + strlen("_ExternalStateVariables") + 1);
        strcpy(name1, symbol);
        strcat(name1, model);
        strcat(name1, "_ExternalStateVariables");
        error_symbol_not_found(libname, name1);
    }

//     char** test_char = libsymb_get_symbol(DLL_DICT, libname, symbname);
    *pliesv = (ASTERINTEGER)libsymb_get_symbol(DLL_DICT, libname, symbname);

    mfront_name(libname, symbol, model, "_nExternalStateVariables", &symbname);
    if ( symbname == NULL ) {
        name1 = (char *)malloc(strlen(symbol) + strlen(model)\
                               + strlen("_ExternalStateVariables") + 1);
        strcpy(name1, symbol);
        strcat(name1, model);
        strcat(name1, "_ExternalStateVariables");
        error_symbol_not_found(libname, name1);
    }
//     int* test_int = libsymb_get_symbol(DLL_DICT, libname, symbname);
    *pnbesv = (ASTERINTEGER)libsymb_get_symbol(DLL_DICT, libname, symbname);
    if ( symbname == NULL ) {
        error_symbol_not_found(libname, symbname);
    }

    // may be used for performance reason: pointers in a cache
    *pmatprop = (ASTERINTEGER)0;
    *pnbprop = (ASTERINTEGER)0;

    FreeStr(libname);
    FreeStr(model);
    FreeStr(symbol);
    FreeStr(symbname);
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

void DEFSSSP(MFRONT_SET_OUTOFBOUNDS_POLICY,
             mfront_set_outofbounds_policy,
    char* nomlib, STRING_SIZE lnomlib, char* nomsub, STRING_SIZE lnomsub,
    char* nommod, STRING_SIZE lnommod, ASTERINTEGER* value)
{
#ifdef ASTER_PLATFORM_POSIX
    char *libname, *symbol, *model, *symbname=NULL, *name1;
    int retour = 0;
    FUNC_MFRONT_SET_OUTOFBOUNDS_POLICY(f_mfront) = NULL;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    libname = MakeCStrFromFStr(nomlib, lnomlib);
    symbol = MakeCStrFromFStr(nomsub, lnomsub);
    model = MakeCStrFromFStr(nommod, lnommod);

    mfront_name(libname, symbol, model, "_setOutOfBoundsPolicy", &symbname);
    if ( symbname == NULL )
    {
        name1 = (char *)malloc(strlen(symbol) + strlen(model)\
                               + strlen("_setOutOfBoundsPolicy") + 1);
        strcpy(name1, symbol);
        strcat(name1, model);
        strcat(name1, "_setOutOfBoundsPolicy");
        error_symbol_not_found(libname, name1);
    }

    f_mfront = (FUNC_MFRONT_SET_OUTOFBOUNDS_POLICY())libsymb_get_symbol(DLL_DICT,libname,symbname);
    CALLMFRONTSETOUTOFBOUNDSPOLICY(f_mfront, (unsigned short)(*value));

    FreeStr(libname);
    FreeStr(model);
    FreeStr(symbol);
    FreeStr(symbname);
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

void DEFSSSPP(MFRONT_GET_NBVARI, mfront_get_nbvari,
    char* nomlib, STRING_SIZE lnomlib, char* nomsub, STRING_SIZE lnomsub,
    char* nommod, STRING_SIZE lnommod,
    ASTERINTEGER* ndim, ASTERINTEGER* nbvari)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper
    */
    char *libname, *symbol, *model, *symbname=NULL, *name1;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    libname = MakeCStrFromFStr(nomlib, lnomlib);
    symbol = MakeCStrFromFStr(nomsub, lnomsub);
    model = MakeCStrFromFStr(nommod, lnommod);

    mfront_name(libname, symbol, model, "_InternalStateVariablesTypes", &symbname);
    if ( symbname == NULL ) {
        name1 = (char *)malloc(strlen(symbol) + strlen(model)
                               + strlen("_InternalStateVariablesTypes") + 1);
        strcpy(name1, symbol);
        strcat(name1, model);
        strcat(name1, "_InternalStateVariablesTypes");
        error_symbol_not_found(libname, name1);
    }

    int* int_var = (int*)libsymb_get_symbol(DLL_DICT, libname, symbname);

    mfront_name(libname, symbol, model, "_nInternalStateVariables", &symbname);
    if ( symbname == NULL ) {
        name1 = (char *)malloc(strlen(symbol) + strlen(model)
                               + strlen("_nInternalStateVariables") + 1);
        strcpy(name1, symbol);
        strcat(name1, model);
        strcat(name1, "_nInternalStateVariables");
        error_symbol_not_found(libname, name1);
    }
    unsigned short* nb_int_var = (unsigned short*)libsymb_get_symbol(DLL_DICT, libname, symbname);

    *nbvari = 0;
    unsigned short i;
    for ( i = 0; i < *nb_int_var; ++i )
    {
        if ( int_var[i] == 0 )
        {
            ++(*nbvari);
        }
        else if ( int_var[i] == 1 )
        {
            if ( *ndim == 2 )
            {
                (*nbvari) += 4;
            }
            else if ( *ndim == 3 )
            {
                (*nbvari) += 6;
            }
            else
            {
                AS_ASSERT( *ndim == 2 || *ndim == 3 );
            }
        }
        else if ( int_var[i] == 3 )
        {
                 (*nbvari) += 9;
        }
        else
        {
            AS_ASSERT( int_var[i] == 0 || int_var[i] == 1 || int_var[i] == 3);
        }
    }

    FreeStr(libname);
    FreeStr(symbol);
    FreeStr(model);
    FreeStr(symbname);
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

void DEFPPPPPPPPPPPPPPPPPP(MFRONT_BEHAVIOUR, mfront_behaviour,
    ASTERINTEGER* pfcmfr, ASTERDOUBLE* stress, ASTERDOUBLE* statev,
    ASTERDOUBLE* ddsdde, ASTERDOUBLE* stran, ASTERDOUBLE* dstran,
    ASTERDOUBLE* dtime, ASTERDOUBLE* temp, ASTERDOUBLE* dtemp,
    ASTERDOUBLE* predef, ASTERDOUBLE* dpred, ASTERINTEGER* ntens,
    ASTERINTEGER* nstatv, ASTERDOUBLE* props, ASTERINTEGER* nprops,
    ASTERDOUBLE* drot, ASTERDOUBLE* pnewdt, ASTERINTEGER* nummod)
{
#ifdef ASTER_PLATFORM_POSIX
    /* MFRONT Wrapper : wrapper to the MFRONT function through the function pointer
     * Load the library if necessary (at the first call).
    */
    FUNC_MFRONT(f_mfront) = NULL;

    f_mfront = (FUNC_MFRONT())(*pfcmfr);

    CALLMFRONTBEHAVIOUR(*f_mfront,
        stress, statev, ddsdde, stran, dstran,
        dtime, temp, dtemp, predef, dpred,
        ntens, nstatv, props, nprops, drot,
        pnewdt, nummod);
#else
    printf("Not available under Windows.\n");
    abort();
#endif
}

/**
 * \brief Fill the array of the material properties names
 * @param pmatprop  Pointer on the data in the library
 * @param nbval     Number of values of material properties
 * @param txval     Array of strings
 */
void DEFSPS(MFRONT_GET_MATER_PROP,
            mfront_get_mater_prop,
             _IN char* rela, STRING_SIZE lrela,
            _OUT ASTERINTEGER* nbval,
            _OUT char* txval, STRING_SIZE ltx)
{
#ifdef ASTER_HAVE_MFRONT
    /* MFRONT Wrapper
    */
    char *crela;
    char library[] = "lib" ASTER_BEHAVIOUR_LIB;
    unsigned int i, size;
    char **props;
    AS_ASSERT(ltx == 16);

    crela = MakeCStrFromFStr(rela, lrela);
    props = getTridimMaterialPropertiesNames(crela, &size);
    for (i = 0; i < size; ++i) {
        SetTabFStr( txval, i, props[i], 16 );
        free(props[i]);
    }
    *nbval = (ASTERINTEGER)size;
    free(props);
    FreeStr(crela);
#else
    printf("MFront library is required for this functionnality.\n");
    abort();
#endif
}

int load_mfront_lib(const char* libname, const char* symbol)
{
    void *mfront_handle;
    char *error;
    char symbol_[256], *valk;
    ASTERINTEGER ibid=0, n0=0, nk=0;
    ASTERDOUBLE rbid=0.;
    FUNC_MFRONT(f_mfront) = NULL;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    AS_ASSERT(strlen(symbol) < 255);
    strcpy(symbol_, symbol);

    DEBUG_DLL_VV("Loading '%s'%s ", libname, "...");
    mfront_handle = dlopen(libname, RTLD_NOW);
    if ( ! mfront_handle ) {
        printf("\n%s\n", dlerror());
        nk = 2;
        valk = MakeTabFStr(nk, VALK_SIZE);
        SetTabFStr(valk, 0, "MFRONT", VALK_SIZE);
        SetTabFStr(valk, 1, (char *)libname, VALK_SIZE);
        CALL_UTMESS_CORE("F", "FERMETUR_13", &nk, valk, &n0,
                         &ibid, &n0, &rbid, &n0, " ");
        FreeStr(valk);  // uncallable
    }
    DEBUG_DLL_VV("searching symbol '%s'%s ", symbol, "...");
    dlerror();    /* Clear any existing error */

    *(void **) (&f_mfront) = dlsym(mfront_handle, symbol);
    if ((error = dlerror()) != NULL)  {
        dlerror();
        strcat(symbol_, "_");
        DEBUG_DLL_VV("trying symbol '%s'%s ", symbol_, "...");
        *(void **) (&f_mfront) = dlsym(mfront_handle, symbol_);
    }

    if ((error = dlerror()) != NULL)  {
        DEBUG_DLL_VV("not found %s%s\n", ":-(", "");
        return 1;
    }
    DEBUG_DLL_VV("found: %s %p", "address", (char *)f_mfront);

    /* register these MFRONT lib */
    if ( libsymb_register(DLL_DICT, libname, symbol,
                            mfront_handle, (FUNC_PTR)f_mfront) ) {
        printf("Registering of '%s' and '%s' failed!\n", libname, symbol);
    }
    return 0;
}

char* test_mfront_symbol(const char* libname, char* name1, char* name2)
{
    int retour = 0;
    PyObject* DLL_DICT;
    DLL_DICT = get_dll_register_dict();

    if ( ! libsymb_is_known(DLL_DICT, libname, name1) )
    {
        retour = load_mfront_lib(libname, name1);
        if ( retour == 0 )
            return name1;
    }
    else
        return name1;
    if ( ! libsymb_is_known(DLL_DICT, libname, name2) )
    {
        retour = load_mfront_lib(libname, name2);
        if ( retour == 0 )
            return name2;
    }
    else
        return name2;
    return NULL;
}

void mfront_name(
         _IN char* libname, _IN char* symbol, _IN char* model,
         _IN char* basename, _OUT char** name)
{
    char *name1, *name2;

    name1 = (char *)malloc(strlen(symbol) + strlen(model) + strlen(basename) + 1);
    strcpy(name1, symbol);
    strcat(name1, model);
    strcat(name1, basename);

    name2 = (char *)malloc(strlen(symbol) + strlen(basename) + 1);
    strcpy(name2, symbol);
    strcat(name2, basename);
    DEBUG_DLL_VV("name1: '%s' name2: '%s'", name1, name2);

    *name = test_mfront_symbol(libname, name1, name2);
    if ( *name == NULL ) {
        DEBUG_DLL_VV(" libname = >%s<%s", libname, " ")
        DEBUG_DLL_VV(" symbol1 = >%s<, symbol2 = >%s<", name1, name2)
        free(name1);
        free(name2);
    }
    else if ( strcmp(*name, name1) == 0 ) {
        free(name2);
    }
    else if ( strcmp(*name, name2) == 0 ) {
        free(name1);
    }
}

void error_symbol_not_found(const char* libname, const char* symbname)
{
    char *valk;
    ASTERINTEGER ibid=0, n0=0, nk=3;
    ASTERDOUBLE rbid=0.;
    valk = MakeTabFStr(nk, VALK_SIZE);
    SetTabFStr(valk, 0, "MFRONT", VALK_SIZE);
    SetTabFStr(valk, 1, (char *)libname, VALK_SIZE);
    SetTabFStr(valk, 2, (char *)symbname, VALK_SIZE);
    CALL_UTMESS_CORE("F", "FERMETUR_14", &nk, valk, &n0,
                     &ibid, &n0, &rbid, &n0, " ");
    FreeStr(valk);  // uncallable
}
