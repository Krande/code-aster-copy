cmake_minimum_required(VERSION 3.28)

project(bibcxx_dummy LANGUAGES CXX C Fortran)

# I use these to generate the .def file(s) for the shared libraries
OPTION(BUILD_BIBCXX "Build the bibcxx library" OFF)
OPTION(BUILD_BIBC "Build the bibc library" OFF)
OPTION(BUILD_BIBFOR "Build the bibfor library" OFF)
OPTION(BUILD_ASTER "Build the aster library" OFF)
OPTION(BUILD_MFRONT "Build the mfront library" OFF)

set(INCL_DIRS
        bibcxx/include
        bibcxx bibc/include
        bibfor/include
        build/std/debug
        $ENV{CONDA_PREFIX}/include
        $ENV{CONDA_PREFIX}/Library/include
        $ENV{CONDA_PREFIX}/Lib/site-packages/numpy/core/include build/std/debug
)

set(DEP_PATHS
        /LIBPATH:$ENV{CONDA_PREFIX}/libs /LIBPATH:$ENV{CONDA_PREFIX}/Library/lib /LIBPATH:$ENV{CONDA_PREFIX}/Library/bin
)
# create a list
set(DEP_LIBS
        pthread.lib medfwrap.lib hdf5.lib metis.lib MFrontGenericInterface.lib
        scotch.lib scotcherr.lib
        mkl_intel_ilp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib libiomp5md.lib
)

# BIBCXX
if (BUILD_BIBCXX)
    # Find all .cxx files recursively in the subdirectory bibcxx
    file(GLOB_RECURSE bibcxx_SOURCES "bibcxx/*.cxx")
    list(REMOVE_ITEM bibcxx_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bibcxx/MFront/MFrontBehaviour.cxx")

    # Create a shared library from the sources
    add_library(bibcxx_shared SHARED ${bibcxx_SOURCES})

    # Add /std:c++17 to the compiler flags
    target_compile_options(bibcxx_shared PRIVATE /std:c++17 /nologo /FS /EHs /permissive- /Zi /DTFEL_ARCH64 /DASTER_HAVE_INTEL_IFORT /DH5_BUILT_AS_DYNAMIC_LIB)

    # Add lib path to the linker flags
    target_link_options(bibcxx_shared PRIVATE ${DEP_PATHS})

    # add nodefaultlib:libcmt.lib to the linker flags
    #target_link_options(bibcxx_shared PRIVATE /nodefaultlib)

    # Add .lib files to the linker flags
    target_link_libraries(bibcxx_shared PRIVATE ${DEP_LIBS})

    # Include the headers located at bibcxx/include
    target_include_directories(bibcxx_shared PUBLIC ${INCL_DIRS})

    # use the WINDOWS_EXPORT_ALL_SYMBOLS to create a .def file
    set_target_properties(bibcxx_shared PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS TRUE
    )
endif (BUILD_BIBCXX)

# MFRONT
if (BUILD_MFRONT)
    # Find all .cxx files recursively in the subdirectory mfront
    file(GLOB_RECURSE mfront_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bibcxx/MFront/*.cxx")

    # Create a shared library from the sources
    add_library(mfront_shared SHARED ${mfront_SOURCES})

    # Add /std:c++17 to the compiler flags
    target_compile_options(mfront_shared PRIVATE /std:c++17 /nologo /DTFEL_ARCH64)

    # Add lib path to the linker flags
    target_link_options(mfront_shared PRIVATE ${DEP_PATHS})

    # Add .lib files to the linker flags
    target_link_libraries(mfront_shared PRIVATE ${DEP_LIBS})

    # Include the headers located at mfront/include
    target_include_directories(mfront_shared PUBLIC ${INCL_DIRS})

    # use the WINDOWS_EXPORT_ALL_SYMBOLS to create a .def file
    set_target_properties(mfront_shared PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS TRUE
    )
endif (BUILD_MFRONT)

# BIBC
if (BUILD_BIBC)
    # Find all .c files recursively in the subdirectory bibc
    file(GLOB_RECURSE bibc_SOURCES "bibc/*.c")
    # remove file name python.c
    list(REMOVE_ITEM bibc_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bibc/supervis/python.c")

    # Create a shared library from the sources
    add_library(bibc_shared SHARED ${bibc_SOURCES})

    # Add /std:c++17 to the compiler flags
    target_compile_options(bibc_shared PRIVATE /FS /Zi /0d /nologo /DTFEL_ARCH64 /DASTER_HAVE_INTEL_IFORT /DH5_BUILT_AS_DYNAMIC_LIB)

    # Add lib path to the linker flags
    target_link_options(bibc_shared PRIVATE ${DEP_PATHS})

    # Add .lib files to the linker flags
    target_link_libraries(bibc_shared PRIVATE ${DEP_LIBS})

    # Include the headers located at bibc/include
    target_include_directories(bibc_shared PUBLIC ${INCL_DIRS})

    # use the WINDOWS_EXPORT_ALL_SYMBOLS to create a .def file
    set_target_properties(bibc_shared PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS TRUE
            #        POSITION_INDEPENDENT_CODE ON
    )
endif (BUILD_BIBC)
# Bibfor
if (BUILD_BIBFOR)
    # Find all .F90 files recursively in the subdirectory bibfor
    file(GLOB_RECURSE bibfor_SOURCES "bibfor/*.F90")

    # Create a shared library from the sources
    add_library(bibfor_shared SHARED ${bibfor_SOURCES})

    # Add fortran flags /names:lowercase /assume:underscore to the compiler flags
    target_compile_options(bibfor_shared PRIVATE /nologo /DTFEL_ARCH64 /names:lowercase /assume:underscore)

    # Add lib path to the linker flags
    target_link_options(bibfor_shared PRIVATE ${DEP_PATHS})

    # Add .lib files to the linker flags
    target_link_libraries(bibfor_shared PRIVATE ${DEP_LIBS})

    # Include the headers located at bibfor/include
    target_include_directories(bibfor_shared PUBLIC ${INCL_DIRS})

    # use the WINDOWS_EXPORT_ALL_SYMBOLS to create a .def file
    set_target_properties(bibfor_shared PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS TRUE
    )
endif (BUILD_BIBFOR)

# Aster
if (BUILD_ASTER)
    set(aster_SOURCES
            "${CMAKE_CURRENT_SOURCE_DIR}/bibc/supervis/python.c"
    )

    # Create a shared library from the sources
    add_library(aster_shared SHARED ${aster_SOURCES})

    # Add /std:c++17 to the compiler flags
    target_compile_options(aster_shared PRIVATE /nologo /DTFEL_ARCH64)

    target_include_directories(aster_shared PUBLIC ${INCL_DIRS})

    # Add lib path to the linker flags
    target_link_options(aster_shared PRIVATE ${DEP_PATHS})

    target_link_libraries(aster_shared PRIVATE ${DEP_LIBS} bibcxx_shared bibc_shared)

    set_target_properties(
            aster_shared PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS TRUE
    )
endif (BUILD_ASTER)