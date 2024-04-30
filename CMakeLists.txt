project(bibcxx_dummy LANGUAGES CXX C)

cmake_minimum_required(VERSION 3.28)

# BIBCXX
# Find all .cxx files recursively in the subdirectory bibcxx
#file(GLOB_RECURSE bibcxx_SOURCES "bibcxx/*.cxx")
#
## Create a shared library from the sources
#add_library(bibcxx_shared SHARED ${bibcxx_SOURCES})
#
## Add /std:c++17 to the compiler flags
#target_compile_options(bibcxx_shared PRIVATE /std:c++17 /nologo /MD /DTFEL_ARCH64)
#
## Add lib path to the linker flags
#target_link_options(bibcxx_shared PRIVATE /LIBPATH:$ENV{CONDA_PREFIX}/libs)
#
## Include the headers located at bibcxx/include
#target_include_directories(bibcxx_shared PUBLIC bibcxx/include bibcxx bibc/include bibfor/include $ENV{CONDA_PREFIX}/include $ENV{CONDA_PREFIX}/Library/include extern/include $ENV{CONDA_PREFIX}/Lib/site-packages/numpy/core/include build/std/debug)
#
## use the WINDOWS_EXPORT_ALL_SYMBOLS to create a .def file
#set_target_properties(bibcxx_shared PROPERTIES
#        WINDOWS_EXPORT_ALL_SYMBOLS TRUE
#        POSITION_INDEPENDENT_CODE ON
#)

# BIBC
# Find all .c files recursively in the subdirectory bibc
file(GLOB_RECURSE bibc_SOURCES "bibc/*.c")
# remove file name python.c
list(REMOVE_ITEM bibc_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/bibc/supervis/python.c")

# Create a shared library from the sources
add_library(bibc_shared SHARED ${bibc_SOURCES})

# Add /std:c++17 to the compiler flags
target_compile_options(bibc_shared PRIVATE /nologo /MD /DTFEL_ARCH64)

# Add lib path to the linker flags
target_link_options(bibc_shared PRIVATE /LIBPATH:$ENV{CONDA_PREFIX}/libs)

# Include the headers located at bibc/include
target_include_directories(bibc_shared PUBLIC bibc/include bibcxx/include bibcxx bibc bibfor/include $ENV{CONDA_PREFIX}/include $ENV{CONDA_PREFIX}/Library/include extern/include $ENV{CONDA_PREFIX}/Lib/site-packages/numpy/core/include build/std/debug)

# use the WINDOWS_EXPORT_ALL_SYMBOLS to create a .def file
set_target_properties(bibc_shared PROPERTIES
        WINDOWS_EXPORT_ALL_SYMBOLS TRUE
#        POSITION_INDEPENDENT_CODE ON
)