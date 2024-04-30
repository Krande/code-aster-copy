project(mylibrary)

cmake_minimum_required(VERSION 3.28)

# Find all .cxx files recursively in the subdirectory bibcxx
file(GLOB_RECURSE mylibrary_SOURCES "bibcxx/*.cxx")

# Create a shared library from the sources
add_library(bibcxx_shared SHARED ${mylibrary_SOURCES})

# Add /std:c++17 to the compiler flags
target_compile_options(bibcxx_shared PRIVATE /std:c++17 /nologo /MD /DTFEL_ARCH64)

# Add lib path to the linker flags
target_link_options(bibcxx_shared PRIVATE /LIBPATH:$ENV{CONDA_PREFIX}/libs)

# Include the headers located at bibcxx/include
target_include_directories(bibcxx_shared PUBLIC bibcxx/include bibcxx bibc/include bibfor/include $ENV{CONDA_PREFIX}/include $ENV{CONDA_PREFIX}/Library/include extern/include $ENV{CONDA_PREFIX}/Lib/site-packages/numpy/core/include)

# use the WINDOWS_EXPORT_ALL_SYMBOLS to create a .def file
set_target_properties(bibcxx_shared PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
