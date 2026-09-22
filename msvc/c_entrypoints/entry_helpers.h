//
// Created by Kristoffer on 19.06.2024.
//

#ifndef CODE_ASTER_SRC_ENTRY_HELPERS_H
#define CODE_ASTER_SRC_ENTRY_HELPERS_H

#include <Python.h>
#include <windows.h>
#include <string>
#include <vector>
#include <iostream>

// Function typedef for PyInit_aster
typedef PyObject* (*PyInit_func_t)(void);

// Retrieve the path of the current DLL
std::string GetDllPath();

// Retrieve the parent directory of the current DLL
std::string GetDllDirectory();

// Split paths in environment variables
std::vector<std::string> SplitPaths(const std::string &paths, char delimiter = ';');

// Load the specified DLL and retrieve the function pointer
HMODULE LoadDllAndGetFunction(const std::string& dllName, const std::string& funcName, PyInit_func_t& funcPtr);


#endif // CODE_ASTER_SRC_ENTRY_HELPERS_H
