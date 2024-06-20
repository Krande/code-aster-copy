#include "entry_helpers.h"

std::string GetDllPath()
{
    HMODULE hModule = nullptr;
    char path[MAX_PATH];

    // Get the handle to the current DLL module
    if (GetModuleHandleEx(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                          reinterpret_cast<LPCSTR>(&GetDllPath), &hModule))
    {
        // Get the path of the DLL
        if (GetModuleFileName(hModule, path, MAX_PATH) > 0)
        {
            return std::string(path);
        }
    }
    return std::string();
}

std::string GetDllDirectory()
{
    std::string dllPath = GetDllPath();
    if (!dllPath.empty())
    {
        size_t pos = dllPath.find_last_of("\\/");
        if (pos != std::string::npos)
        {
            return dllPath.substr(0, pos);
        }
    }
    return std::string();
}

std::vector<std::string> SplitPaths(const std::string &paths, char delimiter)
{
    std::vector<std::string> pathList;
    size_t start = 0;
    size_t end = paths.find(delimiter);
    while (end != std::string::npos)
    {
        pathList.push_back(paths.substr(start, end - start));
        start = end + 1;
        end = paths.find(delimiter, start);
    }
    pathList.push_back(paths.substr(start));
    return pathList;
}

HMODULE LoadDllAndGetFunction(const std::string& dllName, const std::string& funcName, PyInit_func_t& funcPtr)
{
    std::string dllDirectory = GetDllDirectory();
    char dllPath[MAX_PATH];
    snprintf(dllPath, MAX_PATH, "%s\\%s", dllDirectory.c_str(), dllName.c_str());

    // Print the path to the console
    printf("Loading %s from %s\n", dllName.c_str(), dllPath);

    // Add the DLL directory to the search path
    SetDllDirectory(dllDirectory.c_str());

    // Load the DLL
    HMODULE hDll = LoadLibrary(dllPath);
    if (!hDll) {
        DWORD error = GetLastError();
        LPVOID lpMsgBuf;
        FormatMessage(
            FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
            NULL,
            error,
            MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
            (LPTSTR)&lpMsgBuf,
            0, NULL);
        std::cerr << "Could not load " << dllName << ". Error: " << (char *)lpMsgBuf << std::endl;
        LocalFree(lpMsgBuf);
        return nullptr;
    }

    // Get the function pointer
    funcPtr = (PyInit_func_t)GetProcAddress(hDll, funcName.c_str());
    if (!funcPtr)
    {
        // Try using the decorated name if necessary
        std::string decoratedFuncName = "?" + funcName + "@@YAPEAU_object@@XZ";
        funcPtr = (PyInit_func_t)GetProcAddress(hDll, decoratedFuncName.c_str());
        if (!funcPtr)
        {
            std::cerr << "Could not find " << funcName << " in " << dllName << std::endl;
            return nullptr;
        }
    }

    // Reset the DLL directory to avoid side effects
    SetDllDirectory(NULL);

    return hDll;
}
