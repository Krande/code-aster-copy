import ctypes.util
name = ctypes.util.find_library('aster')
lib = ctypes.cdll.LoadLibrary(name)
import aster

