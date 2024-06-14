import os
import pathlib
import subprocess
from dataclasses import dataclass
from typing import List, Optional

CONDA_PREFIX = pathlib.Path(os.getenv('CONDA_PREFIX'))


@dataclass
class ExternalSymbolInfo:
    symbol_name: str
    original_so_filepath: str
    external_library: Optional[str]  # This might be None if the library can't be determined


def find_external_symbols_with_readelf(so_filepath: pathlib.Path, symbol_contains: str, search_all_modules=False) -> List[
    ExternalSymbolInfo]:
    if isinstance(so_filepath, str):
        so_filepath = pathlib.Path(so_filepath)

    # Use readelf to get undefined symbols
    readelf_symbols_cmd = ['readelf', '-Ws', so_filepath]
    symbols_proc = subprocess.run(readelf_symbols_cmd, capture_output=True, text=True, check=True)
    undefined_symbols = [line.split()[-1] for line in symbols_proc.stdout.split('\n') if 'UND' in line]

    # Filter symbols containing the specified substring
    filtered_symbols = [sym for sym in undefined_symbols if symbol_contains in sym]

    # Use readelf to get dynamic section info, including NEEDED libraries
    readelf_dynamic_cmd = ['readelf', '-d', so_filepath]
    dynamic_proc = subprocess.run(readelf_dynamic_cmd, capture_output=True, text=True, check=True)
    if search_all_modules:
        needed_libraries = CONDA_PREFIX.rglob("*.so")
    else:
        needed_libraries = [line.split('[')[-1].rstrip(']') for line in dynamic_proc.stdout.split('\n') if
                            'NEEDED' in line]
        needed_libraries = [CONDA_PREFIX / f"lib/{lib}" for lib in needed_libraries]

    # Attempt to find which library contains each symbol
    symbol_infos = []
    for symbol in filtered_symbols:
        found_library = None
        for lib_path in needed_libraries:
            if lib_path.stem ==so_file_path.stem:
                continue
            # This step is heuristic-based and assumes libraries are in standard locations
            # You may need to adjust the path or search mechanism for your environment
            readelf_lib_cmd = ['readelf', '-Ws', lib_path.as_posix()]
            try:
                readelf_lib_proc = subprocess.run(readelf_lib_cmd, capture_output=True, text=True, check=True)
                if symbol in readelf_lib_proc.stdout:
                    found_library = lib_path
                    break
            except subprocess.CalledProcessError:
                # Handle libraries that cannot be found or read
                pass
        # if found_library is None:
        #     raise ValueError("Unable to resolve the ")
        symbol_info = ExternalSymbolInfo(symbol_name=symbol, original_so_filepath=so_filepath,
                                         external_library=found_library)
        symbol_infos.append(symbol_info)

    return symbol_infos


if __name__ == '__main__':
    # Example usage
    so_file_path = CONDA_PREFIX / "lib/aster/libbibfor.so"
    symbol_contains = "r8prem"
    symbols_info = find_external_symbols_with_readelf(so_file_path, symbol_contains, True)
    for info in symbols_info:
        print(info)
