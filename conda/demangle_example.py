import pathlib

from cpp_demangle import demangle

from conda.manual_create_def_file import iter_symbol_names
from conda.msvc_utils import call_using_env


def main():
    sym_file = pathlib.Path('temp/bibcxx_sym.def').resolve().absolute().as_posix()
    dest_file = sym_file.replace('.def', '_demangled.txt')
    result = call_using_env(['undname', sym_file])
    with open(dest_file, 'w') as f:
        f.write(result.stdout)

    # for sym_name in iter_symbol_names("temp/_raw_bibcxx_sym.txt"):
    #     try:
    #         result = demangle(sym_name)
    #     except ValueError:
    #         continue
    #     print(result)


if __name__ == "__main__":
    main()
