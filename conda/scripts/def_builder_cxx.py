import pathlib
import json
import re


ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent.parent.absolute()
BIBCXX_DIR = ROOT_DIR / "bibcxx"
BIBCXX_DEF = BIBCXX_DIR / "bibcxx.def"


def scan_cpp_files(
        directory: pathlib.Path, cache_file: pathlib.Path = pathlib.Path("temp/cpp_cache.json"), clear_cache=False
) -> dict:
    if cache_file.exists() and not clear_cache:
        with open(cache_file, "r") as cache:
            return json.load(cache)

    # Regex patterns to identify functions and classes in C++
    namespace_pattern = re.compile(r"^\s*namespace\s+(\w+)", re.IGNORECASE)
    func_pattern = re.compile(r"^\s*([a-zA-Z_][\w:]*)\s+(\w+)\s*\(", re.IGNORECASE)
    class_pattern = re.compile(r"^\s*class\s+(\w+)", re.IGNORECASE)
    macro_pattern = re.compile(r"^\s*\w+\s*\(\s*\w+\s*,\s*(\w+)", re.IGNORECASE)  # Match macros like DEFSP

    # Dictionary to store namespaces, classes, and functions
    cpp_symbols = {"global": {"functions": [], "classes": []}}

    # Scan each C++ file in the directory
    for fp in directory.rglob("*"):
        # Match typical C++ file extensions
        if not fp.suffix.lower() in {".h", ".hpp", ".cpp", ".cxx", ".cc"}:
            continue
        if "deprecated" in list(fp.parents):
            continue
        with open(fp, "r", encoding="utf-8") as file:
            current_namespace = "global"
            current_class = None
            for line in file:
                # Check if line defines a namespace
                namespace_match = namespace_pattern.match(line)
                if namespace_match:
                    current_namespace = namespace_match.group(1)
                    if current_namespace not in cpp_symbols:
                        cpp_symbols[current_namespace] = {"functions": [], "classes": []}

                # Check if line defines a class
                class_match = class_pattern.match(line)
                if class_match:
                    current_class = class_match.group(1)
                    cpp_symbols[current_namespace]["classes"].append(current_class)

                # Check if line defines a function (ignoring method bodies and declarations)
                func_match = func_pattern.match(line)
                if func_match:
                    return_type, function_name = func_match.groups()
                    if current_class:
                        full_name = f"{current_class}::{function_name}"
                    else:
                        full_name = function_name
                    cpp_symbols[current_namespace]["functions"].append(full_name)

                # Check if line defines a macro containing function-like symbols
                macro_match = macro_pattern.match(line)
                if macro_match:
                    macro_function = macro_match.group(1)
                    cpp_symbols[current_namespace]["functions"].append(macro_function)

    # Save to cache
    cache_file.parent.mkdir(exist_ok=True, parents=True)
    with open(cache_file, "w") as cache:
        json.dump(cpp_symbols, cache, indent=4)

    return cpp_symbols


def create_def_file(cpp_symbols: dict, output_file: pathlib.Path, suffix="_", make_lower=True, make_upper=False):
    output_file.parent.mkdir(exist_ok=True, parents=True)

    functions_only = []
    for ns, data in cpp_symbols.items():
        functions_only.extend(data["functions"])

    functions_sorted = sorted(set(functions_only))

    with open(output_file, "w") as def_file:
        def_file.write("LIBRARY libcxx\n")
        def_file.write("EXPORTS\n")
        for function in functions_sorted:
            if make_lower:
                def_file.write(f"    {function.lower()}{suffix}\n")
            elif make_upper:
                def_file.write(f"    {function.upper()}{suffix}\n")
            else:
                def_file.write(f"    {function}{suffix}\n")


def main():
    output_file = BIBCXX_DEF.parent / "bibcxx_v2.def"
    cpp_symbols = scan_cpp_files(BIBCXX_DIR, clear_cache=True)

    create_def_file(cpp_symbols, output_file)
    print(f".def file created at: {output_file}")


if __name__ == "__main__":
    main()
