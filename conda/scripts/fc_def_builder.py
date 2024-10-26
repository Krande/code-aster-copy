import pathlib
import json
import re

BIBFOR_DIR = pathlib.Path("../../bibfor").resolve().absolute()
BIBFOR_DEF = BIBFOR_DIR / "bibfor.def"

def scan_fortran_files(
    directory: pathlib.Path, cache_file: pathlib.Path = pathlib.Path("temp/fc_cache.json"), clear_cache=False
) -> dict:
    if cache_file.exists() and not clear_cache:
        with open(cache_file, "r") as cache:
            return json.load(cache)

    # Regex patterns to identify modules and subroutines/functions
    module_pattern = re.compile(r"^\s*module\s+(\w+)", re.IGNORECASE)
    func_pattern = re.compile(r"^\s*(subroutine|function)\s+(\w+)", re.IGNORECASE)

    # Dictionary to store modules and associated functions
    modules = {"na": []}

    # Scan each Fortran file in the directory
    for fp in directory.rglob("*"):
        if not fp.name.lower().endswith((".f90", ".f", ".for")):  # Common Fortran extensions
            continue

        with open(fp, "r", encoding="utf-8") as file:
            current_module = 'na'
            for line in file:
                if line.startswith("!"):
                    continue
                # Check if line defines a module
                module_match = module_pattern.match(line)
                if module_match:
                    current_module = module_match.group(1)
                    modules[current_module] = []

                # Check if line defines a function/subroutine
                func_match = func_pattern.match(line)
                if func_match:
                    function_name = func_match.group(2)
                    modules[current_module].append(function_name)
    # save to cache
    cache_file.parent.mkdir(exist_ok=True, parents=True)
    with open(cache_file, "w") as cache:
        json.dump(modules, cache, indent=4)

    return modules


def create_def_file(modules: dict, output_file: pathlib.Path, suffix="_", make_lower=True, make_upper=False):
    output_file.parent.mkdir(exist_ok=True, parents=True)
    functions_only = set([i.lower() for x in modules.values() for i in x])
    functions_sorted = sorted(functions_only)
    with open(output_file, "w") as def_file:
        def_file.write("LIBRARY bibfor\n")
        def_file.write("EXPORTS\n")
        for function in functions_sorted:
            if make_lower:
                def_file.write(f"    {function.lower()}{suffix}\n")
            elif make_upper:
                def_file.write(f"    {function.upper()}{suffix}\n")
            else:
                def_file.write(f"    {function}{suffix}\n")

def compare_with_existing_def(existing_def_file: pathlib.Path, modules: dict):
    functions_only = set([f"{i.lower()}_" for x in modules.values() for i in x])
    existing = set()
    for existing_def in existing_def_file.read_text().splitlines():
        if existing_def.startswith("LIBRARY") or existing_def.startswith("EXPORTS"):
            continue
        existing.add(existing_def.strip().lower())

    result_dict = {}
    missing = existing - functions_only
    print(f"Missing symbols: {len(missing)}")
    result_dict["missing"] = list(missing)

    unnecessary = functions_only - existing
    print(f"unnecessary symbols: {len(unnecessary)}")
    result_dict["unnecessary"] = list(unnecessary)

    matching = functions_only & existing
    print(f"Matching symbols: {len(matching)}")
    result_dict["matching"] = list(matching)

    with open("temp/fc_def_compare.json", "w") as f:
        json.dump(result_dict, f, indent=4)

def main():
    output_file = pathlib.Path("temp/output.def")  # Output .def file path
    modules = scan_fortran_files(BIBFOR_DIR, clear_cache=False)
    compare_with_existing_def(BIBFOR_DEF, modules)
    #create_def_file(modules, output_file)
    print(f".def file created at: {output_file}")


if __name__ == "__main__":
    main()
