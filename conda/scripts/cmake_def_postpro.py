import pathlib
from config import ROOT_DIR
from conda.scripts.run_cmake import CMAKE_BUILD_DIR
import re

RE_CAPTURE = re.compile(r"^(\s*)\?\?.*?@", re.MULTILINE)


def check_uniquness(data: str):
    unique_exports = set()
    for line in RE_CAPTURE.finditer(data):
        group = line.group(0)
        unique_exports.add(group)

    for unique_exp in unique_exports:
        print(unique_exp)


def main(exports_def, filtered_def):
    if not filtered_def.parent.exists():
        filtered_def.parent.mkdir(parents=True, exist_ok=True)

    data = exports_def.read_text(encoding="utf-8")

    # make a filtered copy of data
    with open(filtered_def, "w") as f:
        f.write("LIBRARY bibcxx\n")
        f.write("EXPORTS\n")

        for line in sorted(data.splitlines()):
            line_strip = line.strip()
            if line_strip == "EXPORTS" or line_strip.startswith("LIBRARY"):
                continue
            if line_strip.startswith("??$") or line_strip.startswith("?$"):
                continue
            if "@detail@" in line_strip and "@pybind11@" in line_strip:
                continue
            if "??_R0?AV?$_Ref_count_obj2" in line_strip and "?$basic_string" in line_strip:
                continue
            if "??_R0X@8" in line_strip:
                continue
            if line_strip.startswith("("):
                continue
            f.write(line)
            f.write("\n")


if __name__ == '__main__':
    # main(CMAKE_BUILD_DIR / "CMakeFiles" / "bibcxx_shared.dir" / "exports.def", ROOT_DIR / "bibcxx/bibcxx.def")
    main(pathlib.Path("temp") / "bibcxx_sym.def", pathlib.Path("temp") / "bibcxx_sym_filtered.def")
