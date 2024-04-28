# Based on https://learn.microsoft.com/en-us/previous-versions/visualstudio/visual-studio-2008/kkt2hd12(v=vs.90)?redirectedfrom=MSDN
# and https://stackoverflow.com/questions/362830/circular-dependencies-between-dlls-with-visual-studio

from config import CAMod, BUILD_DIR, TMP_DIR, THIS_DIR, CONDA_PREFIX_DIR, LIB_RAW_PREFIX, DEFOption

from msvc_utils import call_using_env


def run_lib(lib_name: CAMod | str, def_opt: DEFOption, use_wx=False):
    """Link the library using the LIB.exe command from the Visual Studio compiler
    https://learn.microsoft.com/en-us/previous-versions/visualstudio/visual-studio-2008/0xb6w1f8(v=vs.90)
    """
    if isinstance(lib_name, str):
        lib_name = CAMod(lib_name)

    module_dir = BUILD_DIR / str(lib_name.value)
    if lib_name == CAMod.BIBASTER:
        module_dir.mkdir(parents=True, exist_ok=True)
        files = [BUILD_DIR / "bibc" / "supervis" / "python.c.2.o"]
    elif lib_name == CAMod.BIBC:
        files = list(module_dir.rglob("*.o"))
        # remove the python.o file
        files_clean = [x for x in files if "python.c.2.o" not in str(x)]
        if len(files) - len(files_clean) != 1:
            raise ValueError(f"python.c.2.o not found in {files}")
        files = files_clean
    elif lib_name == CAMod.BIBCXX:
        files = set(module_dir.rglob("*.o"))
        if use_wx:
            to_be_removed = {
                "ConstantFieldOnCells.cxx.2.o",
                "ElementaryTerm.cxx.2.o",
                "FieldOnCells.cxx.2.o",
                "BehaviourDefinition.cxx.2.o",
                "ElementaryModeling.cxx.2.o",
                "FieldOnNodes.cxx.2.o",
            }
            files_cleaned = set([x for x in files if not any(y in str(x) for y in to_be_removed)])
            result = files_cleaned.intersection(to_be_removed)
            if len(result) > 0:
                raise ValueError(f"These files where not removed {result} from source")
            files = files_cleaned
    else:
        files = list(module_dir.rglob("*.o"))

    if def_opt == DEFOption.USE_BLANK_DEF:
        lib_file_name = f"{LIB_RAW_PREFIX}{lib_name.value}_blankdef"
    elif def_opt == DEFOption.USE_DEF:
        lib_file_name = f"{LIB_RAW_PREFIX}{lib_name.value}_def"
    else:  # NO_DEF
        lib_file_name = f"{LIB_RAW_PREFIX}{lib_name.value}_nodef"

    txt_file = TMP_DIR / "inputs" / f"{lib_file_name}_ofiles.txt"
    txt_file.parent.mkdir(exist_ok=True, parents=True)

    out_file = TMP_DIR / f"{lib_file_name}.lib"

    if len(files) == 0:
        raise ValueError(f"No files found in {module_dir}")

    with open(txt_file, "w") as f:
        f.write("\n".join(map(str, files)))

    cmd = [
        "LIB.exe",
        f"/LIBPATH:{CONDA_PREFIX_DIR}/libs",
        f"@{txt_file.as_posix()}",
        f"/OUT:{out_file.as_posix()}",
        "/MACHINE:X64",
        "/VERBOSE",
    ]
    if use_wx:
        cmd.append("/WX")

    if def_opt == DEFOption.USE_DEF:
        def_file = THIS_DIR / f"{lib_name.value}.def"
        if not def_file.exists():
            raise FileNotFoundError(f"{def_file} does not exist")
        cmd.append(f"/DEF:{def_file.as_posix()}")
    elif def_opt == DEFOption.USE_BLANK_DEF:
        cmd.append("/DEF")

    print(" ".join(cmd))

    result = call_using_env(cmd)
    if result.returncode != 0:
        # Error is reported in result.stdout
        stdout_multiline = "\n".join(result.stdout.splitlines())
        raise ValueError(f"Error linking {lib_name} due to {result.stderr=}, {stdout_multiline=}")


def main(lib_option: CAMod, def_opt: DEFOption):
    for lib_name in [CAMod.BIBFOR, CAMod.BIBCXX, CAMod.BIBC, CAMod.BIBASTER]:
        if lib_option != CAMod.ALL and lib_name != lib_option:
            continue
        run_lib(lib_name, def_opt)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Manually Link Code Aster libraries")
    parser.add_argument("--bibc", action="store_true", help="Create .lib for bibc")
    parser.add_argument("--bibfor", action="store_true", help="Create .lib for bibfor")
    parser.add_argument("--bibcxx", action="store_true", help="Create .lib for bibcxx")
    parser.add_argument("--bibaster", action="store_true", help="Create .lib for bibaster")
    parser.add_argument("--all", action="store_true", help="Create .lib for all")
    parser.add_argument("--use-def", action="store_true", help="Export specific symbols using LIB.exe")
    parser.add_argument("--use-blank-def", action="store_true", help="Export symbols using blank /DEF LIB.exe")

    args = parser.parse_args()

    if not any([args.bibc, args.bibfor, args.bibcxx, args.bibaster, args.all]):
        parser.error("No action requested, add --bibc, --bibfor, --bibcxx, or --all")

    if args.use_def:
        def_option = DEFOption.USE_DEF
    elif args.use_blank_def:
        def_option = DEFOption.USE_BLANK_DEF
    else:
        def_option = DEFOption.NO_DEF

    if args.all:
        main(CAMod.ALL, def_option)
    elif args.bibc:
        main(CAMod.BIBC, def_option)
    elif args.bibfor:
        main(CAMod.BIBFOR, def_option)
    elif args.bibcxx:
        main(CAMod.BIBCXX, def_option)
    elif args.bibaster:
        main(CAMod.BIBASTER, def_option)
    else:
        raise ValueError("Invalid state reached")
