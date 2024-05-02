# Based on https://learn.microsoft.com/en-us/previous-versions/visualstudio/visual-studio-2008/kkt2hd12(v=vs.90)?redirectedfrom=MSDN
# and https://stackoverflow.com/questions/362830/circular-dependencies-between-dlls-with-visual-studio

from config import (
    CAMod,
    BUILD_DIR,
    TMP_DIR,
    THIS_DIR,
    CONDA_PREFIX_DIR,
    LIB_RAW_PREFIX,
    DEFOption,
    get_obj_list_path,
    get_lib_file,
    get_bibc_compile_files,
    get_bibaster_compile_files,
    get_bibcxx_compile_files,
    get_bibfor_compile_files,
    DEF_FILE_MAP,
    CompileStage,
)

from msvc_utils import call_using_env


def run_lib(lib_name: CAMod | str, def_opt: DEFOption, use_wx=False):
    """Link the library using the LIB.exe command from the Visual Studio compiler
    https://learn.microsoft.com/en-us/previous-versions/visualstudio/visual-studio-2008/0xb6w1f8(v=vs.90)
    """
    if isinstance(lib_name, str):
        lib_name = CAMod(lib_name)

    module_dir = BUILD_DIR / str(lib_name.value)
    if lib_name == CAMod.LIBASTER:
        files = get_bibaster_compile_files()
    elif lib_name == CAMod.BIBC:
        files = get_bibc_compile_files()
    elif lib_name == CAMod.BIBCXX:
        files = get_bibcxx_compile_files()
    else:
        files = get_bibfor_compile_files()

    txt_file = get_obj_list_path(str(lib_name.value), CompileStage.LIB)
    txt_file.parent.mkdir(exist_ok=True, parents=True)

    output_lib_file = get_lib_file(str(lib_name.value), def_opt)

    if len(files) == 0:
        raise ValueError(f"No files found in {module_dir}")

    cmd = [
        "LIB.exe",
        f"@{txt_file.as_posix()}",
    ]
    args = [
        "/MACHINE:X64",
        "/VERBOSE",
    ]
    if use_wx:
        args.append("/WX")

    if def_opt == DEFOption.USE_DEF:
        def_file = DEF_FILE_MAP.get(lib_name)
        if not def_file.exists():
            raise FileNotFoundError(f"{def_file} does not exist")
        args.append(f"/DEF:{def_file.as_posix()}")
    elif def_opt == DEFOption.USE_BLANK_DEF:
        args.append("/DEF")

    print(" ".join(cmd))
    with open(txt_file, "w") as f:
        f.write(f"/LIBPATH:{CONDA_PREFIX_DIR}/libs\n")
        f.write(f"/OUT:{output_lib_file.as_posix()}\n")
        f.write("\n".join(map(str, args)))
        f.write("\n")
        f.write("\n".join(map(str, files)))

    result = call_using_env(cmd)
    if result.returncode != 0:
        # Error is reported in result.stdout
        stdout_multiline = "\n".join(result.stdout.splitlines())
        raise ValueError(f"Error linking {lib_name} due to {result.stderr=}, {stdout_multiline=}")


def main(lib_option: CAMod, def_opt: DEFOption):
    for lib_name in [CAMod.BIBFOR, CAMod.BIBCXX, CAMod.BIBC, CAMod.LIBASTER]:
        if lib_option != CAMod.ALL and lib_name != lib_option:
            continue
        run_lib(lib_name, def_opt)


def cli():
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
        main(CAMod.LIBASTER, def_option)
    else:
        raise ValueError("Invalid state reached")


def manual():
    main(CAMod.ALL, DEFOption.USE_DEF)


if __name__ == "__main__":
    # cli()
    manual()
