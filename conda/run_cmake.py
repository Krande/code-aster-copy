from msvc_utils import call_using_env, ROOT_DIR

CMAKE_BUILD_DIR = ROOT_DIR / "build-cmake"
if not CMAKE_BUILD_DIR.exists():
    CMAKE_BUILD_DIR.mkdir(exist_ok=True, parents=True)


def run_cmake(module_name: str):
    cmd = ["cmake", "..", "-G", "Ninja", f"-DBUILD_{module_name.upper()}=ON", "-DCMAKE_BUILD_TYPE=release", "--fresh"]
    result = call_using_env(cmd, cwd=CMAKE_BUILD_DIR)
    if result.returncode != 0:
        print(result.stderr)
        print(result.stdout)
        raise ValueError(f"Error running cmake due to {result.stderr=}")
    else:
        print("CMake ran successfully")

    cmd = ["ninja"]
    result = call_using_env(cmd, cwd=CMAKE_BUILD_DIR)
    if result.returncode != 0:
        print(result.stderr)
        print(result.stdout)
        raise ValueError(f"Error running ninja due to {result.stderr=}")
    else:
        print("Ninja ran successfully")


def cli():
    import argparse

    parser = argparse.ArgumentParser(description="Run CMake and Ninja")
    parser.add_argument("--module-name", type=str, help="Name of the module")
    args = parser.parse_args()
    run_cmake(args.module_name)


if __name__ == "__main__":
    cli()
