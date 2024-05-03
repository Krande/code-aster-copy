from msvc_utils import call_using_env, ROOT_DIR

CMAKE_BUILD_DIR = ROOT_DIR / "build-cmake"
if not CMAKE_BUILD_DIR.exists():
    CMAKE_BUILD_DIR.mkdir()


def main():
    cmd = ["cmake", "..", "-G", "Ninja", "-DBUILD_BIBCXX=ON", "-DCMAKE_BUILD_TYPE=Release", "--fresh"]
    result = call_using_env(cmd, cwd=CMAKE_BUILD_DIR)
    if result.returncode != 0:
        print(result.stderr)
        print(result.stdout)
        raise ValueError(f"Error running cmake due to {result.stderr=}")

    cmd = ["ninja"]
    result = call_using_env(cmd, cwd=CMAKE_BUILD_DIR)
    if result.returncode != 0:
        print(result.stderr)
        print(result.stdout)
        raise ValueError(f"Error running ninja due to {result.stderr=}")


if __name__ == "__main__":
    main()
