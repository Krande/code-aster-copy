# Use rglob(*.lib) on the conda_prefix/lib directory and create a symbol map for each library
import json
import logging
import os
import pathlib
import re
import subprocess
from dotenv import load_dotenv

load_dotenv()

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
CONDA_PREFIX_DIR = pathlib.Path(os.getenv(("CONDA_PREFIX")))
THIS_DIR = pathlib.Path(__file__).resolve().parent


def get_env():
    """Launch vcvarsall.bat and read the settings from its environment"""
    conda_env_bat = ROOT_DIR / "conda_env.bat"
    popen = subprocess.Popen(
        [conda_env_bat.as_posix(), "print"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ROOT_DIR
    )

    stdout, stderr = popen.communicate()
    if popen.wait() != 0:
        logging.error(stderr.decode("mbcs"))

    result = {}
    stdout = stdout.decode("mbcs")
    # use regex to find all lines after [vcvarsall.bat] Environment initialized for: 'x64'
    for line in stdout.split("\n"):
        if "=" not in line:
            continue
        line = line.strip()
        if line.lower().startswith(("rem", "(")):
            continue
        key, value = line.split("=", 1)
        result[key] = value

    return result


def process_dumpbin_output_static(output: str) -> dict:
    result = {}
    start_capturing = False  # Flag to start capturing after "Exports" section
    for line in output.splitlines():
        if line.strip() == "Exports":
            start_capturing = True
        elif start_capturing:
            parts = line.split()
            if len(parts) == 2 and parts[0].isdigit():
                # This line is assumed to have the format "ordinal name"
                # Adjust if the actual format varies
                symbol_name = parts[1]
                result[symbol_name] = {"ordinal": int(parts[0])}
            elif len(parts) == 1:
                # Handle symbols without ordinals
                symbol_name = parts[0]
                result[symbol_name] = {"ordinal": None}
    return result


def process_dumpbin_output_dynamic(output: str) -> dict:
    result = {}
    capturing = False  # Flag to start capturing symbols after "COFF SYMBOL TABLE"
    for line in output.splitlines():
        if "COFF SYMBOL TABLE" in line:
            capturing = True  # Start capturing after seeing "COFF SYMBOL TABLE"
        elif capturing and "|" in line:
            parts = line.split("|", 1)  # Split line at the first occurrence of '|'
            if len(parts) == 2:
                symbol_info, symbol_name = parts
                symbol_name = symbol_name.strip()  # Clean up symbol name
                if symbol_name not in result:
                    result[symbol_name] = []  # Initialize list if symbol is new
                result[symbol_name].append(symbol_info.strip())  # Add symbol info
    return result


def main():
    env_cache_file = THIS_DIR / ".env_cache.json"
    if not env_cache_file.exists():
        env = get_env()
        with open(env_cache_file, "w") as f:
            json.dump(env, f, indent=4)
    else:
        with open(env_cache_file, "r") as f:
            env = json.load(f)

    static_glob_map = {}
    dynamic_glob_map = {}
    for lib in CONDA_PREFIX_DIR.rglob("*.lib"):
        result_static = subprocess.run(["dumpbin", "/exports", lib], shell=True, capture_output=True, text=True, env=env)
        result_dynamic = subprocess.run(["dumpbin", "/symbols", lib], shell=True, capture_output=True, text=True,
                                       env=env)
        static_map_obj = process_dumpbin_output_static(result_static.stdout)
        dynamic_map_obj = process_dumpbin_output_dynamic(result_dynamic.stdout)
        static_glob_map[lib.as_posix()] = static_map_obj
        dynamic_glob_map[lib.as_posix()] = dynamic_map_obj

    # save glob_map to a json file
    print(static_glob_map)

    with open("lib_symbol_map_static.json", "w") as f:
        json.dump(static_glob_map, f, indent=4)

    with open("lib_symbol_map_dynamic.json", "w") as f:
        json.dump(dynamic_glob_map, f, indent=4)


if __name__ == "__main__":
    main()
