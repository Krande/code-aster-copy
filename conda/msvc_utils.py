import pathlib
import subprocess

THIS_DIR = pathlib.Path(__file__).resolve().parent
ROOT_DIR = THIS_DIR.parent
BAT_CALL_FILE = ROOT_DIR / "call_compile.bat"


def call_using_env(args: list[str]) -> subprocess.CompletedProcess:
    command = f"call {BAT_CALL_FILE} {' '.join(map(str, args))}"

    # Execute the command
    result = subprocess.run(command, cwd=ROOT_DIR, shell=True, capture_output=True, text=True, encoding='utf-8')
    return result
