import os
import pathlib
from dataclasses import dataclass
from datetime import datetime

from error_patterns import find_error_pattern, ErrorPattern


@dataclass
class TestStats:
    name: str
    mess_file_path: pathlib.Path
    error_pattern: ErrorPattern


def check_mess_file(mess_file: str | pathlib.Path) -> TestStats:
    if isinstance(mess_file, str):
        mess_file = pathlib.Path(mess_file)

    data = mess_file.read_text(errors='replace', encoding='utf-8')

    error_pattern = find_error_pattern(data)

    return TestStats(
        name=mess_file.stem,
        mess_file_path=mess_file,
        error_pattern=error_pattern
    )


def eval_tests(test_dir: str | pathlib.Path):
    if isinstance(test_dir, str):
        test_dir = pathlib.Path(test_dir)

    tot_seq_files = 2204
    failed_test_files = list(test_dir.glob("*.mess"))
    tot_failed = len(failed_test_files)
    perc_passing = 100 - 100 * tot_failed / tot_seq_files
    tot_passing = tot_seq_files - tot_failed

    error_map = {}

    for mess_file in failed_test_files:
        error_data = check_mess_file(mess_file)
        if error_data.error_pattern not in error_map:
            error_map[error_data.error_pattern] = []
        error_map[error_data.error_pattern].append(error_data)

    err_str = f"Total passing tests: {perc_passing:.2f}% [{tot_passing}/{tot_seq_files}]\n"
    err_str += f"Total failed tests: {tot_failed} of {tot_seq_files}\n"
    # sort error_map by ref name of key object
    for error in sorted(error_map.keys(), key=lambda x: x.ref):
        failing_tests = error_map[error]
        perc_err = len(failing_tests) / tot_failed * 100
        files = '|'.join([f"{err.name}" for err in failing_tests])
        err_str += f"{error.ref} ('{error.description}'):\n  {len(failing_tests)} [{perc_err:.2f}%]\n  {files}\n"

    # save to file with todays date
    os.makedirs('results', exist_ok=True)
    today_str = datetime.now().strftime("%Y-%m-%d")
    with open(f"results/{today_str}.txt", "w") as f:
        f.write(err_str)


if __name__ == '__main__':
    eval_tests("../../temp/seq-debug")
