import pathlib
from dataclasses import dataclass
from enum import Enum


class FailCause(str, Enum):
    OverflowInt = "OverflowInt"
    Unknown = "Unknown"


@dataclass
class TestStats:
    name: str
    mess_file_path: pathlib.Path
    error_reason: FailCause


def check_mess_file(mess_file: str | pathlib.Path) -> TestStats:
    with open(mess_file, "r") as f:
        data = f.read()

    if "OverflowError: can't convert negative int to unsigned" in data:
        error_reason = FailCause.OverflowInt
    else:
        error_reason = FailCause.Unknown

    return TestStats(
        name=mess_file.stem,
        mess_file_path=mess_file,
        error_reason=error_reason
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
        if error_data.error_reason not in error_map:
            error_map[error_data.error_reason] = []
        error_map[error_data.error_reason].append(error_data)

    print(f"Total passing tests: {perc_passing:.2f}% [{tot_passing}/{tot_seq_files}]")
    overflow_errors = error_map.get(FailCause.OverflowInt)
    if overflow_errors:
        perc_overflow = len(overflow_errors) / tot_failed * 100
        print(f"Overflow errors: {len(overflow_errors)} [{perc_overflow:.2f}%]")
        print('|'.join([f"{err.name}" for err in overflow_errors]))
    else:
        print("No overflow errors")

    unknown_errors = error_map.get(FailCause.Unknown)
    if unknown_errors:
        perc_unknown = len(unknown_errors) / tot_failed * 100
        print(f"Unknown errors: {len(unknown_errors)} [{perc_unknown:.2f}%]")
        print('|'.join([f"{err.name}" for err in unknown_errors]))
    else:
        print("No unknown errors")


if __name__ == '__main__':
    eval_tests("../../temp/seq")
