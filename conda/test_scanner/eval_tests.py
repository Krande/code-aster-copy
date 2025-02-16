import os
import pathlib
from dataclasses import dataclass
from datetime import datetime
import xml.etree.ElementTree as ET

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

@dataclass
class XmlLogStats:
    name: str
    tot_num_jobs: int
    tot_failures: int
    time: float
    failure_data: list[str] = None


def get_xml_log_stats(xml_file: pathlib.Path) -> XmlLogStats:
    # Parsing the XML
    root = ET.fromstring(xml_file.read_text(encoding="utf-8"))
    name = root.get("name")
    time = float(root.get("time"))
    tot_num_jobs = int(root.get("tests"))
    tot_failures = int(root.get("failures"))

    # Extracting test cases with failures
    testcases_with_failure = [testcase for testcase in root.findall("testcase") if testcase.find("failure") is not None]


    return XmlLogStats(name=name, tot_num_jobs=tot_num_jobs, tot_failures=tot_failures, time=time)

def eval_tests(test_dir: str | pathlib.Path, results_dir: str | pathlib.Path = "results", set_passing_env_var: bool = False):

    if isinstance(test_dir, str):
        test_dir = pathlib.Path(test_dir)

    test_cases_xml = get_xml_log_stats(test_dir / "run_testcases.xml")

    tot_seq_files = test_cases_xml.tot_num_jobs
    failed_test_files = list(test_dir.glob("*.mess"))
    tot_failed = len(failed_test_files)
    perc_passing = 100 - 100 * tot_failed / tot_seq_files
    if set_passing_env_var:
        os.environ["PASSED_TESTS"] = str(perc_passing)
        print(f"Set environment variable PASSED_TESTS to {perc_passing:.2f}%")
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
    with open(f"{results_dir}/{today_str}.txt", "w") as f:
        f.write(err_str)

def scan_cli():
    import argparse
    parser = argparse.ArgumentParser(description='Evaluate test results')
    parser.add_argument('test_dir', type=str, help='Directory containing test files')
    parser.add_argument('--output', type=str, help='Output directory for results')
    parser.add_argument('--set-passing-env-var', action='store_true', help='Set environment variable PASSED_TESTS which represents the percentage (0-100) of passing tests')
    args = parser.parse_args()
    eval_tests(args.test_dir, args.output, args.set_passing_env_var)

if __name__ == '__main__':
    #eval_tests("../../test_output")
    scan_cli()