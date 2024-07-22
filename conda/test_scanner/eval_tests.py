import os
import pathlib
from dataclasses import dataclass
from datetime import datetime
from enum import Enum


class TestResult(str, Enum):
    ErrOverflowInt = "OverflowInt"
    ErrUnknown = "Unknown"
    Passing = "Passing"
    MED_mpfprw = "MED_mpfprw"
    MED_mlclow = "MED_mlclow"
    MED_mfiope = "MED_mfiope"
    MED_ACCESS_ERROR = "MED_ACCESS_ERROR"
    CALC_MATE_HOMO_ACCESS_ERROR = "CALC_MATE_HOMO_ACCESS_ERROR"
    CATANESS_3_ALGELINE_1 = "CATANESS_3_ALGELINE_1"
    JEVEUX1_55 = "JEVEUX1_55"
    JEVEUX_26_ORDR = "JEVEUX_26_ORDR"
    FACTOR_55_MUMPS = "FACTOR_55_MUMPS"
    FLOAT_INT_ERROR = "FLOAT_INT_ERROR"
    PERMISSION_ERROR_WIN32 = "PERMISSION_ERROR_WIN32"
    PERMISSION_ERROR_macr_lign_coupe_ops = "PERMISSION_ERROR_macr_lign_coupe_ops"
    LAPACK_ERROR = "LAPACK_ERROR"
    CALC_ESSAI_GEOMECA_PERM_ERROR = "CALC_ESSAI_GEOMECA_PERM_ERROR"
    CALC_ESSAI_GEOMECA_ACCESS_ERROR = "CALC_ESSAI_GEOMECA_ACCESS_ERROR"
    CALC_ESSAI_GEOMECA_INDEX_ERROR = "CALC_ESSAI_GEOMECA_INDEX_ERROR"
    DVP1_jxveri = "DVP1_jxveri"
    DEFI_CONT_OPS = "DEFI_CONT_OPS"
    PROC_NA = "PROC_NA"

    NA_XMGRACE = "NA_XMGRACE"
    NA_HOMARD = "NA_HOMARD"


@dataclass
class TestStats:
    name: str
    mess_file_path: pathlib.Path
    error_reason: TestResult


def advanced_categorize(test_name: str, data: str) -> TestResult:
    if "OverflowError: can't convert negative int to unsigned" in data:
        return TestResult.ErrOverflowInt
    elif "Erreur signalée dans la bibliothèque MED" in data:
        if "nom de l'utilitaire : mpfprw" in data:
            return TestResult.MED_mpfprw
        elif "nom de l'utilitaire : mlclow" in data:
            return TestResult.MED_mlclow
        elif "nom de l'utilitaire : mfiope" in data:
            return TestResult.MED_mfiope
        else:
            raise ValueError(f"Unknown MED error in {test_name}")
    elif "<JEVEUX1_55>" in data:
        if "Un écrasement aval est détecté, la zone mémoire" in data:
            return TestResult.JEVEUX1_55
        else:
            raise ValueError(f"Unknown JEVEUX1_55 error in {test_name}")
    elif "TypeError: 'float' object cannot be interpreted as an integer" in data:
        return TestResult.FLOAT_INT_ERROR
    elif "PermissionError: [WinError 32]" in data:
        return TestResult.PERMISSION_ERROR_WIN32
    elif "erreur LAPACK (ou BLAS) au niveau de la routine  DLASWP" in data:
        return TestResult.LAPACK_ERROR
    elif "CALC_ESSAI_GEOMECA" in data:
        if "PermissionError" in data:
            return TestResult.CALC_ESSAI_GEOMECA_PERM_ERROR
        elif "Windows fatal exception: access violation" in data:
            return TestResult.CALC_ESSAI_GEOMECA_ACCESS_ERROR
        elif "IndexError" in data:
            return TestResult.CALC_ESSAI_GEOMECA_INDEX_ERROR
        else:
            raise ValueError(f"Unknown CALC_ESSAI_GEOMECA error in {test_name}")
    elif "<EXCEPTION> <DVP_1>" in data and "jxveri.F90, ligne 98" in data:
        return TestResult.DVP1_jxveri
    elif "defi_cont_ops.py" in data and "ABORT - exit code 17" in data:
        return TestResult.DEFI_CONT_OPS
    elif "<EXCEPTION> <JEVEUX_26>" in data and "Objet inexistant dans les bases ouvertes : 0000000a           .ORDR" in data:
        return TestResult.JEVEUX_26_ORDR
    elif "<FACTOR_55>" in data and "Problème ou alarme dans le solveur MUMPS" in data:
        return TestResult.FACTOR_55_MUMPS
    elif "PermissionError: [Errno 13] Permission denied" in data and "macr_lign_coupe_ops.py" in data:
        return TestResult.PERMISSION_ERROR_macr_lign_coupe_ops
    elif "Windows fatal exception: access violation" in data and "Création du fichier au format MED 3.3.1" in data:
        return TestResult.MED_ACCESS_ERROR
    elif "Le fichier xmgrace n'existe pas." in data:
        return TestResult.NA_XMGRACE
    elif "Le fichier homard est inconnu" in data:
        return TestResult.NA_HOMARD
    elif "<CATAMESS_3>" in data and "L'alarme ALGELINE_1 est aggravée en erreur" in data:
        return TestResult.CATANESS_3_ALGELINE_1
    elif "No such file or directory: '/proc/23436/status'" in data:
        return TestResult.PROC_NA
    elif "CALC_MATE_HOMO" in data and "Windows fatal exception: access violation" in data:
        return TestResult.CALC_MATE_HOMO_ACCESS_ERROR
    else:
        return TestResult.ErrUnknown


def check_mess_file(mess_file: str | pathlib.Path) -> TestStats:
    if isinstance(mess_file, str):
        mess_file = pathlib.Path(mess_file)

    with open(mess_file, "r", errors='replace', encoding='utf-8') as f:
        data = f.read()

    error_reason = advanced_categorize(mess_file.stem, data)

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

    overflow_errors = error_map.get(TestResult.ErrOverflowInt)
    err_str = f"Total passing tests: {perc_passing:.2f}% [{tot_passing}/{tot_seq_files}]\n"
    err_str += f"Total failed tests: {tot_failed} of {tot_seq_files}\n"
    if overflow_errors:
        perc_overflow = len(overflow_errors) / tot_failed * 100
        err_str += f"Overflow errors: {len(overflow_errors)} [{perc_overflow:.2f}%]\n"
        err_str += '|'.join([f"{err.name}" for err in overflow_errors]) + '\n'
    else:
        err_str += "No overflow errors\n"

    unknown_errors = error_map.get(TestResult.ErrUnknown)
    if unknown_errors:
        perc_unknown = len(unknown_errors) / tot_failed * 100
        err_str += f"Unknown errors: {len(unknown_errors)} [{perc_unknown:.2f}%]\n"
        err_str += '|'.join([f"{err.name}" for err in unknown_errors]) + '\n'
    else:
        err_str += "No unknown errors\n"

    err_str += "Specific errors:\n"
    for key, value in error_map.items():
        if not value or key in [TestResult.ErrOverflowInt, TestResult.ErrUnknown]:
            continue
        perc_err = len(value) / tot_failed * 100
        files = '|'.join([f"{err.name}" for err in value])
        err_str += f"{key.value}: {len(value)} [{perc_err:.2f}%] | {files}\n"

    # save to file with todays date
    os.makedirs('results', exist_ok=True)
    today_str = datetime.now().strftime("%Y-%m-%d")
    with open(f"results/{today_str}.txt", "w") as f:
        f.write(err_str)


if __name__ == '__main__':
    eval_tests("../../temp/seq-debug2")
