import re
from dataclasses import dataclass, field


@dataclass(frozen=True)
class ErrorPattern:
    ref: str
    contains: tuple[str | re.Pattern, ...] = field(default_factory=tuple)
    description: str = ""


ERROR_PATTERNS = [
    ErrorPattern(
        "OverflowError",
        ("OverflowError: can't convert negative int to unsigned",),
        "Overflowing integer"
    ),
    ErrorPattern(
        "MED_mpfprw",
        ("Erreur signalée dans la bibliothèque MED", "nom de l'utilitaire : mpfprw")
    ),
    ErrorPattern(
        "MED_mlclow",
        ("Erreur signalée dans la bibliothèque MED", "nom de l'utilitaire : mlclow")
    ),
    ErrorPattern(
        "MED_mfiope",
        ("Erreur signalée dans la bibliothèque MED", "nom de l'utilitaire : mfiope")
    ),
    ErrorPattern(
        "MED_18_SANS_PARTITIONNEUR",
        ("<MED_18>", "mot-clé SANS dans PARTITIONNEUR")
    ),
    ErrorPattern(
        "MED_57_RESU_READ_MED_F90",
        ("présence du champ demandé dans le fichier", "resuReadMed.F90")
    ),
    ErrorPattern(
        "MED_57_MEDFIELDN_C",
        ("présence du champ demandé dans le fichier", "MEDfieldnComponentByName.c", "Erreur à l'ouverture du groupe")
    ),
    ErrorPattern(
        "MED_57_LRVEMA",
        ("<MED_57>", "lrvema.F90")
    ),
    ErrorPattern(
        "DYNA_VIBRA_TYPE_ERROR_FLOAT_INT",
        (
            re.compile(
                r"DYNA_VIBRA[^\#\n]*(?:\n[^\#\n]*)*^TypeError: 'float' object cannot be interpreted as an integer",
                re.MULTILINE | re.DOTALL),
        ),
    ),
    ErrorPattern(
        "CALC_FERRAILLAGE_TYPE_ERROR_FLOAT_INT",
        (
            re.compile(
                r"CALC_FERRAILLAGE[^\#\n]*(?:\n[^\#\n]*)*^TypeError: 'float' object cannot be interpreted as an integer",
                re.MULTILINE | re.DOTALL),
        )
    ),
    ErrorPattern(
        "CALC_FERRAILLAGE_WIN_FATAL",
        (
            re.compile(
                r"CALC_FERRAILLAGE[^\#\n]*(?:\n[^\#\n]*)*^Windows fatal exception",
                re.MULTILINE | re.DOTALL),
        )
    ),
    ErrorPattern(
        "JEVEUX1_55_AFFE_CHAR_THER",
        (
            re.compile(r"AFFE_CHAR_THER[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in AFFE_CHAR_THER"
    ),
    ErrorPattern(
        "JEVEUX1_55_AFFE_CHAR_MECA",
        (
            re.compile(r"AFFE_CHAR_MECA[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in AFFE_CHAR_MECA"
    ),
    ErrorPattern(
        "JEVEUX1_55_POST_ELEM",
        (
            re.compile(r"POST_ELEM[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in POST_ELEM"
    ),
    ErrorPattern(
        "JEVEUX1_55_AFFE_MATERIAU",
        (
            re.compile(r"AFFE_MATERIAU[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in AFFE_MATERIAU"
    ),
    ErrorPattern(
        "JEVEUX1_55_MODI_MAILLAGE",
        (
            re.compile(r"MODI_MAILLAGE[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in MODI_MAILLAGE"
    ),
    ErrorPattern(
        "JEVEUX1_55_AFFE_MODELE",
        (
            re.compile(r"AFFE_MODELE[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in AFFE_MODELE"
    ),
    ErrorPattern(
        "JEVEUX1_55_AFFE_CARA_ELEM",
        (
            re.compile(r"AFFE_CARA_ELEM[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in AFFE_CARA_ELEM"
    ),
    ErrorPattern(
        "JEVEUX1_55_CREA_MAILLAGE",
        (
            re.compile(r"CREA_MAILLAGE[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in CREA_MAILLAGE"
    ),
    ErrorPattern(
        "JEVEUX1_55_AFFE_CHAR_CINE",
        (
            re.compile(r"AFFE_CHAR_CINE[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in AFFE_CHAR_CINE"
    ),
    ErrorPattern(
        "JEVEUX1_55_DEFI_GROUP",
        (
            re.compile(r"DEFI_GROUP[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in DEFI_GROUP"
    ),
    ErrorPattern(
        "JEVEUX1_55_IMPR_RESU",
        (
            re.compile(r"IMPR_RESU[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in IMPR_RESU"
    ),
    ErrorPattern(
        "JEVEUX1_55_POST_MAIL_XFEM",
        (
            re.compile(r"POST_MAIL_XFEM[^\#\n]*(?:\n[^\#\n]*)*<JEVEUX1_55>", re.MULTILINE | re.DOTALL),
        ),
        "JEVEUX1_55 error in POST_MAIL_XFEM"
    ),
    ErrorPattern(
        "JEVEUX_26_ORDR",
        ("<JEVEUX_26>", "Objet inexistant dans les bases ouvertes", ".ORDR")
    ),
    ErrorPattern(
        "ASRUN_PERMISSION_ERROR_WIN32",
        ("PermissionError: [WinError 32]", "asrun01b")
    ),
    ErrorPattern(
        "LAPACK_DLASWP_ERROR",
        ("erreur LAPACK (ou BLAS) au niveau de la routine  DLASWP",)
    ),
    ErrorPattern(
        "LAPACK_ZUNMQR_ERROR",
        ("erreur LAPACK (ou BLAS) au niveau de la routine  ZUNMQR",)
    ),
    ErrorPattern(
        "ASTERERROR_MATERIALFIELD_EMPTY",
        ("Raising C++ AsterError with id 'ValueError: MaterialField is empty'...",)
    ),
    ErrorPattern(
        "CALC_ESSAI_GEOMECA_PERM_ERROR",
        ("CALC_ESSAI_GEOMECA", "PermissionError",),
        "Permission denied error in CALC_ESSAI_GEOMECA"
    ),
    ErrorPattern(
        "CALC_ESSAI_GEOMECA_ACCESS_ERROR",
        ("CALC_ESSAI_GEOMECA", "Windows fatal exception: access violation",),
        "Access violation error in CALC_ESSAI_GEOMECA"
    ),
    ErrorPattern(
        "CALC_ESSAI_GEOMECA_INDEX_ERROR",
        ("CALC_ESSAI_GEOMECA", "IndexError",),
        "Index error in CALC_ESSAI_GEOMECA"
    ),
    ErrorPattern(
        "DVP1_jxveri",
        ("<DVP_1>", "jxveri.F90, ligne 98"),
        "DVP_1 jxveri exception"
    ),
    ErrorPattern(
        "DVP1_indlia",
        ("<DVP_1>", "indlia.F90, ligne 284"),
    ),
    ErrorPattern(
        "DEFI_CONT_OPS_EXIT_17",
        ("defi_cont_ops.py", "ABORT - exit code 17",),
        "Abort error in defi_cont_ops"
    ),
    ErrorPattern(
        "CREA_MAILLAGE_EXIT_17",
        ("crea_maillage.py", "ABORT - exit code 17",),
        "Abort error in crea_maillage"
    ),
    ErrorPattern(
        "UMAT_NO_WIN_SUPPORT",
        ("UMAT: Not available under Windows",)
    ),
    ErrorPattern(
        "FACTOR_55_MUMPS",
        ("<FACTOR_55>", "Problème ou alarme dans le solveur MUMPS",),
        "FACTOR_55 MUMPS problem or alarm"
    ),
    ErrorPattern(
        "MUMPS_WIN_FATAL",
        ("CALCUL MUMPS CENTRALISE", "Windows fatal exception: access violation")
    ),
    ErrorPattern(
        "PERMISSION_ERROR_macr_lign_coupe_ops",
        ("PermissionError: [Errno 13] Permission denied", "macr_lign_coupe_ops.py",),
        "Permission denied error in macr_lign_coupe_ops"
    ),
    ErrorPattern(
        "MED_ACCESS_ERROR",
        ("Windows fatal exception: access violation", "Création du fichier au format MED 3.3.1",),
        "Access violation error in MED file creation"
    ),
    ErrorPattern(
        "CATANESS_3_ALGELINE_1",
        ("<CATAMESS_3>", "L'alarme ALGELINE_1 est aggravée en erreur",),
        "ALGELINE_1 error in CATAMESS_3"
    ),
    ErrorPattern(
        "PROC_NA",
        ("No such file or directory: '/proc/23436/status'",),
        "No such file or directory error"
    ),
    ErrorPattern(
        "CALC_MATE_HOMO_ACCESS_ERROR",
        ("CALC_MATE_HOMO", "Windows fatal exception: access violation",),
        "Access violation error in CALC_MATE_HOMO"
    ),
    ErrorPattern(
        "MODI_MAILLAGE_WIN_FATAL",
        ("MODI_MAILLAGE", "Windows fatal exception: access violation")
    ),
    ErrorPattern(
        "MISSING_XMGRACE",
        (
            "Le fichier xmgrace n'existe pas.",
        ),
        "Missing xmgrace"
    ),
    ErrorPattern(
        "MISSING_HOMARD",
        (
            "Le fichier homard est inconnu",
        ),
        "Missing homard"
    ),
    ErrorPattern(
        "MISSING_SCIPY",
        ("ModuleNotFoundError: No module named 'scipy'",)
    ),
    ErrorPattern(
        "MISSING_MISS3D",
        ("run_miss3d", "The process cannot access the file")
    ),
    ErrorPattern(
        "MISSING_ASRUN",
        ("ModuleNotFoundError: No module named 'asrun'",)
    ),
    ErrorPattern(
        "FILE_EXISTS_WIN183",
        ("FileExistsError: [WinError 183]",)
    ),
    ErrorPattern(
        "NOOK_TEST_RESU",
        ("DIAGNOSTIC JOB : NOOK_TEST_RESU",)
    ),
    ErrorPattern(
        "WIN_FATAL_CODE",
        ("Windows fatal exception: code 0xc0000374",)
    ),
    ErrorPattern(
        "POST_ELEM_WIN_FATAL",
        ("POST_ELEM", "Windows fatal exception: access violation")
    ),
    ErrorPattern(
        "LIRE_RESU_WIN_FATAL",
        ("LIRE_RESU", "Windows fatal exception: access violation"),
        "Windows fatal exception: access violation error in LIRE_RESU"
    ),
    ErrorPattern(
        "FACTOR_11_MECA_STATIQUE",
        ("Problème : la matrice est singulière ou presque singulière",)
    ),
    ErrorPattern(
        "MECA_STATIQUE_WIN_FATAL",
        ("MECA_STATIQUE", "Windows fatal exception: access violation",)
    ),
    ErrorPattern(
        "REST_SPEC_PHYS",
        ("REST_SPEC_PHYS", "Windows fatal exception: access violation"),
    ),
    ErrorPattern(
        "STAT_NON_LINE_WIN_FATAL",
        ("STAT_NON_LINE", "Windows fatal exception: access violation"),
    ),
    ErrorPattern(
        "SIMU_POINT_MAT_WIN_FATAL",
        ("SIMU_POINT_MAT", "Windows fatal exception: access violation"),
    ),
    ErrorPattern(
        "SIMU_POINT_MAT_INDEX_ERROR",
        ("SIMU_POINT_MAT", "IndexError: stol argument out of range"),
    ),
    ErrorPattern(
        "NO_TEST_RESU",
        ("DIAGNOSTIC JOB : NO_TEST_RESU", "Table.destr")
    ),
    ErrorPattern(
        "PROC_ERROR",
        ("No such file or directory: '/proc/6980/status'",)
    ),
    ErrorPattern(
        "DEFI_BASE_REDUITE_WIN_FATAL",
        ("DEFI_BASE_REDUITE", "Windows fatal exception: access violation"),
    )
]


def find_error_pattern(data: str) -> ErrorPattern | None:
    for error_pattern in ERROR_PATTERNS:
        has_match = False
        contains_all = True
        for pattern in error_pattern.contains:
            if isinstance(pattern, re.Pattern):
                result = pattern.search(data)
                if result:
                    has_match = True
                else:
                    contains_all = False
            else:
                if pattern in data:
                    has_match = True
                else:
                    contains_all = False

        if has_match and contains_all:
            return error_pattern

    return ErrorPattern("Unknown", description="Unknown error pattern")
