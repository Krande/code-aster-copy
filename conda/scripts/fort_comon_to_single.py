import re

from config import ROOT_DIR


# COMMON Block
# Beginning with the Intel Fortran Compiler (ifx) version 2023.2.0 and the Intel Fortran Complier Classic (ifort) version 2021.10.0,
# Programs that declare a COMMON block, instead of individual COMMON block variables, in an OpenMP data sharing clause cause a runtime failure, i.e. segmentation fault or incorrect result. The workaround is to declare the individual COMMON block variables.
# https://www.intel.com/content/www/us/en/developer/articles/release-notes/oneapi-fortran-compiler-release-notes.html

def iter_f90_files(common_block: str):
    bibfor_dir = ROOT_DIR / "bibfor"

    for fi in bibfor_dir.rglob("*.f90"):
        yield fi, fi.read_text('utf-8')


def replace_common_statements(input_text: str, target_line: str) -> str:
    # Regular expression to match the common statement
    pattern = re.escape(target_line)

    def replacer(match):
        common_block = match.group(0)
        common_name = common_block.split('/')[1]
        variables = common_block.split('/')[2].split(',')

        # Create a new common statement for each variable
        result = []
        for i, var in enumerate(variables):
            var = var.strip()
            new_line = f"common/{common_name}/{var}"
            if i > 0:
                new_line = ' ' * 4 + new_line
            result.append(new_line)

        return '\n'.join(result)

    updated_text = re.sub(pattern, replacer, input_text)
    return updated_text


def find_and_copy(common_block: str):
    for fi, txt in iter_f90_files(common_block):
        if common_block not in txt:
            continue
        updated_txt = replace_common_statements(txt, common_block)
        fi.write_text(updated_txt, encoding='utf-8')
        print(f"Updated file: {fi}")


if __name__ == '__main__':
    find_and_copy("common/ienvje/lbis, lois, lols, lor8, loc8")
