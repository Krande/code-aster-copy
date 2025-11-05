import pathlib
import re

THIS_DIR = pathlib.Path(__file__).parent
ROOT_DIR = THIS_DIR.parent.parent


def add_bind_c_directive(code: str) -> str:
    # Regular expression to find subroutine definitions
    regex = r"(subroutine\s+(\w+)\s*\(([^)]*)\))"

    # Function to replace each found subroutine definition with the updated one including BIND(C)
    def replacer(match):
        subroutine_declaration = match.group(1)  # Full declaration e.g., subroutine mfafai(fid, maa, ind, ...)
        subroutine_name = match.group(2).upper()  # Subroutine name, e.g., mfafai
        # Return the modified subroutine declaration
        print(f"Found subroutine: {subroutine_name}")
        return f"""
#ifdef ASTER_PLATFORM_MSVC64
    {subroutine_declaration} BIND(C, name='{subroutine_name}')
#else
    {subroutine_declaration}
#endif"""

    # Apply the regular expression to add BIND(C) to all subroutine declarations
    updated_code = re.sub(regex, replacer, code, flags=re.DOTALL)
    return updated_code


def main():
    bibfor_med_incl = ROOT_DIR / 'bibfor' / 'include' / "med"

    for header in bibfor_med_incl.rglob('*.h'):
        with open(header, 'r') as file:
            code = file.read()

        updated_code = add_bind_c_directive(code)

        if updated_code != code:
            with open(header, 'w') as file:
                file.write(updated_code)


if __name__ == '__main__':
    main()
