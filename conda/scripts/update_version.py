import os
import pathlib
from time import strftime

THIS_DIR = pathlib.Path(__file__).parent


def main():
    src_dir = THIS_DIR.parent.parent
    # Commit hash and branch name set to n/a for now
    chash = 'n/a'
    bname = 'n/a'

    version = os.getenv('PKG_VERSION', "17.0.99")
    version_tuple = tuple([int(x) for x in version.split('.')])

    os.environ['_CA_VERSION'] = version

    with open(src_dir / 'code_aster/pkginfo.py', 'w') as f:
        curr_time = strftime("%d/%m/%Y")
        f.write(f'pkginfo =({version_tuple}, "{chash}", "{bname}", "{curr_time}", "n/a", 1, ["no source repository"],)')


if __name__ == '__main__':
    main()
