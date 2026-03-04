import os
import pathlib
from time import strftime

THIS_DIR = pathlib.Path(__file__).parent


def main():
    src_dir = pathlib.Path(os.environ.get("SRC_DIR", THIS_DIR.parent.parent))
    chash = 'n/a'
    bname = 'n/a'

    version = os.getenv('PKG_VERSION', "17.0.99")
    # Only use first 3 components (major, minor, patch) — drop win_build suffix
    version_tuple = tuple([int(x) for x in version.split('.')[:3]])

    os.environ['_CA_VERSION'] = version

    with open(src_dir / 'code_aster/pkginfo.py', 'w') as f:
        curr_time = strftime("%d/%m/%Y")
        f.write(f'pkginfo =({version_tuple}, "{chash}", "{bname}", "{curr_time}", "n/a", 1, ["no source repository"],)')


if __name__ == '__main__':
    main()
