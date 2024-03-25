import pathlib
import subprocess

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / 'build' / "debug"


def bibc():
    bibc_objects = set(
        (BUILD_DIR / x).absolute().as_posix() for x in (THIS_DIR / 'bibc_default_order.txt').read_text().splitlines())
    bibfor_objects = set(
        (BUILD_DIR / f'{x}.1.o').absolute().as_posix() for x in (THIS_DIR / 'cshlib.txt').read_text().splitlines())

    # add all objects from jevaux
    jeveux_files = (BUILD_DIR / "bibfor/jeveux").rglob('*.o')
    bibfor_objects.update(set(x.absolute().as_posix() for x in jeveux_files))

    # add all objects from bibfor/utilifor
    utilifor_files = (BUILD_DIR / "bibfor/utilifor").rglob('*.o')
    bibfor_objects.update(set(x.absolute().as_posix() for x in utilifor_files))

    bibc_obj_list_path = ROOT_DIR / 'tmp_bibc_objects.txt'
    bibfor_obj_list_path = ROOT_DIR / 'tmp_bibfor_objects.txt'
    link_bat_path = ROOT_DIR / 'tmp_blink.bat'

    with open(bibc_obj_list_path, 'w') as f:
        f.write('\n'.join(bibc_objects))
    with open(bibfor_obj_list_path, 'w') as f:
        f.write('\n'.join(bibfor_objects))

    args = ['call_link.bat', '/nologo',
            '/MANIFEST',
            '/nologo',
            '/MANIFEST',
            '/subsystem:console',
            '/IMPLIB:bibc\\bibc.lib',
            '/DLL',
            '/LIBPATH:C:\\work\\mambaforge\\envs\\codeaster-deps/libs',
            '/LIBPATH:C:\\work\\mambaforge\\envs\\codeaster-deps/include',
            f"@{bibc_obj_list_path.name}",
            f"@{bibfor_obj_list_path.name}",
            '/OUT:C:\\work\\code\\krande-code-aster-src\\build\\debug\\bibc\\bibc.dll',
            '/LIBPATH:C:/work/mambaforge/envs/codeaster-deps/Library/lib',
            'esmumps.lib',
            'scotch.lib',
            'scotcherr.lib',
            'metis.lib',
            'medC.lib',
            'hdf5.lib',
            'pthread.lib',
            '/DEBUG',
            # '/INCREMENTAL'
            ]
    # create a batch file that calls the linker
    with open(link_bat_path, 'w') as f:
        f.write('@echo off\n')
        f.write('setlocal\n')
        f.write(' '.join(args))
        f.write('\n')
        f.write('endlocal\n')

    # Now call that batch file
    subprocess.run(link_bat_path.name, shell=True, cwd=ROOT_DIR)


def bibfor():
    ...


if __name__ == '__main__':
    bibc()
