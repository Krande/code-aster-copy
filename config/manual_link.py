import pathlib
import subprocess

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / 'build' / "debug"


def bibc():
    bibc_objects = set(
        (BUILD_DIR / x).absolute() for x in (THIS_DIR / 'bibc_default_order.txt').read_text().splitlines())
    bibfor_objects = set(
        (BUILD_DIR / f'{x}.1.o').absolute() for x in (THIS_DIR / 'cshlib.txt').read_text().splitlines())

    # add all objects from bibfor/utilifor
    utilifor_files = (BUILD_DIR / "bibfor").rglob('*.o')
    bibfor_objects.update(set(x.absolute() for x in utilifor_files))

    rem_obj_by_name = ["ar_dgeqr2.F90", "ar_dtrevc.F90", "ar_dtrsen.F90"]
    rem_objects = set()
    for fp in bibfor_objects:
        if fp.name.startswith('ar_d'):
            rem_objects.add(fp)
        for obj in rem_obj_by_name:
            if obj in fp.name:
                rem_objects.add(fp)
    for robj in rem_objects:
        bibfor_objects.remove(robj)

    bibc_obj_list_path = ROOT_DIR / 'tmp_bibc_objects.txt'
    bibfor_obj_list_path = ROOT_DIR / 'tmp_bibfor_objects.txt'
    link_bat_path = ROOT_DIR / 'tmp_blink.bat'

    with open(bibc_obj_list_path, 'w') as f:
        f.write('\n'.join(map(str, bibc_objects)))
    with open(bibfor_obj_list_path, 'w') as f:
        f.write('\n'.join(map(str, bibfor_objects)))

    args = ['call_link.bat', '/nologo',
            '/MANIFEST',
            '/subsystem:console',
            '/IMPLIB:bibc\\bibc.lib',
            '/DLL',
            '/LIBPATH:C:\\work\\mambaforge\\envs\\codeaster-deps/libs',
            '/LIBPATH:C:\\work\\mambaforge\\envs\\codeaster-deps/include',
            '/LIBPATH:C:\\work\\code\\krande-code-aster-src\\bibfor\\include',
            f"@{bibc_obj_list_path.name}",
            f"@{bibfor_obj_list_path.name}",
            '/OUT:C:\\work\\code\\krande-code-aster-src\\build\\debug\\bibc\\bibc.dll',
            '/LIBPATH:C:/work/mambaforge/envs/codeaster-deps/Library/lib',
            '/LIBPATH:C:/work/mambaforge/envs/codeaster-deps/Library/bin',
            'MFrontGenericInterface.lib',
            'esmumps.lib',
            'lapack.lib',
            'smumps.lib',
            'scotch.lib',
            'scotcherr.lib',
            'metis.lib',
            'medC.lib',
            'hdf5.lib',
            'pthread.lib',
            '/DEBUG',
            # '/INCREMENTAL'
            ]

    # Now call that batch file
    subprocess.run(args, shell=True, cwd=ROOT_DIR)


def bibfor():
    ...


if __name__ == '__main__':
    bibc()
