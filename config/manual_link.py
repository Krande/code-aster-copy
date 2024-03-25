import pathlib
import subprocess

ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent
THIS_DIR = pathlib.Path(__file__).resolve().parent
BUILD_DIR = ROOT_DIR / 'build' / "debug"

def main():
    # use subprocess to run batch file that activates the environment which will set the correct environment variables
    # and then run the build command
    bibc_objects = [(BUILD_DIR / x).absolute().as_posix() for x in (THIS_DIR / 'bibc_default_order.txt').read_text().splitlines()]
    bibfor_objects = [(BUILD_DIR / f'{x}.1.o').absolute().as_posix() for x in (THIS_DIR / 'cshlib.txt').read_text().splitlines()]
    args = ['call_link.bat', '/nologo',
            '/MANIFEST',
            '/nologo',
            '/MANIFEST',
            '/subsystem:console',
            '/IMPLIB:bibc\\bibc.lib',
            '/DLL',
            '/LIBPATH:C:\\work\\mambaforge\\envs\\codeaster-deps/libs',
            '/LIBPATH:C:\\work\\mambaforge\\envs\\codeaster-deps/include',
            *bibc_objects,
            *bibfor_objects,
            '/OUT:C:\\work\\code\\krande-code-aster-src\\build\\debug\\bibc\\bibc.dll',
            '/LIBPATH:C:/work/mambaforge/envs/codeaster-deps/Library/lib',
            'esmumps.lib',
            'scotch.lib',
            'scotcherr.lib',
            'metis.lib',
            'medC.lib',
            'hdf5.lib',
            'pthread.lib',
            '/DEBUG']
    # print(' '.join(args))
    subprocess.run(args, shell=True, cwd=ROOT_DIR)


if __name__ == '__main__':
    main()
