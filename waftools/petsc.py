# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

import os.path as osp
import re
from functools import partial
from waflib import Options, Configure, Logs, Utils, Errors

def options(self):
    group = self.add_option_group("Petsc library options")
    group.add_option('--disable-petsc', dest='enable_petsc',
                   action='store_false', default=None,
                   help='Disable PETSC support')
    group.add_option('--enable-petsc', dest='enable_petsc',
                   action='store_true', default=None,
                   help='Force PETSC support')
    group.add_option('--petsc-libs', type='string',
                   dest='petsc_libs', default=None,
                   help='petsc librairies used when linking')
    group.add_option('--embed-petsc', dest='embed_petsc',
                    default=False, action='store_true',
                    help='Embed PETSC libraries as static library')


def configure(self):
    if not self.env.BUILD_MPI:
        self.define('_DISABLE_PETSC', 1)
        self.undefine('HAVE_PETSC')
        self.define('_DISABLE_PETSC4PY', 1)
        self.define('_HAVE_PETSC4PY',0)
        return
    try:
        self.env.stash()
        self.check_petsc()
    except Errors.ConfigurationError:
        self.reset_msg()
        self.env.revert()
        self.define('_DISABLE_PETSC', 1)
        self.undefine('HAVE_PETSC')
        self.define('_DISABLE_PETSC4PY', 1)
        self.define('_HAVE_PETSC4PY',0)
        if self.options.enable_petsc:
            raise
    else:
        self.define('_HAVE_PETSC', 1)
        self.define('HAVE_PETSC', 1)
        self.check_petsc4py()

###############################################################################
@Configure.conf
def check_petsc(self):
    opts = self.options
    if opts.enable_petsc is False:
        raise Errors.ConfigurationError('PETSC disabled')

    optlibs = None
    if opts.petsc_libs is None:
        opts.petsc_libs = 'petsc'
        # add optional libs
        optlibs ='ml HYPRE superlu stdc++'
    if opts.petsc_libs:
        self.check_petsc_libs(optlibs)

    self.check_petsc_headers("petsc.h")
    self.check_petsc_headers("petscconf.h")
    self.check_petsc_version()
    self.check_sizeof_petsc_int()
    self.check_petsc_conf("PETSC_HAVE_ML", "HAVE_PETSC_ML")
    self.check_petsc_conf("PETSC_HAVE_HYPRE", "HAVE_PETSC_HYPRE")
    self.check_petsc_conf("PETSC_HAVE_SUPERLU", "HAVE_PETSC_SUPERLU")
    self.check_petsc_conf("PETSC_HAVE_MUMPS", "HAVE_PETSC_MUMPS")

@Configure.conf
def check_petsc_libs(self, optlibs):
    opts = self.options
    keylib = ('st' if opts.embed_all or opts.embed_scotch else '') + 'lib'
    for lib in Utils.to_list(optlibs or ''):
        self.check_cc(uselib_store='PETSC', use='MPI', uselib='PETSC',
                      mandatory=False, **{ keylib: lib})
    for lib in Utils.to_list(opts.petsc_libs):
        self.check_cc(uselib_store='PETSC', use='MPI', uselib='PETSC',
                      mandatory=True, **{ keylib: lib})

@Configure.conf
def check_petsc_headers(self, filename):
    check = partial(self.check, header_name=filename, use='PETSC MPI', uselib='SCOTCH Z',
                    uselib_store='PETSC')

    self.start_msg('Checking for header ' + filename)
    try:
        if not check(mandatory=False):
            if not check(includes=[osp.join(self.env.INCLUDEDIR, 'petsc')], mandatory=False):
                check(includes=[osp.join(self.env.OLDINCLUDEDIR, 'petsc')], mandatory=True)
    except:
        self.end_msg('no', 'YELLOW')
        raise
    else:
        self.end_msg('yes')

@Configure.conf
def check_petsc_version(self):
    fragment = r'''
#include <stdio.h>
#include <petsc.h>
int main(void){
#if defined(PETSC_VERSION_MAJOR) && defined(PETSC_VERSION_MINOR) && defined(PETSC_VERSION_SUBMINOR) && defined(PETSC_VERSION_PATCH)
    printf("PETSCVER: %d.%d.%d.%d", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, PETSC_VERSION_SUBMINOR, PETSC_VERSION_PATCH);
    return 0;
#endif
/* unexpected */
    return 1;
}'''
    self.start_msg('Checking petsc version')
    try:
        ret = self.check_cc(fragment=fragment, use='PETSC MPI', uselib='SCOTCH Z',
                            mandatory=True, execute=True, define_ret=True)
        mat = re.search('PETSCVER: *(?P<vers>[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)', ret)
        vers = mat and mat.group('vers')
        major, minor, sub, patch = [int(i) for i in vers.split('.')]
        vers = '%d.%d.%dp%d' % (major, minor, sub, patch)
        if major < 3 or (major == 3 and minor < 9):
            self.end_msg("unsupported petsc version: {0} "
                         "(expected 3.9.* or newer)".format(vers), 'RED')
            raise Errors.ConfigurationError
        self.define('ASTER_PETSC_VERSION', vers)
    except:
        self.end_msg('can not get version', 'RED')
        raise
    else:
        self.end_msg(vers)

@Configure.conf
def check_petsc4py(self):
    try:
        self.check_python_module('petsc4py')
        pymodule_path = self.get_python_variables(['petsc4py.get_include()'],
                                                  ['import petsc4py'])[0]
        self.env.append_unique('CYTHONFLAGS', '-I{0}'.format(pymodule_path))
    except Errors.ConfigurationError:
        self.define('_DISABLE_PETSC4PY', 1)
        self.define('_HAVE_PETSC4PY',0)
    else:
        self.undefine('_DISABLE_PETSC4PY')
        self.define('_HAVE_PETSC4PY', 1)

@Configure.conf
def check_sizeof_petsc_int(self):
    fragment = r'''
#include <stdio.h>
#include <petscconf.h>
int main(void) {
    printf("%d", (int)PETSC_SIZEOF_INT);
    return 0;
}'''
    self.code_checker('PETSC_INT_SIZE', self.check_cc, fragment,
                      'Checking size of PETSc integer',
                      'unexpected value for sizeof(petsc_int): %s',
                      into=(4, 8), use='PETSC')

@Configure.conf
def check_petsc_conf(self, petsc_var, aster_var):
    fragment = r'''
#include <stdio.h>
#include <petscconf.h>
int main(void) {{
#ifdef {name}
    printf("%d", {name});
    return 0;
#else
    printf("0");
    return 1;
#endif
}}'''.format(name=petsc_var)
    self.code_checker(aster_var, self.check_cc, fragment,
                      'Checking value of ' + aster_var, 'failure',
                      mandatory=False, setbool=True, use='PETSC')
