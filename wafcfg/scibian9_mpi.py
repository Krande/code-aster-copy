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

"""
Configuration for Scibian 9  MPI

To be used from container, use before:
. /opt/public/scibian9_mpi.sh

or with native prerequisites:
. $HOME/dev/codeaster/devtools/etc/env_unstable_mpi.sh

./waf_mpi configure
./waf_mpi install
"""

import os
import official_programs


CT = True
if not os.getenv("SINGULARITY_NAME"):
    CT = False


def configure(self):
    if not CT:
        configure_legacy(self)
    else:
        configure_ct(self)


def configure_legacy(self):
    import scibian9_std
    YAMMROOT = scibian9_std.YAMMROOT
    opts = self.options

    opts.parallel = True
    scibian9_std.configure(self)
    self.env['ADDMEM'] = 500

    self.env.prepend_value('LIBPATH', [
        YAMMROOT + '/prerequisites/Parmetis_aster-403_aster3/lib',
        YAMMROOT + '/prerequisites/Scotch_aster-604_aster7/MPI/lib',
        YAMMROOT + '/prerequisites/Mumps-512_consortium_aster3/MPI/lib',
        YAMMROOT + '/prerequisites/Petsc_mpi-394_aster/lib',
    ])

    self.env.prepend_value('INCLUDES', [
        YAMMROOT + '/prerequisites/Parmetis_aster-403_aster3/include',
        YAMMROOT + '/prerequisites/Scotch_aster-604_aster7/MPI/include',
        YAMMROOT + '/prerequisites/Mumps-512_consortium_aster3/MPI/include',
        YAMMROOT + '/prerequisites/Petsc_mpi-394_aster/include',
    ])

    opts.enable_petsc = True
    self.env.append_value('LIB_METIS', ('parmetis'))
    self.env.append_value('LIB_SCOTCH', ('ptscotch','ptscotcherr','ptscotcherrexit'))


def configure_ct(self):
    PREREQ_PATH = "/opt/public/20190513-med41/gcc8-openblas-ompi2"
    opts = self.options

    opts.parallel = True

    official_programs.configure(self)
    opts.with_prog_xmgrace = False
    official_programs.check_prerequisites_package(self, PREREQ_PATH, '20190513')

    self.env['ADDMEM'] = 500

    TFELHOME = PREREQ_PATH + '/mfront-3.2.1'
    TFELVERS = '3.2.1'
    self.env.TFELHOME = TFELHOME
    self.env.TFELVERS = TFELVERS

    self.env.append_value('LIBPATH', [
        PREREQ_PATH + '/hdf5-1.10.3/lib',
        PREREQ_PATH + '/med-4.1.0/lib',
        PREREQ_PATH + '/metis-5.1.0_aster4/lib',
        PREREQ_PATH + '/parmetis-4.0.3_aster3/lib',
        PREREQ_PATH + '/scotch-6.0.4_aster7/lib',
        PREREQ_PATH + '/mumps-5.1.2_consortium_aster5/lib',
        PREREQ_PATH + '/mfront-3.2.1/lib',
        PREREQ_PATH + '/petsc-3.9.4_aster/lib',
        PREREQ_PATH + '/scalapack-2.1.0/lib',
    ])

    self.env.append_value('INCLUDES', [
        PREREQ_PATH + '/hdf5-1.10.3/include',
        PREREQ_PATH + '/med-4.1.0/include',
        PREREQ_PATH + '/metis-5.1.0_aster4/include',
        PREREQ_PATH + '/parmetis-4.0.3_aster3/include',
        PREREQ_PATH + '/scotch-6.0.4_aster7/include',
        PREREQ_PATH + '/mumps-5.1.2_consortium_aster5/include',
        PREREQ_PATH + '/mfront-3.2.1/include',
        PREREQ_PATH + '/petsc-3.9.4_aster/include',
    ])

    opts.maths_libs = ['openblas', 'scalapack']
    self.env.append_value('LIB_METIS', ('parmetis'))
    self.env.append_value('LIB', ('pthread', 'util'))
    self.env.append_value('LIB_SCOTCH', ('ptscotch','ptscotcherr','ptscotcherrexit'))
    # to fail if not found
    opts.enable_hdf5 = True
    opts.enable_med = True
    opts.enable_metis = True
    opts.enable_mumps = True
    opts.enable_scotch = True
    opts.enable_mfront = True

    opts.enable_petsc = True
