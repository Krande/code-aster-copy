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
Configuration for Scibian 9

To be used from container, use before:
. /opt/public/scibian9_std.sh

or with native prerequisites:
. $HOME/dev/codeaster/devtools/etc/env_unstable.sh

./waf configure
./waf install
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
    YAMMROOT = os.environ['ROOT_SALOME']
    opts = self.options

    official_programs.configure(self)
    official_programs.check_prerequisites_package(self, YAMMROOT, '20190513')

    self.env.append_value('CXXFLAGS', ['-D_GLIBCXX_USE_CXX11_ABI=0'])
    self.env['ADDMEM'] = 350

    TFELHOME = YAMMROOT + '/prerequisites/Mfront-TFEL321_aster'
    TFELVERS = '3.2.1'
    self.env.TFELHOME = TFELHOME
    self.env.TFELVERS = TFELVERS

    self.env.append_value('LIBPATH', [
        '/opt/hdf5/1.10.3/lib',
        YAMMROOT + '/prerequisites/Medfichier-400/lib',
        YAMMROOT + '/prerequisites/Metis_aster-510_aster4/lib',
        YAMMROOT + '/prerequisites/Scotch_aster-604_aster7/SEQ/lib',
        YAMMROOT + '/prerequisites/Mumps-512_consortium_aster3/SEQ/lib',
        TFELHOME + '/lib',
    ])

    self.env.append_value('INCLUDES', [
        '/opt/hdf5/1.10.3/include',
        YAMMROOT + '/prerequisites/Medfichier-400/include',
        YAMMROOT + '/prerequisites/Metis_aster-510_aster4/include',
        YAMMROOT + '/prerequisites/Scotch_aster-604_aster7/SEQ/include',
        YAMMROOT + '/prerequisites/Mumps-512_consortium_aster3/SEQ/include',
        YAMMROOT + '/prerequisites/Mumps-512_consortium_aster3/SEQ/include_seq',
        TFELHOME + '/include',
    ])

    self.env.append_value('LIB', ('pthread', 'util'))
    self.env.append_value('LIB_SCOTCH', ('scotcherrexit'))
    # to fail if not found
    opts.enable_hdf5 = True
    opts.enable_med = True
    opts.enable_metis = True
    opts.enable_mumps = True
    opts.enable_scotch = True
    opts.enable_mfront = True

    opts.enable_petsc = False


def configure_ct(self):
    PREREQ_PATH = "/opt/public/20190513-med41/gcc8-openblas-seq"
    opts = self.options

    official_programs.configure(self)
    opts.with_prog_xmgrace = False
    official_programs.check_prerequisites_package(self, PREREQ_PATH, '20190513')

    self.env['ADDMEM'] = 350

    TFELHOME = PREREQ_PATH + '/mfront-3.2.1'
    TFELVERS = '3.2.1'
    self.env.TFELHOME = TFELHOME
    self.env.TFELVERS = TFELVERS

    self.env.append_value('LIBPATH', [
        PREREQ_PATH + '/hdf5-1.10.3/lib',
        PREREQ_PATH + '/med-4.1.0/lib',
        PREREQ_PATH + '/metis-5.1.0_aster4/lib',
        PREREQ_PATH + '/scotch-6.0.4_aster7/lib',
        PREREQ_PATH + '/mumps-5.1.2_consortium_aster5/lib',
        PREREQ_PATH + '/mfront-3.2.1/lib',
    ])

    self.env.append_value('INCLUDES', [
        PREREQ_PATH + '/hdf5-1.10.3/include',
        PREREQ_PATH + '/med-4.1.0/include',
        PREREQ_PATH + '/metis-5.1.0_aster4/include',
        PREREQ_PATH + '/scotch-6.0.4_aster7/include',
        PREREQ_PATH + '/mumps-5.1.2_consortium_aster5/include',
        PREREQ_PATH + '/mumps-5.1.2_consortium_aster5/include_seq',
        PREREQ_PATH + '/mfront-3.2.1/include',
    ])

    self.env.append_value('LIB', ('pthread', 'util'))
    self.env.append_value('LIB_SCOTCH', ('scotcherrexit'))
    # to fail if not found
    opts.enable_hdf5 = True
    opts.enable_med = True
    opts.enable_metis = True
    opts.enable_mumps = True
    opts.enable_scotch = True
    opts.enable_mfront = True

    opts.enable_petsc = False
