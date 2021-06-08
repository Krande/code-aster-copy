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
Configuration for cronos

./waf configure
./waf install -p
"""

import os
import official_programs
ROOT = os.environ['PREREQ_PATH']


def configure(self):
    opts = self.options

    official_programs.configure(self)
    official_programs.check_prerequisites_package(self, ROOT, "20190513-med41")
    opts.with_prog_salome = False
    opts.with_prog_europlexus = False
    opts.with_prog_xmgrace = False

    self.env['ADDMEM'] = 700
    self.env.append_value('OPT_ENV', [
        '. /etc/profile.d/modules.sh',
        'module use /home/G79848/env/modulefiles /fscronos/software/shared/modules/all /usr/share/easybuild/modules/all',
        'module load imkl/2020.4.304',
        '',
        'export LD_PRELOAD=\\',
        '/fscronos/software/shared/easybuild/software/imkl/2020.4.304/mkl/lib/intel64/libmkl_gf_lp64.so:\\',
        '/fscronos/software/shared/easybuild/software/imkl/2020.4.304/mkl/lib/intel64/libmkl_gnu_thread.so:\\',
        '/fscronos/software/shared/easybuild/software/imkl/2020.4.304/mkl/lib/intel64/libmkl_core.so:\\',
        '/lib64/libgomp.so.1',
        ])
    self.env.append_value('OPT_ENV_FOOTER', [
        'export PATH=$PATH:' + ROOT + '/grace-0.0.1/bin'
    ])

    TFELHOME = ROOT + '/mfront-3.2.1'
    TFELVERS = '3.2.1'
    self.env.TFELHOME = TFELHOME
    self.env.TFELVERS = TFELVERS

    self.env.append_value('LIBPATH', [
        ROOT + '/hdf5-1.10.3/lib',
        ROOT + '/med-4.1.0/lib',
        ROOT + '/metis-5.1.0_aster4/lib',
        ROOT + '/scotch-6.0.4_aster7/lib',
        ROOT + '/mumps-5.1.2_consortium_aster5/lib',
        TFELHOME + '/lib',
    ])

    self.env.append_value('INCLUDES', [
        ROOT + '/hdf5-1.10.3/include',
        ROOT + '/med-4.1.0/include',
        ROOT + '/metis-5.1.0_aster4/include',
        ROOT + '/scotch-6.0.4_aster7/include',
        ROOT + '/mumps-5.1.2_consortium_aster5/include',
        ROOT + '/mumps-5.1.2_consortium_aster5/include_seq',
        TFELHOME + '/include',
    ])

    self.env.append_value('LIB_SCOTCH', ('scotcherrexit'))

    # to fail if not found
    opts.enable_hdf5 = True
    opts.enable_med = True
    opts.enable_metis = True
    opts.enable_mumps = True
    opts.enable_scotch = True
    opts.enable_mfront = True

    opts.enable_petsc = False
