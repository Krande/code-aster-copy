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

import code_aster
from code_aster.Commands import *
from code_aster import MPI

code_aster.init("--test")

test = code_aster.TestCase()

#parallel=False
parallel = True

if parallel:
    rank=MPI.COMM_WORLD.Get_rank()
    MAIL = code_aster.ParallelMesh()
    MAIL.readMedFile("petsc04g/%d.med"%rank, True)
else:
    MAIL = code_aster.Mesh()
    MAIL.readMedFile("petsc04a.mmed")

model = AFFE_MODELE(
                    AFFE=_F(MODELISATION=('3D_INCO_UP', ), PHENOMENE='MECANIQUE', TOUT='OUI'),
                    MAILLAGE=MAIL,
                )

mater = DEFI_MATERIAU(ELAS=_F(E=1.0, NU=0.4999))

fieldmat = AFFE_MATERIAU(
                        AFFE=_F(MATER=(mater, ), TOUT='OUI'), MAILLAGE=MAIL, MODELE=model
                    )

myOptions=(
# Pressure KSP
""" -prefix_push fieldsplit_PRES_ """ +
""" -ksp_max_it 5 """ +
""" -ksp_monitor """ +
""" -ksp_converged_reason """ +
""" -ksp_rtol 1.e-3 """ +
""" -ksp_type fgmres """ +
""" -pc_type jacobi """ +
""" -prefix_pop """ +
# Displ KSP
""" -prefix_push fieldsplit_DXDYDZ_ """ +
""" -ksp_monitor """ +
""" -ksp_converged_reason """ +
""" -ksp_rtol 1e-3 """ +
""" -ksp_max_it 20 """ +
""" -ksp_type gmres """ +
""" -pc_type hypre  -pc_hypre_boomeramg_coarsen_type HMIS """ +
""" -prefix_pop """ +
# Global KSP
""" -ksp_monitor """ +
""" -ksp_converged_reason  """ +
""" -ksp_type fgmres """ +
""" -pc_fieldsplit_schur_factorization_type upper """ +
""" -pc_fieldsplit_schur_precondition a11 """ +
""" -pc_fieldsplit_type schur """ +
""" -log_view """
)

BC = AFFE_CHAR_CINE(
                    MECA_IMPO=(
                                _F(DZ=0.0, GROUP_MA=('Zinf', )),
                                _F(DY=0.0, GROUP_MA=('Yinf', 'Ysup')),
                                _F(DX=0.0, GROUP_MA=('Xsup', 'Xinf')),
                                _F(DX=1.0, DZ=0.0, GROUP_MA=('Zsup', )),
                                _F(PRES=0.0, GROUP_NO='N_test',),
                    ),
                    MODELE=model
                )

resnonl = MECA_STATIQUE(
                        CHAM_MATER=fieldmat,
                        EXCIT=_F(CHARGE=BC, ),
                        MODELE=model,
                        SOLVEUR=_F(METHODE='PETSC', PRE_COND='FIELDSPLIT',RESI_RELA=1.e-8,
                                   NOM_CMP=('DX','DY','DZ','PRES'), PARTITION_CMP=(3,1), OPTION_PETSC=myOptions),
                        # SOLVEUR=_F(METHODE='MUMPS',),
                        INFO=2,
                    )


# ajouter TEST_RESU comme petsc04c
TEST_RESU(
       RESU=_F(
       CRITERE='ABSOLU',
       GROUP_NO='N_test',
       NOM_CHAM='DEPL',
       NOM_CMP='DX',
       NUME_ORDRE=1,
       PRECISION=1.e-6,
       REFERENCE='AUTRE_ASTER',
        RESULTAT=resnonl,
        VALE_CALC=-0.19138178062664846,
        VALE_REFE=-0.19138178062664846,
    ))

TEST_RESU(
       RESU=_F(
       CRITERE='ABSOLU',
       GROUP_NO='N_test2',
       NOM_CHAM='DEPL',
       NOM_CMP='DX',
       NUME_ORDRE=1,
       PRECISION=1.e-6,
       REFERENCE='AUTRE_ASTER',
        RESULTAT=resnonl,
        VALE_CALC=0.1920978452775024,
        VALE_REFE=0.1920978452775024,
    ))


# at least it pass here!

test.printSummary()

# if parallel:
#     rank = MPI.COMM_WORLD.Get_rank()
#     resnonl.printMedFile('/tmp/par_%d.resu.med'%rank)
# else:
#     resnonl.printMedFile('/tmp/seq.resu.med')

FIN()
