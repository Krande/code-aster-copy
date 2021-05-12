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

import os
import code_aster
from code_aster.Commands import *
from code_aster.Utilities import haveMPI
test = code_aster.TestCase()

code_aster.init("--test")

nProc = code_aster.getMPINumberOfProcs()
rank = code_aster.getMPIRank()


pMesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR="PTSCOTCH", INFO_MED=1)
pMesh = MODI_MAILLAGE(reuse=pMesh, MAILLAGE=pMesh,
                    ORIE_PEAU=_F(GROUP_MA_PEAU=('HAUT','BAS','DROITE', 'GAUCHE',),),)


pmodel = AFFE_MODELE(MAILLAGE=pMesh,
                    AFFE=(_F(MODELISATION='D_PLAN', PHENOMENE='MECANIQUE', TOUT='OUI',),
                          _F(MODELISATION='D_PLAN_INCO_UPG', PHENOMENE='MECANIQUE', GROUP_MA='F3'),))



MA=DEFI_MATERIAU(ELAS=_F(E=210000e6,
                         NU=0.3))

pMATE=AFFE_MATERIAU(MAILLAGE=pMesh,
                   AFFE=_F(TOUT='OUI',
                           MATER=MA))

pload0 = AFFE_CHAR_MECA(MODELE=pmodel,
                    DDL_IMPO=(_F(GROUP_MA='GAUCHE', DX=-1.0, ),),
                    )

pload1 = AFFE_CHAR_MECA(MODELE=pmodel,
                    LIAISON_GROUP=(_F(GROUP_MA_1='GAUCHE', DDL_1="DX",
                              GROUP_MA_2='DROITE', DDL_2="DX",
                              COEF_MULT_1=1.0, COEF_MULT_2=1.0, COEF_IMPO=0.),),
                      )

pload2 = AFFE_CHAR_MECA(MODELE=pmodel,
                    PRES_REP=(_F(GROUP_MA='HAUT', PRES=1.0,),
                      ),)


pload3 = AFFE_CHAR_CINE(MODELE=pmodel,
                    MECA_IMPO=(_F(GROUP_MA='BAS', DY=0.0,  ),
                      ),)



pRESU=MECA_STATIQUE(MODELE=pmodel,
                        CHAM_MATER=pMATE,
                        EXCIT=(_F(CHARGE=pload0),
                                _F(CHARGE=pload1),
                                _F(CHARGE=pload2), _F(CHARGE=pload3),),
                        INST=1.0,
                        SOLVEUR=_F(METHODE="PETSC",RESI_RELA=1e-9,),
                        )


# os.system("cp fort.%d /home/C00976/tmp/bug_affe/fort.%d"%(130+rank, 130+rank))
# os.system("cp fort.%d /home/C00976/tmp/bug_affe/fort.%d"%(190+rank, 190+rank))
# os.system("cp fort.%d /home/C00976/tmp/bug_affe/fort.%d"%(601+rank, 601+rank))


TEST_RESU(
                RESU=(_F(
                    GROUP_NO=('GN0', ),
                    NOM_CHAM='DEPL',
                    NOM_CMP='DX',
                    INST=1.0,
                    PRECISION=0.00001,
                    REFERENCE='AUTRE_ASTER',
                    RESULTAT=pRESU,
                    CRITERE='ABSOLU',
                    VALE_CALC=-0.33333333333333204,
                    VALE_REFE=-0.33333333333333204,
                ),
            _F(
                    GROUP_NO=('GN0', ),
                    NOM_CHAM='DEPL',
                    NOM_CMP='DY',
                    INST=1.0,
                    REFERENCE='AUTRE_ASTER',
                    RESULTAT=pRESU,
                    ORDRE_GRANDEUR=1e-6,
                    CRITERE='ABSOLU',
                    VALE_CALC=2.6988334396258498E-17,
                    VALE_REFE=0.0,
                ),
            ),)


TEST_RESU(
                RESU=(_F(
                    GROUP_NO=('GN1', ),
                    NOM_CHAM='DEPL',
                    NOM_CMP='DX',
                    INST=1.0,
                    PRECISION=0.00001,
                    REFERENCE='AUTRE_ASTER',
                    RESULTAT=pRESU,
                    CRITERE='RELATIF',
                    VALE_CALC=0.33333333333333515,
                    VALE_REFE=0.33333333333333515,
                ),
            _F(
                    CRITERE=('RELATIF', ),
                    GROUP_NO=('GN1', ),
                    NOM_CHAM='DEPL',
                    NOM_CMP='PRES',
                    INST=1.0,
                    PRECISION=0.00001,
                    REFERENCE='AUTRE_ASTER',
                    RESULTAT=pRESU,
                    VALE_CALC=19999999999.380936,
                    VALE_REFE=19999999999.380936,
                ),
            ),)

if False:
#if haveMPI():
    # pour one mpi
    deeg_1 = [-3, -1, -9, 0, 1, 1, 1, 2, -1, -1, -4, -1, -7, 0, 2, 1, 2, 2, -2, -1, 3, 1, 3, 2, 3, 9, 3, 97, -10, 0, 4, 1, 4, 2, 4, 9, 4, 97, -8, 0, 5, 1, 5, 2, 6, 1, 6, 2, 6, 9, 6, 97, 7, 1, 7, 2, 8, 1, 8, 2, 8, 9, 8, 97, -5, -1, -11, 0, 9, 1, 9, 2, -6, -1, 10, 1, 10, 2, 11, 1, 11, 2, 12, 1, 12, 2, 13, 1, 13, 2, -12, 0, 14, 1, 14, 2, 15, 1, 15, 2, 16, 1, 16, 2, 17, 1, 17, 2, 18, 1, 18, 2]

    nuls_1 = [46, 52, 0, 1, 44, 47, 50, 2, 3, 45, 4, 5, 6, 7, 53, 8, 9, 10, 11, 51, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 48, 54, 24, 25, 49, 26, 27, 28, 29, 30, 31, 32, 33, 55, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43]

    test.assertEqual(sum(deeg_1), 778)
    test.assertEqual(sum(nuls_1), 1540)
    test.assertEqual(len(nuls_1), len(set(nuls_1)))

    # on construit une liste nueq -> nume_noeud (elle doit être indépendante de
    # la découpe des domaines)
    pos = 0
    nueq_nuno_1 = {}
    for nueq in nuls_1:
        nueq_nuno_1[nueq] = [deeg_1[pos], deeg_1[pos+1]]
        pos += 2

    print(nueq_nuno_1, flush=True)



    pnumeddl = NUME_DDL(MODELE=pmodel, CHARGE=(pload0, pload1, pload2))

    pdeeg = pnumeddl.getDEEG()
    pnuls = pnumeddl.getNULS()

    # on vérifie que le numerotation reste toujours la même
    pos = 0
    for nueq in pnuls:
        test.assertSequenceEqual(nueq_nuno_1[nueq], [pdeeg[pos], pdeeg[pos+1]])
        pos += 2


test.printSummary()

FIN()
