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

DEBUT(CODE=_F(NIV_PUB_WEB='INTERNET',),
      DEBUG=_F(SDVERI='OUI',),
      INFO=1,)

test = code_aster.TestCase()

mesh=LIRE_MAILLAGE(FORMAT='MED', UNITE=20, PARTITIONNEUR='PTSCOTCH')

model=AFFE_MODELE(MAILLAGE=mesh,
                    AFFE=_F(TOUT='OUI',
                         PHENOMENE='MECANIQUE',
                         MODELISATION='3D',),)

# Very high elasticity limit to simulate elasticity
acier=DEFI_MATERIAU(ELAS=_F(E=200000.,
                            NU=0.3,),
                    ECRO_LINE=_F(D_SIGM_EPSI=2000.,
                                 SY=200000.,),)

mater=AFFE_MATERIAU(MAILLAGE=mesh,
                   AFFE=_F(TOUT='OUI',
                           MATER=acier,),)


encast=AFFE_CHAR_MECA(MODELE=model,
                     DDL_IMPO=(_F(GROUP_MA='BAS',
                                   DX=0, DY=0.0, DZ=0.0,
                                    ),),)

depl=AFFE_CHAR_MECA(MODELE=model,
                     DDL_IMPO=(_F(GROUP_MA='HAUT',
                                    DZ=1.0,
                                    ),),)

LIST=DEFI_LIST_REEL(DEBUT=0.0,
                    INTERVALLE=_F(JUSQU_A=1.0,
                                  NOMBRE=2,),)

RAMPE=DEFI_FONCTION(NOM_PARA='INST',
                        VALE=(0., 0.,
                            1000., 1000.))

# STAT_NON_LINE DE REFERENCE
SOLUT=STAT_NON_LINE(MODELE=model,
                    CHAM_MATER=mater,
                    EXCIT=(_F(CHARGE=encast,FONC_MULT=RAMPE),
                           _F(CHARGE=depl,FONC_MULT=RAMPE),),
                    COMPORTEMENT=_F(RELATION='VMIS_ISOT_LINE',),
                    CONVERGENCE=_F(RESI_GLOB_MAXI=1e-8,),
                    NEWTON=_F(REAC_INCR=1,
                              PREDICTION='ELASTIQUE',
                              MATRICE='TANGENTE',
                              REAC_ITER=1,),
                    INCREMENT=_F(LIST_INST=LIST,),
                    INFO=1)

print("Field in original SNL:", flush=True)
SOLUT.printListOfFields()

nbRank = SOLUT.getNumberOfRanks()


# New result
SOLUN = code_aster.NonLinearResult()

SOLUN.allocate(nbRank)

for rank in range(nbRank):
    SOLUN.setModel( SOLUT.getModel(rank), rank )
    SOLUN.setMaterialField( SOLUT.getMaterialField(rank), rank )
    SOLUN.setField( SOLUT.getFieldOnNodesReal("DEPL", rank), "DEPL", rank )
    SOLUN.setField( SOLUT.getFieldOnCellsReal("SIEF_ELGA", rank), "SIEF_ELGA", rank )
    SOLUN.setField( SOLUT.getFieldOnCellsReal("VARI_ELGA", rank), "VARI_ELGA", rank )
    compor = SOLUT.getConstantFieldOnCellsChar16("COMPORTEMENT", rank)
    SOLUN.setField( compor, "COMPORTEMENT", rank )
    SOLUN.setTimeValue( SOLUT.getTimeValue(rank), rank)

list_field = []
list_field += SOLUN.getFieldsOnNodesNames()
list_field += SOLUN.getFieldsOnCellsNames()
list_field += SOLUN.getConstantFieldsOnCellsNames()

test.assertSequenceEqual( sorted(list_field), sorted(["DEPL", "VARI_ELGA", "SIEF_ELGA", "COMPORTEMENT"]))


#=======================================================
#             IMPR_RESU et LIRE_RESU
#=======================================================

IMPR_RESU(
        FORMAT='MED',UNITE=81,
        RESU=_F(RESULTAT=SOLUT,),
    )

IMPR_RESU(
        FORMAT='MED',UNITE=80,
        RESU=_F(RESULTAT=SOLUN, NOM_RESU_MED = "SOLUN",),
    )

SOLUT = LIRE_RESU(MODELE=model,
                 FORMAT='MED',
                 UNITE=81,
                 TYPE_RESU='EVOL_NOLI',
                 COMPORTEMENT=_F(RELATION='VMIS_ISOT_LINE',),
                 CHAM_MATER=mater,
                 EXCIT=(_F(CHARGE=encast,FONC_MULT=RAMPE),
                           _F(CHARGE=depl,FONC_MULT=RAMPE),),
                 FORMAT_MED=(_F(NOM_RESU="SOLUT", NOM_CHAM='DEPL'),
                             _F(NOM_RESU="SOLUT", NOM_CHAM='SIEF_ELGA'),
                             _F(NOM_RESU="SOLUT", NOM_CHAM='VARI_ELGA')),
                 TOUT_ORDRE="OUI")

SOLUR = LIRE_RESU(MODELE=model,
                 FORMAT='MED',
                 UNITE=80,
                 TYPE_RESU='EVOL_NOLI',
                 COMPORTEMENT=_F(RELATION='VMIS_ISOT_LINE',),
                 CHAM_MATER=mater,
                 EXCIT=(_F(CHARGE=encast,FONC_MULT=RAMPE),
                           _F(CHARGE=depl,FONC_MULT=RAMPE),),
                 FORMAT_MED=(_F(NOM_RESU="SOLUN", NOM_CHAM='DEPL'),
                             _F(NOM_RESU="SOLUN", NOM_CHAM='SIEF_ELGA'),
                             _F(NOM_RESU="SOLUN", NOM_CHAM='VARI_ELGA')),
                 TOUT_ORDRE="OUI")

#=========================================================
#            REALISATION DES TESTS
#=========================================================

test.assertEqual(SOLUT.getNumberOfRanks(), SOLUR.getNumberOfRanks())
test.assertSequenceEqual(SOLUT.getRanks(), [0, 1, 2])
test.assertSequenceEqual(SOLUR.getRanks(), [0, 1, 2])

# ON EXTRAIT LES CHAMPS A TESTER au dernier instant

for rank in range(SOLUT.getNumberOfRanks()):

    DEPL_REF = SOLUT.getFieldOnNodesReal("DEPL", rank)
    SIGMA_REF= SOLUT.getFieldOnCellsReal("SIEF_ELGA", rank)
    VARI_REF=SOLUT.getFieldOnCellsReal("VARI_ELGA", rank)

    DEPL= SOLUR.getFieldOnNodesReal("DEPL", rank)
    SIGMA=SOLUR.getFieldOnCellsReal("SIEF_ELGA", rank)
    VARI=SOLUR.getFieldOnCellsReal("VARI_ELGA", rank)

    DIF_DEPL = DEPL_REF - DEPL

    # DIF_SIG = SIG_REF - SIGM
    DIF_SIG=CREA_CHAMP(OPERATION='ASSE',
                    MODELE=model,
                    TYPE_CHAM='ELGA_SIEF_R',
                    ASSE=(_F(CHAM_GD=SIGMA_REF,
                                TOUT='OUI',
                                CUMUL='OUI',
                                COEF_R=1.),
                            _F(CHAM_GD=SIGMA,
                                TOUT='OUI',
                                CUMUL='OUI',
                                COEF_R=-1.),
                            ))


    # DIF_VAR = VAR_REF - VARI
    DIF_VAR=CREA_CHAMP(OPERATION='ASSE',
                    MODELE=model,
                    TYPE_CHAM='ELGA_VARI_R',
                    ASSE=(_F(CHAM_GD=VARI_REF,
                                TOUT='OUI',
                                CUMUL='OUI',
                                COEF_R=1.),
                            _F(CHAM_GD=VARI,
                                TOUT='OUI',
                                CUMUL='OUI',
                                COEF_R=-1.),
                            ))

    TEST_RESU(CHAM_ELEM=(_F(CRITERE='ABSOLU',
                            REFERENCE='AUTRE_ASTER',
                            PRECISION=1.E-08,
                            TYPE_TEST='MIN',
                            CHAM_GD=DIF_SIG,
                            VALE_CALC= 1.5916157281026244E-12,
                            VALE_REFE=0.0,
                            VALE_ABS='OUI',),
                        _F(CRITERE='ABSOLU',
                            REFERENCE='AUTRE_ASTER',
                            PRECISION=1.E-08,
                            TYPE_TEST='MAX',
                            CHAM_GD=DIF_SIG,
                            VALE_CALC= 1.1368683772161603E-12,
                            VALE_REFE=0.0,
                            VALE_ABS='OUI',),
                        _F(CRITERE='ABSOLU',
                            REFERENCE='AUTRE_ASTER',
                            ORDRE_GRANDEUR=1.E-08,
                            TYPE_TEST='MIN',
                            CHAM_GD=DIF_VAR,
                            VALE_CALC= 0.0,
                            VALE_REFE=0.0,
                            VALE_ABS='OUI',),
                        _F(CRITERE='ABSOLU',
                            REFERENCE='AUTRE_ASTER',
                            ORDRE_GRANDEUR=1.0e-8,
                            TYPE_TEST='MAX',
                            CHAM_GD=DIF_VAR,
                            VALE_CALC= 0.0,
                            VALE_REFE=0.0,
                            VALE_ABS='OUI',),
                        ),
            )


    TEST_RESU(CHAM_NO=(_F(CRITERE='ABSOLU',
                            REFERENCE='AUTRE_ASTER',
                            PRECISION=1.E-08,
                            TYPE_TEST='MIN',
                            CHAM_GD=DIF_DEPL,
                            VALE_CALC= 6.984919309616089E-10,
                            VALE_REFE=0.0,
                            VALE_ABS='OUI',),
                        _F(CRITERE='ABSOLU',
                            REFERENCE='AUTRE_ASTER',
                            PRECISION=1.E-08,
                            TYPE_TEST='MAX',
                            CHAM_GD=DIF_DEPL,
                            VALE_CALC= 9.313225746154785E-10,
                            VALE_REFE=0.0,
                            VALE_ABS='OUI',),
            ))

FIN()
