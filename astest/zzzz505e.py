# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

mesh = LIRE_MAILLAGE(FORMAT='MED', UNITE=20)

model = AFFE_MODELE(MAILLAGE=mesh,
                    AFFE=_F(TOUT='OUI',
                            PHENOMENE='MECANIQUE',
                            MODELISATION='3D',),)

# Very high elasticity limit to simulate elasticity
acier = DEFI_MATERIAU(ELAS=_F(E=200000.,
                              NU=0.3,),
                      ECRO_LINE=_F(D_SIGM_EPSI=2000.,
                                   SY=200000.,),)

mater = AFFE_MATERIAU(MAILLAGE=mesh,
                      AFFE=_F(TOUT='OUI',
                              MATER=acier,),)


encast = AFFE_CHAR_MECA(MODELE=model,
                        DDL_IMPO=(_F(GROUP_MA='BAS',
                                     DX=0, DY=0.0, DZ=0.0,
                                     ),),)

depl = AFFE_CHAR_MECA(MODELE=model,
                      DDL_IMPO=(_F(GROUP_MA='HAUT',
                                   DZ=1.0,
                                   ),),)

LIST = DEFI_LIST_REEL(DEBUT=0.0,
                      INTERVALLE=_F(JUSQU_A=1.0,
                                    NOMBRE=2,),)

RAMPE = DEFI_FONCTION(NOM_PARA='INST',
                      VALE=(0., 0.,
                            1000., 1000.))

# STAT_NON_LINE DE REFERENCE
SOLUT = STAT_NON_LINE(MODELE=model,
                      CHAM_MATER=mater,
                      EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE),
                             _F(CHARGE=depl, FONC_MULT=RAMPE),),
                      COMPORTEMENT=_F(RELATION='VMIS_ISOT_LINE',),
                      CONVERGENCE=_F(RESI_GLOB_MAXI=1e-8,),
                      NEWTON=_F(REAC_INCR=1,
                                PREDICTION='ELASTIQUE',
                                MATRICE='TANGENTE',
                                REAC_ITER=1,),
                      INCREMENT=_F(LIST_INST=LIST,),
                      INFO=1)

nbRank = SOLUT.getNumberOfRanks()

# New result
SOLUN = code_aster.NonLinearResult()

SOLUN.allocate(nbRank)

for rank in (0, 1, 2, 2, 1):
    print("Rank: ", rank, flush=True)
    SOLUN.setModel(SOLUT.getModel(rank), rank)
    SOLUN.setMaterialField(SOLUT.getMaterialField(rank), rank)
    SOLUN.setField(SOLUT.getField("DEPL", rank), "DEPL", rank)
    SOLUN.setField(SOLUT.getField("SIEF_ELGA", rank), "SIEF_ELGA", rank)
    SOLUN.setField(SOLUT.getField("VARI_ELGA", rank), "VARI_ELGA", rank)
    compor = SOLUT.getField("COMPORTEMENT", rank)
    SOLUN.setField(compor, "COMPORTEMENT", rank)
    SOLUN.setTimeValue(SOLUT.getTimeValue(rank), rank)

    depl = CREA_CHAMP(OPERATION="EXTR",
                      TYPE_CHAM="NOEU_DEPL_R",
                      NOM_CHAM="DEPL",
                      RESULTAT=SOLUN,
                      NUME_ORDRE=rank)

    sief = CREA_CHAMP(OPERATION="EXTR",
                      TYPE_CHAM="ELGA_SIEF_R",
                      NOM_CHAM="SIEF_ELGA",
                      RESULTAT=SOLUN,
                      NUME_ORDRE=rank)

# =========================================================
#            REALISATION DES TESTS
# =========================================================

test.assertEqual(SOLUT.getNumberOfRanks(), SOLUN.getNumberOfRanks())
test.assertSequenceEqual(SOLUT.getRanks(), [0, 1, 2])
test.assertSequenceEqual(SOLUN.getRanks(), [0, 1, 2])

# ON EXTRAIT LES CHAMPS A TESTER au dernier instant

for rank in range(SOLUT.getNumberOfRanks()):

    DEPL_REF = SOLUT.getField("DEPL", rank)
    SIGMA_REF = SOLUT.getField("SIEF_ELGA", rank)
    VARI_REF = SOLUT.getField("VARI_ELGA", rank)

    DEPL = SOLUN.getField("DEPL", rank)
    SIGMA = SOLUN.getField("SIEF_ELGA", rank)
    VARI = SOLUN.getField("VARI_ELGA", rank)

    DIF_DEPL = DEPL_REF - DEPL

    # DIF_SIG = SIG_REF - SIGM
    DIF_SIG = CREA_CHAMP(OPERATION='ASSE',
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
    DIF_VAR = CREA_CHAMP(OPERATION='ASSE',
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
                            VALE_CALC=1.5916157281026244E-12,
                            VALE_REFE=0.0,
                            VALE_ABS='OUI',),
                         _F(CRITERE='ABSOLU',
                            REFERENCE='AUTRE_ASTER',
                            PRECISION=1.E-08,
                            TYPE_TEST='MAX',
                            CHAM_GD=DIF_SIG,
                            VALE_CALC=1.1368683772161603E-12,
                            VALE_REFE=0.0,
                            VALE_ABS='OUI',),
                         _F(CRITERE='ABSOLU',
                            REFERENCE='AUTRE_ASTER',
                            ORDRE_GRANDEUR=1.E-08,
                            TYPE_TEST='MIN',
                            CHAM_GD=DIF_VAR,
                            VALE_CALC=0.0,
                            VALE_REFE=0.0,
                            VALE_ABS='OUI',),
                         _F(CRITERE='ABSOLU',
                            REFERENCE='AUTRE_ASTER',
                            ORDRE_GRANDEUR=1.0e-8,
                            TYPE_TEST='MAX',
                            CHAM_GD=DIF_VAR,
                            VALE_CALC=0.0,
                            VALE_REFE=0.0,
                            VALE_ABS='OUI',),
                         ),
              )

    TEST_RESU(CHAM_NO=(_F(CRITERE='ABSOLU',
                          REFERENCE='AUTRE_ASTER',
                          PRECISION=1.E-08,
                          TYPE_TEST='MIN',
                          CHAM_GD=DIF_DEPL,
                          VALE_CALC=6.984919309616089E-10,
                          VALE_REFE=0.0,
                          VALE_ABS='OUI',),
                       _F(CRITERE='ABSOLU',
                          REFERENCE='AUTRE_ASTER',
                          PRECISION=1.E-08,
                          TYPE_TEST='MAX',
                          CHAM_GD=DIF_DEPL,
                          VALE_CALC=9.313225746154785E-10,
                          VALE_REFE=0.0,
                          VALE_ABS='OUI',),
                       ))

# =========================================================
#            TEST RESIZE METHOD
# =========================================================

nbRank = SOLUN.getNumberOfRanks()
SOLUN.resize(nbRank + 1)
SOLUN.setModel(SOLUN.getModel(nbRank-1), nbRank)
SOLUN.setMaterialField(SOLUN.getMaterialField(nbRank-1), nbRank)
SOLUN.setField(SOLUN.getFieldOnNodesReal("DEPL", nbRank-1), "DEPL", nbRank)

# CHECK THAT the fields are properly copied
test.assertEqual(SOLUN.getModel(nbRank-1), SOLUN.getModel(nbRank))
test.assertEqual(SOLUN.getMaterialField(nbRank-1),
                 SOLUN.getMaterialField(nbRank-1))
test.assertEqual(SOLUN.getFieldOnNodesReal("DEPL", nbRank-1).getValues(),
                 SOLUN.getFieldOnNodesReal("DEPL", nbRank).getValues())


# =========================================================
#            TEST NORM METHODS
# =========================================================

sig = SOLUN.getFieldOnCellsReal("SIEF_ELGA", nbRank - 1)
test.assertAlmostEqual(sig.norm("NORM_2"), 24162.483694014656, 8)
test.assertAlmostEqual(sig.norm("NORM_1"), 640625.2493377978, 8)
test.assertAlmostEqual(sig.norm("NORM_INFINITY"), 1317.268981963458, 8)


# Fix 31780
depl = SOLUT.getField("DEPL", 0)
depl.printMedFile("depl.med")

chvga = SOLUT.getField("VARI_ELGA", 0)
chvga.printMedFile("vari_elga.med")


chmaxsig = CREA_CHAMP(CRITERE='RELATIF',
                      INFO=1,
                      NUME_ORDRE=0,
                      NOM_CHAM='SIEF_ELGA',
                      OPERATION='EXTR',
                      PRECISION=1e-06,
                      PROL_ZERO='NON',
                      RESULTAT=SOLUT,
                      TYPE_CHAM='ELGA_SIEF_R',
                      TYPE_MAXI='MAXI',
                      TYPE_RESU='VALE')

print("CHMAX1: ", chmaxsig.getDescription())

chprint = CREA_RESU(AFFE=_F(CHAM_GD=chmaxsig,
                            CRITERE='RELATIF',
                            INST=0.0,
                            PRECISION=0.0),
                    NOM_CHAM='SIEF_ELGA',
                    OPERATION='AFFE',
                    TYPE_RESU='EVOL_NOLI',
                    VERI_VARI='OUI')

chmaxsig2 = CREA_CHAMP(CRITERE='RELATIF',
                      INFO=1,
                      INST=0.0,
                      NOM_CHAM='SIEF_ELGA',
                      OPERATION='EXTR',
                      PRECISION=1e-06,
                      PROL_ZERO='NON',
                      RESULTAT=chprint,
                      TYPE_CHAM='ELGA_SIEF_R',
                      TYPE_MAXI='MAXI',
                      TYPE_RESU='VALE')

print("CHMAX2: ", chmaxsig2.getDescription())


chprint2 = CREA_RESU(AFFE=_F(CHAM_GD=chmaxsig2,
                            CRITERE='RELATIF',
                            INST=0.0,
                            PRECISION=0.0),
                    NOM_CHAM='SIEF_ELGA',
                    OPERATION='AFFE',
                    TYPE_RESU='EVOL_NOLI',
                    VERI_VARI='OUI')

FIN()
