# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

# ----------------------------------------------------------------
# 1) ON VEUT FABRIQUER UN CHAM_ELEM (ELGA) DE SIEF_R   : SIGINI
# QUI POURRAIT SERVIR D'ETAT INITIAL POUR UN CALCUL DE MECANIQUE
# DES SOLS :
# ON VEUT :
#     SIZZ = RHO*G*Z
#     SIYY = SIXX = KP*SIZZ
#
# 2) ON VEUT FABRIQUER UN CHAM_NO DE TEMP_R    : TEMP0
# ON VEUT :
#     TEMP = 2*X +3*Y +3*Z +4*INST
#
# 3) ON VEUT FABRIQUER UN cham_elem_VARI_R  (ELGA)
#   pouvant servir de champ de variables internes initiales pour STAT_NON_LINE
# ----------------------------------------------------------------

from code_aster.Commands import *
from code_aster import CA
from code_aster.Helpers.debugging import DebugChrono


DEBUT(CODE="OUI", IGNORE_ALARM=("MECANONLINE5_37", "MECANONLINE2_37"), DEBUG=_F(SDVERI="OUI"))

# IGNORE_ALARME: ARRET='NON' AUTORISE POUR VALIDATION CREA_RESU

MA = LIRE_MAILLAGE(FORMAT="ASTER")

MO = AFFE_MODELE(MAILLAGE=MA, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))


# ===========================================================
# 1) Construction de sigm_init en passant par un champ simple
# ===========================================================

with DebugChrono.measure("sigma"):
    RHO = 1000.0
    G = 10.0
    KP = 3.0

    # Récupération des coordonnées aux points de Gauss
    coor_elga = CALC_CHAM_ELEM(MODELE=MO, OPTION="COOR_ELGA").toSimpleFieldOnCells()
    valz = coor_elga.Z

    # Initialisation pour avoir la description de SIEF_ELGA
    stress = CA.FieldOnCellsReal(MO, "ELGA", "SIEF_R").toSimpleFieldOnCells()
    sixx = stress.SIXX

    # COOR_ELGA est défini sur toutes les mailles du modèle, pas SIEF_ELGA
    valz = valz.onSupportOf(sixx)

    # Calcul des composantes
    vz = RHO * G * valz
    vx = KP * vz
    vy = vx

    # Affectation des composantes
    stress.setComponentValues("SIXX", vx)
    stress.setComponentValues("SIYY", vy)
    stress.setComponentValues("SIZZ", vz)
    del vx, vy, vz, sixx

    # Conversion en FieldOnCells
    fed = MO.getFiniteElementDescriptor()
    sigm_init = stress.toFieldOnCells(fed)

TEST_RESU(
    CHAM_ELEM=(
        _F(
            REFERENCE="ANALYTIQUE",
            POINT=2,
            NOM_CMP="SIZZ",
            GROUP_MA="CUBE",
            CHAM_GD=sigm_init,
            VALE_CALC=5.77350269e03,
            VALE_REFE=5773.5,
        ),
        _F(
            REFERENCE="ANALYTIQUE",
            POINT=2,
            NOM_CMP="SIYY",
            GROUP_MA="CUBE",
            CHAM_GD=sigm_init,
            VALE_CALC=1.73205081e04,
            VALE_REFE=1.73205e4,
        ),
    )
)

# =======================================
# 2) Construction du champ de température
# =======================================

with DebugChrono.measure("temp"):
    # Récupération des coordonnées aux noeuds
    coor = MA.getCoordinatesAsSimpleFieldOnNodes()
    valx = coor.X
    valy = coor.Y
    valz = coor.Z

    inst = 1.0
    temp0 = 2.0 * valx + 3.0 * valy + 4.0 * valz + 5.0 * inst

    # Création du champ simple
    tempf = CA.SimpleFieldOnNodesReal(MA, "TEMP_R", ["TEMP"], prol_zero=True)

    # Affectation de la composante TEMP
    tempf.setComponentValues("TEMP", temp0)

    # Conversion en FieldOnNodes
    temp_init = tempf.toFieldOnNodes()
    del valx, valy, valz, temp0

TEST_RESU(
    CHAM_NO=_F(
        GROUP_NO="N7",
        REFERENCE="ANALYTIQUE",
        NOM_CMP="TEMP",
        CHAM_GD=temp_init,
        VALE_CALC=14.000000000,
        VALE_REFE=14.0,
    )
)

# ===============================================================================================
# 3) Creation d'un champ de variables internes initiales non nulles  :
# -------------------------------------------------------------------
#
#  On veut que :
#     STAT_NON_LINE :
#         COMPORTEMENT=(_F( GROUP_MA='MASSIF',RELATION = 'CJS'),
#                    _F( GROUP_MA='BETON', RELATION = 'ENDO_ISOT_BETON'),),
#
#     pour la relation de comportement CJS (16 variables internes), on veut affecter :
#           - V1 = 1.0   et  V9 = 9.0
#
#     pour la relation de comportement ENDO_ENDO_ISOT_BETON (2 variables internes), on veut affecter :
#           - V2 = 2.0
#
# ===============================================================================================


# 3.1  calcul non lineaire bidon pour avoir un modele du champ de variables internes
# ----------------------------------------------------------------------------------
BETON = DEFI_MATERIAU(ELAS=_F(E=20000.0, NU=0.0), BETON_ECRO_LINE=_F(SYT=6.0, D_SIGM_EPSI=-10000.0))


MASSIF = DEFI_MATERIAU(
    ELAS=_F(E=35.0e3, NU=0.15),
    CJS=_F(
        BETA_CJS=-0.55,
        GAMMA_CJS=0.82,
        PA=-100.0,
        RM=0.289,
        N_CJS=0.6,
        KP=25.5e3,
        RC=0.265,
        A_CJS=0.25,
    ),
)

CHMAT = AFFE_MATERIAU(
    MAILLAGE=MA, AFFE=(_F(GROUP_MA="MASSIF", MATER=MASSIF), _F(GROUP_MA="BETON", MATER=BETON))
)

# 3.2 Définition du champ de variables internes
# ---------------------------------------------

with DebugChrono.measure("create behaviour"):
    # On ne peut pas créer un champ de variables internes sans le comportement.
    # Définition du problème et affectation des comportements
    phys_pb = CA.PhysicalProblem(MO, CHMAT)
    phys_pb.computeBehaviourProperty(
        COMPORTEMENT=(
            _F(GROUP_MA="MASSIF", RELATION="CJS"),
            _F(GROUP_MA="BETON", RELATION="ENDO_ISOT_BETON"),
        )
    )

with DebugChrono.measure("vari.init"):
    vari_elga = CA.FieldOnCellsReal(MO, "ELGA", "VARI_R", phys_pb.getBehaviourProperty())
    # Mise à zero de toutes les composantes
    vari_elga.setValues(0.0)

with DebugChrono.measure("build vari_elga"):

    with DebugChrono.measure("vari.toSimple"):
        # Conversion en champ simple
        vari_elga = vari_elga.toSimpleFieldOnCells()

    with DebugChrono.measure("vari.oper"):
        # Liste des mailles pour chaque zone d'intérêt
        beton = MA.getCells("BETON")
        massif = MA.getCells("MASSIF")

        # On affecte la même valeur sur la zone considérée
        v1 = vari_elga.V1
        v1.restrict(massif)
        v1.values = 1.0

        v9 = vari_elga.V9
        v9.restrict(massif)
        v9.values = 9.0

        v2 = vari_elga.V2
        v2.restrict(beton)
        v2.values = 2.0

    with DebugChrono.measure("vari.assign"):
        vari_elga.setComponentValues("V1", v1)
        vari_elga.setComponentValues("V2", v2)
        vari_elga.setComponentValues("V9", v9)

    with DebugChrono.measure("vari.toField"):
        # Conversion en FieldOnCells
        fed = MO.getFiniteElementDescriptor()
        vari_init = vari_elga.toFieldOnCells(fed)

with DebugChrono.measure("example"):
    # À titre d'exemple, on montre comment affecter une valeur sur "BETON" et une autre sur "MASSIF"
    v2 = vari_elga.V2
    # On crée deux objets Component indépendants car restrict modifie l'objet en place
    v2b = v2.copy()
    v2b.restrict(beton)
    v2b.values = 2.0

    v2m = v2.copy()
    v2m.restrict(massif)
    v2m.values = 1.0e-6

    # v2b et v2m ne sont pas définies sur les mêmes mailles,
    # il faut donc les transposer sur le même support, celui de v2
    v2tot = v2b.onSupportOf(v2) + v2m.onSupportOf(v2)
    # On affecterait 'v2tot'
    # vari_elga.setComponentValues("V2", v2tot)

    VAIN2 = CREA_CHAMP(
        OPERATION="AFFE",
        TYPE_CHAM="CART_NEUT_R",
        MODELE=MO,
        INFO=1,
        AFFE=(
            _F(GROUP_MA="BETON", NOM_CMP=("X2",), VALE=(2.0,)),
            _F(GROUP_MA="MASSIF", NOM_CMP=("X1", "X3"), VALE=(1.0, 9.0)),
        ),
    )


# 3.3 test des valeurs obtenues
# -----------------------------
TEST_RESU(
    CHAM_ELEM=_F(
        REFERENCE="ANALYTIQUE",
        POINT=8,
        NOM_CMP="V2",
        GROUP_MA="CUBE",
        CHAM_GD=vari_init,
        VALE_CALC=2.000000000,
        VALE_REFE=2.0,
    )
)

TEST_RESU(
    CHAM_ELEM=_F(
        REFERENCE="ANALYTIQUE",
        POINT=3,
        NOM_CMP="V9",
        GROUP_MA="MASSIF",
        CHAM_GD=vari_init,
        VALE_CALC=9.000000000,
        VALE_REFE=9.0,
    )
)

TEST_RESU(
    CHAM_ELEM=_F(
        CRITERE="ABSOLU",
        REFERENCE="ANALYTIQUE",
        POINT=3,
        NOM_CMP="V2",
        GROUP_MA="MASSIF",
        CHAM_GD=vari_init,
        VALE_CALC=0.00000000e00,
        VALE_REFE=0.0,
    )
)

FIN()
