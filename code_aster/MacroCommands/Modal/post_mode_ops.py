# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

# person_in_charge: anaelle.torre at edf.fr

import sys
import numpy as np

import aster
from ...Cata.Syntax import _F
from ...Supervis import CO
from ...Messages import UTMESS
from ...CodeCommands import CREA_TABLE, POST_ELEM, CREA_CHAMP, PROD_MATR_CHAM, PROJ_BASE


def getGravityCenter(model, cara, mater, group_ma):
    """
    Calculation of the gravity center

    Arguments:
        model (model): model associated to the mode_meca object
        cara (cara_elem): Elementary caracteristics associated to the mode_meca object
        mater (cham_mater): Material field associated to the mode_meca object
        group_ma (grma) : groups of elements
    Returns:
        float(3) : coordinates of gravity center
    """

    t_list = isinstance(group_ma, list)
    t_tuple = isinstance(group_ma, tuple)
    if (not t_list) & (not t_tuple):
        group_ma = (group_ma,)

    if cara is not None:
        Mass_Iner = POST_ELEM(
            MODELE=model, CARA_ELEM=cara, CHAM_MATER=mater, MASS_INER=_F(GROUP_MA=group_ma)
        )
    else:
        Mass_Iner = POST_ELEM(
            MODELE=model, CARA_ELEM=cara, CHAM_MATER=mater, MASS_INER=_F(GROUP_MA=group_ma)
        )

    Mass_Tot = Mass_Iner.EXTR_TABLE().values()

    if len(group_ma) > 1:
        for i1, lieu in enumerate(Mass_Tot["LIEU"]):
            if lieu == "UNION_GROUP_MA":
                ind = i1
                break
    else:
        ind = 0

    # Masse = Mass_Tot['MASSE'][ind]
    CG = [Mass_Tot["CDG_%s" % d][ind] for d in ["X", "Y", "Z"]]
    # Iner = {}
    # for k in ['IX_G', 'IY_G', 'IZ_G', 'IXY_G', 'IXZ_G', 'IYZ_G'] :
    #     Iner.setdefault(k, Mass_Tot[k][ind])

    return CG


def getNodesFromGroups(mesh, group_ma):
    """
    Get list of nodes associated to the groupes of elements

    Arguments:
        mesh (mesh): mesh associated to the mode_meca object
        group_ma (grma) : groups of elements
    Returns:
        Nodes (list of int) : list of nodes
    """
    Nodes = []
    for nam in group_ma:
        Nodes += mesh.getNodesFromCells(group_name=nam)
    return Nodes


def getRigidBodyModes(Nume, group_ma, returnType="CHAM", CG=None, cara=None, mater=None):
    """

    Build the 6 rigid body modes around the point CG

    Arguments:
        Nume : NUME_DDL
        group_ma (grma) : groups of elements
    returnType : 'CHAM' pour récuperer la liste des champs aux noeuds
                 'RESU' pour recuperer un concept MODE_MECA avec les 6 modes
        CG : coordinates of the rotation center
        cara (cara_elem): Elementary caracteristics associated to the mode_meca object
        mater (cham_mater): Material field associated to the mode_meca object

    Returns:
        Cham : list of 6 FieldsOnNodes
    """
    model = Nume.getModel()
    mesh = model.getMesh()

    if CG is None:
        CG = getGravityCenter(model, cara, mater, group_ma)

    # nodes coordinates
    Nb_no = mesh.getNumberOfNodes()
    Coord = np.asarray(mesh.getCoordinates().getValues())
    Coord = Coord.reshape((Nb_no, 3))

    XG = CG[0]
    YG = CG[1]
    ZG = CG[2]

    # Build rigid modes on all nodes and 6 componants (even if componant doesn't existe)
    ModRig = np.zeros((6 * Nb_no, 6))
    ModRig[::6, 0] += 1.0  # -- Trans. X

    ModRig[1::6, 1] += 1.0  # -- Trans. Y
    ModRig[2::6, 2] += 1.0  # -- Trans. Z
    # -- Rot. X
    ModRig[3::6, 3] += 1.0
    ModRig[1::6, 3] -= Coord[:, 2] - ZG
    ModRig[2::6, 3] += Coord[:, 1] - YG
    # -- Rot. Y
    ModRig[4::6, 4] += 1.0
    ModRig[2::6, 4] -= Coord[:, 0] - XG
    ModRig[0::6, 4] += Coord[:, 2] - ZG
    # -- Rot. Z
    ModRig[5::6, 5] += 1.0
    ModRig[1::6, 5] += Coord[:, 0] - XG
    ModRig[0::6, 5] -= Coord[:, 1] - YG

    # Filter rigid modes on groups

    NoGrMa = getNodesFromGroups(mesh, group_ma)
    No = mesh.getNodes()

    Eqn = Nume.getEquationNumbering()
    if isinstance(Eqn, tuple) | isinstance(Eqn, list):
        Equ = Eqn[0]

    Dic = Eqn.getDOFFromNodeAndComponentId()
    # {(NodeId,CompId) : eqId}
    # First DOF Id is 0
    # Fisrt NodeId is 0
    # Fisrt Componant Id is 1

    # associate pseudo DOF (100*no + icomp) to Mogrig Numbering
    AllDDL = {100 * no + j1 + 1: 6 * i1 + j1 for i1, no in enumerate(No) for j1 in range(6)}
    # list of pseudo DOF (100*no + icomp) of the groups
    GrpDDL = [100 * no + j1 + 1 for i1, no in enumerate(NoGrMa) for j1 in range(6)]
    # associate pseudo DOF to real DOF Id (from EquationNumbering)
    MdlDDL = {100 * k[0] + k[1]: Dic[k] for k in Dic.keys()}

    indgrpDDL = []
    indallDDL = []

    # build the list of select groups dof Id : indgrpDDL
    # build the list of select groups Mogrig Id : indallDDL
    for ddl in GrpDDL:
        # try :
        #     #-- indice du DDL du groupe de maille dans le modele
        #     indgrpDDL.append( MdlDDL[ddl] )
        #     #-- indice du DDL du groupe de maille dans les champs rigides "6 DDLs"
        #     indallDDL.append(AllDDL[ddl])
        # except :
        #     pass
        if ddl in MdlDDL:
            indgrpDDL.append(MdlDDL[ddl])
            indallDDL.append(AllDDL[ddl])

    NbDDL = Nume.getNumberOfDOFs()
    # list of dofId not on the groups
    indZero = list(set(np.arange(NbDDL)).difference(indgrpDDL))
    indZero.sort()

    ChamNoRig = np.zeros((NbDDL, 6))
    ChamNoRig[indgrpDDL, :] = ModRig[indallDDL, :]

    Cham = [None] * 6
    List_F = []
    for i1 in range(6):
        # create a fieldOnNodes
        Cham[i1] = CREA_CHAMP(
            MODELE=model,
            NUME_DDL=Nume,
            PROL_ZERO="OUI",
            OPERATION="AFFE",
            TYPE_CHAM="NOEU_DEPL_R",
            AFFE=(_F(TOUT="OUI", NOM_CMP="DX", VALE=0.0),),
        )
        # add values
        Cham[i1].setValues(ChamNoRig[:, i1])
        List_F.append(
            {
                "CHAM_GD": Cham[i1],
                "MODELE": model,
                "NOM_CHAM": "DEPL",
                "NUME_MODE": i1 + 1,
                "FREQ": 0.0,
            }
        )

    if returnType == "CHAM":
        return Cham, indZero
    # elif returnType == 'RESU' :
    #     rigidBodyResu = CREA_RESU( OPERATION='AFFE',
    #                         TYPE_RESU='MODE_MECA',
    #                         AFFE=List_F,
    #                         )

    #     return rigidBodyResu


def getModalMassIner(MODES, rigidBodyChams, indZero, MASSE):
    """

    Calcul des masses et inerties modales 'locales' des modes MODES par rapport à la matrice de masse MASSE,
      autour du point CG

    MODES : MODE_MECA

    Returns:
        table

    """

    # --
    # --  Calculs des produits Phi^T * M * Ud, où Ud sont les modes rigides (R5.01.03)
    # --

    # MASSE = MODES.getMassMatrix()
    nume_ddl = MODES.getDOFNumbering()
    model = nume_ddl.getModel()

    nbRigidModes = 6
    # List_F = []
    # Product M * Ud
    M_Ud = [None] * nbRigidModes
    for n in range(nbRigidModes):
        M_Ud[n] = PROD_MATR_CHAM(MATR_ASSE=MASSE, CHAM_NO=rigidBodyChams[n])
        # List_F.append({'CHAM_GD' : rigidBodyChams[n],
        #                 'MODELE' : model,
        #                 'NOM_CHAM':'DEPL',
        #                 'NUME_MODE' : n+1,
        #                 'FREQ' : 0.0})

    # y a t il un intéret à créer un résultat
    # rigidBodyResu = CREA_RESU( OPERATION='AFFE',
    #                     TYPE_RESU='MODE_MECA',
    #                     AFFE=List_F,
    #                     )

    # --
    # -- Calcul des masses et inerties autour du point CG donné
    # --

    MasseInertie = np.zeros(6)

    for n in range(nbRigidModes):
        # fi = rigidBodyResu.getField('DEPL',n+1).getValues()
        fi = rigidBodyChams[n].getValues()
        ve = M_Ud[n].getValues()
        print("vve", ve[:20])
        MasseInertie[n] = np.dot(fi, ve)

    print("MasseInertie", MasseInertie)

    # print('#---------------------#')
    # print('#--  Groupes de mailles :', group_ma)
    # print('#--  Masse :', MasseInertie[0:3])
    # print('#--  Inerties :', MasseInertie[3:])

    # -- Re normalisation des modes, des fois qu'ils soient pas normés par rapport à la masse

    PROJ_BASE(
        BASE=MODES,
        STOCKAGE="DIAG",  # -- pas la peine de calculer des caleurs croiées, on ne veux que la diagonale
        MATR_ASSE_GENE=(_F(MATRICE=CO("MRed"), MATR_ASSE=MASSE),),
    )
    NormMode = np.sqrt(np.diag(MRed.toNumpy()))

    # --
    # --    Calcul des produits Phi^T * (M * Ud)
    # --

    # --
    # -- Projection et calcul
    # --

    nbModes = MODES.getNumberOfIndexes()

    valeEffeUnit = np.zeros((nbModes, 6))
    valeEffe = np.zeros((nbModes, 6))
    for n in range(nbRigidModes):
        MUd = np.asarray(M_Ud[n].getValues())
        MUd[indZero] = 0.0
        Vale = np.zeros(nbModes)
        for m in range(nbModes):
            phi = MODES.getField("DEPL", m + 1).getValues()
            Vale[m] = np.dot(MUd, phi) / NormMode[m]

        valeEffe[:, n] = [Val**2 for Val in Vale]
        valeEffeUnit[:, n] = valeEffe[:, n] / MasseInertie[n]

    table = CREA_TABLE(
        LISTE=(
            _F(LISTE_I=MODES.getAccessParameters()["NUME_ORDRE"], PARA="NUME_ORDRE"),
            _F(LISTE_I=MODES.getAccessParameters()["NUME_MODE"], PARA="NUME_MODE"),
            _F(LISTE_R=MODES.getAccessParameters()["FREQ"], PARA="FREQ"),
            _F(LISTE_R=valeEffe[:, 0], PARA="MASS_EFFE_DX"),
            _F(LISTE_R=valeEffe[:, 1], PARA="MASS_EFFE_DY"),
            _F(LISTE_R=valeEffe[:, 2], PARA="MASS_EFFE_DZ"),
            _F(LISTE_R=valeEffe[:, 3], PARA="INER_EFFE_DX"),
            _F(LISTE_R=valeEffe[:, 4], PARA="INER_EFFE_DY"),
            _F(LISTE_R=valeEffe[:, 5], PARA="INER_EFFE_DZ"),
            _F(LISTE_R=valeEffeUnit[:, 0], PARA="MASS_EFFE_UN_DX"),
            _F(LISTE_R=valeEffeUnit[:, 1], PARA="MASS_EFFE_UN_DY"),
            _F(LISTE_R=valeEffeUnit[:, 2], PARA="MASS_EFFE_UN_DZ"),
            _F(LISTE_R=valeEffeUnit[:, 3], PARA="INER_EFFE_UN_DX"),
            _F(LISTE_R=valeEffeUnit[:, 4], PARA="INER_EFFE_UN_DY"),
            _F(LISTE_R=valeEffeUnit[:, 5], PARA="INER_EFFE_UN_DZ"),
        )
    )

    return table


def post_mode_ops(self, MODE, GROUP_MA, MATR_MASS, CARA_ELEM=None):

    """
    Macro-command POST_MODE
    """

    # args = _F(args)
    nume_ddl = MODE.getDOFNumbering()
    model = nume_ddl.getModel()
    massMatrix = MODE.getMassMatrix()

    mater = None
    cara = CARA_ELEM
    if massMatrix.getNumberOfElementaryMatrix() != 0:
        mater = massMatrix.getMaterialField()

    # mater = MODE.getMaterialField() : ne fonctionne pas si passage par calc_modes_multibandes
    if mater is None:
        UTMESS("F", "MODAL_25")

    group_ma = GROUP_MA
    t_list = isinstance(group_ma, list)
    t_tuple = isinstance(group_ma, tuple)
    if (not t_list) & (not t_tuple):
        group_ma = (group_ma,)

    gravityCenter = getGravityCenter(model, cara, mater, group_ma)
    rigidBodyChams, indZero = getRigidBodyModes(
        nume_ddl, group_ma, CG=gravityCenter, cara=cara, mater=mater
    )
    table = getModalMassIner(MODE, rigidBodyChams, indZero, MATR_MASS)

    return table
