# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
import numpy as np
from ..Messages import UTMESS
from ..Cata.Syntax import _F
from ..Supervis import CO
from ..Objects import Material, ExternalVariableTraits
from ..Commands import (
    DEFI_GROUP,
    AFFE_MATERIAU,
    AFFE_MODELE,
    STAT_NON_LINE,
    CALC_CHAMP,
    CREA_CHAMP,
    CREA_RESU,
    CREA_MAILLAGE,
    MACR_ADAP_MAIL,
    POST_RELEVE_T,
    POST_ELEM,
    CREA_TABLE,
    FORMULE,
    DETRUIRE,
)


def calc_stab_pente_ops(self, **args):
    if args["METHODE_STAB"] == "SRM":
        TABFS = calc_srm(self, args)
    else:
        if args["METHODE_LEM"] in ["BISHOP", "FELLENIUS"]:
            Solver = Surf_Circ_Solver(self, args)
        else:
            Solver = Surf_Non_Circ_Solver(self, args)

        TABFS = Solver.run()

    return TABFS


def calc_srm(self, args):
    """SRM solver

    Args:
        args (dict): Input parameters.

    Returns:
        Table Aster: Table containing the output results.
    """

    chmat = args["CHAM_MATER"]
    modele = args["MODELE"]

    # Vérification du maillage en entrée
    if not chmat.getMesh() == modele.getMesh():
        UTMESS("F", "CALCSTABPENTE_1")
    assert chmat.getMesh() == modele.getMesh()

    # Récupération du concept de maillage, de l'affectation des matériaux et des variables de commande
    mesh = chmat.getMesh()
    Para_Mat = get_para_mat(chmat)
    List_Affe_ExtVari = get_syntax_affevarc(chmat)

    # Vérification et analyse de l'affectation des matériaux dans la zone SRM
    dict_zone_deg = None
    has_zone_srm = False
    lgrma_cree = []

    if args.get("GROUP_MA") is not None:
        has_zone_srm = True
        lgrma_deg = args["GROUP_MA"]
        mesh, dict_zone_deg, lgrma_cree = check_srm_zone(lgrma_deg, Para_Mat)
        Para_Mat["MAILLAGE"] = mesh

    # Mots clé INVARIANT de l'AFFE_MATERIAU
    KW_Chmat_Fix = {}
    KW_Chmat_Fix["MAILLAGE"] = mesh
    if len(List_Affe_ExtVari) > 0:
        KW_Chmat_Fix["AFFE_VARC"] = List_Affe_ExtVari

    # Mots clé INVARIANT du STAT_NON_LINE
    KW_Solveur_Fix = {}
    KW_Solveur_Fix["MODELE"] = modele

    mc_increment = args.get("INCREMENT")
    KW_Solveur_Fix["INCREMENT"] = mc_increment
    if mc_increment["INST_FIN"] is not None:
        tfin = mc_increment["INST_FIN"]
    else:
        linst = mc_increment["LIST_INST"].getValues()
        if mc_increment["NUME_INST_FIN"] is not None:
            tfin = linst[mc_increment["NUME_INST_FIN"]]
        else:
            tfin = linst[-1]

    KW_Solveur_Fix["CONVERGENCE"] = args.get("CONVERGENCE")

    mc_cpt = args.get("COMPORTEMENT")

    CPT_UTIL = check_srm_cpt(mc_cpt, mesh, args)

    KW_Solveur_Fix["COMPORTEMENT"] = mc_cpt
    KW_Solveur_Fix["EXCIT"] = args.get("EXCIT")

    # Paramètres qui contrôlent la convergence de l'algorithme SRM
    para_fs = args["FS"][0]
    fact = para_fs["FS_INIT"]
    prec = para_fs["INCR_INIT"]
    prec_fin = para_fs["RESI_MAXI"]
    num_iter_max = para_fs["ITER_MAXI"]
    loi_vari = para_fs["METHODE"]

    if fact < 0:
        UTMESS("F", "CALCSTABPENTE_13")

    if prec_fin > prec:
        UTMESS("A", "CALCSTABPENTE_2", valr=prec)
        prec_fin = prec

    if loi_vari == "LINEAIRE":
        d_prec = (prec - prec_fin) / para_fs["ITER_RAFF_LINE"]
    else:
        d_prec = None

    num_iter = 0
    NC_fact = -1  # Indicateur de fact "divergant" qu'on vient de trouver

    __DISPT = FORMULE(NOM_PARA=("DX", "DY"), VALE="sqrt(DX*DX + DY*DY)")
    lfos = []
    lmaxdisp = []

    while True:
        if num_iter == num_iter_max:
            # Non-convergence
            UTMESS("F", "CALCSTABPENTE_3")
            break

        KW_Mat_Deg = gene_deg_kwmat(CPT_UTIL, Para_Mat, dict_zone_deg, fact, args)

        __CM_DEG = AFFE_MATERIAU(AFFE=KW_Mat_Deg, **KW_Chmat_Fix)

        try:
            __RESU = STAT_NON_LINE(CHAM_MATER=__CM_DEG, **KW_Solveur_Fix)
        except:
            if num_iter == 0:
                # FS_INIT trop grand.
                UTMESS("F", "CALCSTABPENTE_4")

            DETRUIRE(CONCEPT=_F(NOM=(__CM_DEG)))

            if prec - prec_fin > 1e-10:
                # Raffiner la précision suivant la loi indiquée
                NC_fact = fact
                fact, prec = refinement(fact, prec, d_prec)
                continue
            else:
                break

        # En cas de convergence --> Post-traitement

        __RESU = CALC_CHAMP(reuse=__RESU, RESULTAT=__RESU, VARI_INTERNE=("VARI_NOEU",))

        __CHDEP = CREA_CHAMP(
            INST=tfin, NOM_CHAM="DEPL", OPERATION="EXTR", RESULTAT=__RESU, TYPE_CHAM="NOEU_DEPL_R"
        )
        __CHFONC = CREA_CHAMP(
            AFFE=_F(NOM_CMP=("X1",), TOUT="OUI", VALE_F=(__DISPT,)),
            MAILLAGE=mesh,
            OPERATION="AFFE",
            TYPE_CHAM="NOEU_NEUT_F",
        )
        __CHDT = CREA_CHAMP(
            CHAM_F=__CHFONC, CHAM_PARA=(__CHDEP,), OPERATION="EVAL", TYPE_CHAM="NOEU_NEUT_R"
        )
        if has_zone_srm:
            # Calcul du depl_tot_maxi dans la zone srm
            __TABDT = POST_RELEVE_T(
                ACTION=_F(
                    CHAM_GD=__CHDT,
                    INTITULE="XX",
                    GROUP_MA=lgrma_deg,
                    OPERATION=("EXTREMA",),
                    TOUT="OUI",
                )
            )
        else:
            # Calcul du depl_tot_maxi dans le modèle entier
            __TABDT = POST_RELEVE_T(
                ACTION=_F(CHAM_GD=__CHDT, INTITULE="XX", OPERATION=("EXTREMA",), TOUT="OUI")
            )
        tabdt = __TABDT.EXTR_TABLE().values()
        lmaxdisp.append(tabdt["VALE"][2])
        lfos.append(fact)

        # MAJ des iterateurs
        fact += prec
        num_iter += 1
        if np.abs(NC_fact - fact) < prec:
            fact, prec = refinement(fact, prec, d_prec)
        DETRUIRE(CONCEPT=_F(NOM=(__CHDEP, __CHFONC, __CHDT, __TABDT)))

    # FIN de l'algorithme SRM
    if args.get("CHAM_DEFO") is not None:
        self.register_result(__RESU, args["CHAM_DEFO"])
    l_num_ord = [i + 1 for i in range(len(lfos))]
    TABFS = CREA_TABLE(
        LISTE=(
            _F(LISTE_I=l_num_ord, PARA="NUMERO"),
            _F(LISTE_R=lfos, PARA="FS"),
            _F(LISTE_R=lmaxdisp, PARA="DISP_TOT_MAXI"),
        )
    )

    # Suppression des modifications sur le maillage
    if len(lgrma_cree) > 0:
        mesh = DEFI_GROUP(DETR_GROUP_MA=_F(NOM=lgrma_cree), MAILLAGE=mesh)

    return TABFS


def refinement(fact, prec, d_prec):
    """Computes the refinement factor and precision.

    Args:
        fact (float): Current reduction factor.
        prec (float): Current reduction factor increment (precision).
        d_prec (float/None Type): Variation of precision if METHODE = 'LINEAIRE', None if not.

    Returns:
        tuple: Tuple containing the reduction factor and precision after refinement.
    """

    fact_raf = fact - prec
    if d_prec is None:
        prec_raf = prec / 2
    else:
        prec_raf = prec - d_prec
    fact_raf += prec_raf

    return fact_raf, prec_raf


def check_srm_cpt(mc_cpt, mesh, args):
    """Verify that there exist only one material behaviour
    among MOHR_COULOMB and DRUCKER_PRAGER in the SRM zone.

    Args:
        mc_cpt (list): List of factor keywords of the keyword 'COMPORTEMENT'.
        mesh (Mesh object): Mesh involved in the input cham_mater.
        args (dict): Input arguments of the macro-command.

    Returns:
        str: Name of the material behaviour used in the SRM algorithm.
    """

    cpt_auto = ["MOHR_COULOMB", "DRUCK_PRAGER"]
    cpt_exis = []
    if args.get("TOUT") is not None:
        for mcf_cpt in mc_cpt:
            if mcf_cpt["RELATION"] not in cpt_auto:
                UTMESS("F", "CALCSTABPENTE_5", valk=mcf_cpt["RELATION"])
            cpt_exis.append(mcf_cpt["RELATION"])
    else:
        lgrma_deg = args["GROUP_MA"]
        for mcf_cpt in mc_cpt:
            if mcf_cpt.get("MAILLE"):
                UTMESS("F", "CALCSTABPENTE_6")

            isGPMA = mcf_cpt.get("GROUP_MA")

            cpt = mcf_cpt["RELATION"]
            if isGPMA is None:
                if cpt not in cpt_auto:
                    UTMESS("F", "CALCSTABPENTE_5", valk=cpt)
                cpt_exis.append(cpt)
            else:
                lgrma_cpt = mcf_cpt["GROUP_MA"]
                for grma in lgrma_deg:
                    lgrma_inter = check_ma_intersec(mesh, grma, lgrma_cpt)
                    if len(lgrma_inter) > 1:
                        if cpt not in cpt_auto:
                            UTMESS("F", "CALCSTABPENTE_5", valk=cpt)
                        cpt_exis.append(cpt)

    cpt_exis = list(set(cpt_exis))

    if len(cpt_exis) != 1:
        UTMESS("F", "CALCSTABPENTE_7")

    return cpt_exis[0]


def check_srm_zone(lgrma_deg, Para_Mat):
    """In case of SRM analysis applied to a part of the model,
    verify the existence of the group_mas defining the SRM zone in the group_mas involved
    in the input cham_mater.
    if exist --> copy the group_ma,
    if not exist --> calculate the sous-groupe_ma if intersection detected.

    Args:
        lgrma_deg (list[str]): List of the group_ma defining the SRM zone.
        Para_Mat (dict): Materials and the related groupe_mas involved in the input cham_mater.

    Returns:
        tuple: Tuple containing the mesh enriched by the sous-group_mas,
            dict with materials in the SRM zone as keys and the concerned group_mas as values,
            and the list of names of created group_mas.
    """

    mesh = Para_Mat["MAILLAGE"]
    mat_zone = Para_Mat["AFFE_MATER"]
    dict_zone_deg = {}
    motcles_creagrma = []
    lgrma_cree = []

    lmat = list(mat_zone.keys())

    for imat in range(len(lmat)):
        mat = lmat[imat]
        if mat_zone[mat] == "TOUT":
            dict_zone_deg[mat] = lgrma_deg
        else:
            lgrma_srm = []
            for grma_deg in lgrma_deg:
                if grma_deg in mat_zone[mat]:
                    lgrma_srm.append(grma_deg)
                else:
                    lgrma_inter = check_ma_intersec(mesh, grma_deg, mat_zone[mat])

                    if len(lgrma_inter) > 1:
                        # mat se trouve dans la zone SRM --> Créer grma d'intersec
                        nom_grma = grma_deg + "_mat" + str(imat) + "_SRM"
                        if len(nom_grma) > 24:
                            UTMESS("F", "CALCSTABPENTE_8", valk=grma_deg)
                        motcles_creagrma.append(_F(INTERSEC=lgrma_inter, NOM=nom_grma))
                        lgrma_srm.append(nom_grma)

            if len(lgrma_srm) > 0:
                dict_zone_deg[mat] = lgrma_srm
                lgrma_cree.extend(lgrma_srm)

    if motcles_creagrma:
        mesh = DEFI_GROUP(CREA_GROUP_MA=motcles_creagrma, MAILLAGE=mesh)

    return mesh, dict_zone_deg, lgrma_cree


def check_ma_intersec(mesh, grma_obj, lgrma):
    """Detect the intersection between a group_ma and a list of group_mas.

    Args:
        mesh (Mesh object): Mesh involved in the input cham_mater.
        grma_obj (str): Name of the group_ma to be checked.
        lgrma (list[str]): List of the groups_mas that grma_obj may intersect with.

    Returns:
        list: List of the group_mas in lgrma presenting common cells with grma_obj.
    """
    grma_all = mesh.getGroupsOfCells()
    lgrma_inter = [grma_obj]

    if grma_obj not in grma_all:
        UTMESS("F", "CALCSTABPENTE_9", valk=grma_obj)

    nugrma_obj = mesh.getCells(grma_obj)

    for grma in lgrma:
        nugrma = mesh.getCells(grma)
        nucom = list(set(nugrma) & set(nugrma_obj))

        if len(nucom) > 0:
            lgrma_inter.append(grma)

    return lgrma_inter


def gene_deg_kwmat(CPT_SRM, Para_Mat, dict_zone_deg, fact, args):
    """Generate the factor keywords 'AFFE' of the operator AFFE_MATERIAU
    which produces the degraded cham_mater.

    Args:
        CPT_SRM (str): Material behaviour used in the SRM algorithm.
        Para_Mat (dict): Input material field parameters.
        dict_zone_deg (dict): Degraded material field parameters in the SRM zone.
        fact (float): Strength reduction factor.
        args (dict): Arguments of CALC_STAB_PENTE.

    Returns:
        dict: Dict of the factor keywords of AFFE_MATERIAU.
    """

    mat_zone = Para_Mat["AFFE_MATER"]
    KW_Mat_Deg = []

    if args.get("TOUT") is not None:
        # La zone SRM est la totalité du maillage
        for mat in list(mat_zone.keys()):
            nom_mats = mat.getMaterialNames()

            # Vérifier l'existance de CPT_SRM dans les cpts du matériau
            if CPT_SRM not in nom_mats:
                UTMESS("F", "CALCSTABPENTE_10")

            # Générer le matériau dégradé
            Matdeg = gene_deg_mat(mat, CPT_SRM, fact)
            if mat_zone[mat] == "TOUT":
                KW_Mat_Deg.append(_F(MATER=(Matdeg,), TOUT="OUI"))
            else:
                KW_Mat_Deg.append(_F(MATER=(Matdeg,), GROUP_MA=mat_zone[mat]))

    elif args.get("GROUP_MA") is not None:
        # Affectation par le principe de surcharge

        # Etape 1 : Reproduction du cham_mater saint
        for mat in list(mat_zone.keys()):
            if mat_zone[mat] == "TOUT":
                KW_Mat_Deg.append(_F(MATER=(mat,), TOUT="OUI"))
            else:
                KW_Mat_Deg.append(_F(MATER=(mat,), GROUP_MA=mat_zone[mat]))

        # Etape 2 : Dégradation dans la zone SRM
        for mat in list(dict_zone_deg.keys()):
            if CPT_SRM not in mat.getMaterialNames():
                UTMESS("F", "CALCSTABPENTE_10")

            Matdeg = gene_deg_mat(mat, CPT_SRM, fact)
            KW_Mat_Deg.append(_F(MATER=(Matdeg,), GROUP_MA=dict_zone_deg[mat]))

    else:
        raise TypeError("At least {0} or {1} is required".format("TOUT", "GROUP_MA"))

    return KW_Mat_Deg


def gene_deg_mat(OriginMat, cpt, fact):
    """Generate the Material object with the degraded material behaviour parameters.

    Args:
        OriginMat (Material object): Undegraded materiau.
        cpt (str): Material behaviour used in the SRM algorithm.
        fact (float): Strength reduction factor.

    Returns:
        Material object: Degraded material object.
    """

    Matdeg = Material(OriginMat, [cpt])

    propval = {}
    if cpt == "MOHR_COULOMB":
        Phi_0 = OriginMat.getValueReal(cpt, "PHI")
        Psi_0 = OriginMat.getValueReal(cpt, "ANGDIL")
        if Phi_0 != Psi_0:
            UTMESS("A", "CALCSTABPENTE_11")
        Phi_deg = np.arctan(np.tan(Phi_0 / 180 * np.pi) / fact) * 180 / np.pi
        propval["ANGDIL"] = Phi_deg
        propval["COHESION"] = OriginMat.getValueReal(cpt, "COHESION") / fact
        propval["PHI"] = Phi_deg

    elif cpt == "DRUCK_PRAGER":
        # Réduction c-phi
        A = OriginMat.getValueReal(cpt, "ALPHA")
        SY = OriginMat.getValueReal(cpt, "SY")
        PSI = OriginMat.getValueReal(cpt, "DILAT")
        propval["P_ULTM"] = OriginMat.getValueReal(cpt, "P_ULTM")
        num_ecroui = OriginMat.getValueReal(cpt, "TYPE_DP")
        if num_ecroui == 1.0:
            ECROUI = "LINEAIRE"
            propval["H"] = OriginMat.getValueReal(cpt, "H")
        elif num_ecroui == 2.0:
            ECROUI = "PARABOLIQUE"
            propval["SY_UTLM"] = OriginMat.getValueReal(cpt, "SY_UTLM")
        propval["ECROUISSAGE"] = ECROUI

        tan_phi_deg = 3 * A / (2 * np.sqrt((2 * A + 1) * (1 - A))) / fact
        c_deg = SY / (2 * np.sqrt((2 * A + 1) * (1 - A))) / fact
        psi_deg = np.arctan(np.tan(PSI / 180 * np.pi) / fact) * 180 / np.pi

        sin = tan_phi_deg / np.sqrt(1 + tan_phi_deg**2)
        cos = 1.0 / np.sqrt(1 + tan_phi_deg**2)

        propval["ALPHA"] = 2 * sin / (3 - sin)
        propval["SY"] = 6 * c_deg * cos / (3 - sin)
        propval["DILAT"] = psi_deg

    Matdeg.addProperties(cpt, **propval)

    return Matdeg


def get_para_mat(chmat):
    """Extract the mesh, materials and the associated group_mas from the input cham_mater.

    Args:
        chmat (MaterialField object): Input material field.

    Returns:
        dict: Dict containing the mesh and the material affectation parameters.
    """

    para_mat = {}
    para_mat["MAILLAGE"] = chmat.getMesh()
    para_mat["AFFE_MATER"] = {}

    # Récupérer les matériaux et les entités de maillage associées
    MatOnMeshEnt = chmat.getMaterialsOnMeshEntities()
    for curIter in MatOnMeshEnt:
        list_mat = curIter[0]
        if len(list_mat) > 1:
            UTMESS("F", "CALCSTABPENTE_12")

        meshEnt = curIter[1]
        if str(meshEnt.getType()) == "EntityType.GroupOfCellsType":
            para_mat["AFFE_MATER"][list_mat[0]] = meshEnt.getNames()
        else:
            para_mat["AFFE_MATER"][list_mat[0]] = "TOUT"

    return para_mat


def get_syntax_affevarc(chmat):
    """Extract the command variables from the input material field.

    Args:
        chmat (MaterialField object): Input material field.

    Returns:
        dict: Dict containing the factor keywords of 'AFFE_VARC'.
    """

    ExtVariAffe = chmat.getExtStateVariablesOnMeshEntities()
    List_Affe_ExtVari = []
    if len(ExtVariAffe) > 0:
        for curIter in ExtVariAffe:
            dict_extvari = {}
            # Reconstruction des syntaxes liées à l'objet ExternalStateVariable
            ExtVari = curIter[0]
            Nom_varc = ExternalVariableTraits.getExternVarTypeStr(ExtVari.getType())
            dict_extvari["NOM_VARC"] = Nom_varc

            inputField = ExtVari.getField()
            evolParam = ExtVari.getEvolutionParameter()

            if inputField:
                dict_extvari["CHAM_GD"] = inputField
            if evolParam:
                dict_extvari["EVOL"] = evolParam.getTransientResult()
                dict_extvari["PROL_GAUCHE"] = evolParam.getLeftExtension()
                dict_extvari["PROL_DROITE"] = evolParam.getRightExtension()
                if evolParam.getFieldName():
                    dict_extvari["NOM_CHAM"] = evolParam.getFieldName()
                if evolParam.getTimeFormula():
                    dict_extvari["FONC_INST"] = evolParam.getTimeFormula()
                if evolParam.getTimeFunction():
                    dict_extvari["FONC_INST"] = evolParam.getTimeFunction()
            if ExtVari.isSetRefe():
                dict_extvari["VALE_REF"] = ExtVari.getReferenceValue()

            # Reconstruction des syntaxes liées aux entités du maillage
            meshEnt = curIter[1]
            if str(meshEnt.getType()) == "EntityType.GroupOfCellsType":
                dict_extvari["GROUP_MA"] = meshEnt.getNames()
            else:
                dict_extvari["TOUT"] = "OUI"

            List_Affe_ExtVari.append(_F(dict_extvari))

    return List_Affe_ExtVari


def get_coord(mesh):
    """Get the nodal coordinates of the input mesh.

    Args:
        mesh (Mesh object): Input mesh.

    Returns:
        list[float]: Array of dimension (nb_nodes, 3) containing the nodal coordinates.
    """
    coord = mesh.getCoordinates().getValues()
    coord = np.array(coord).reshape((-1, 3))
    return coord


def zone_gliss(X, Y, resu_stab):
    """Indicate if the given location belongs to the sliding mass or not.

    Args:
        X (float): Abscissa of the interested location.
        Y (float): Ordinate of the interested location.
        resu_stab (dict): Result dictionnary of the LEM calculation.

    Returns:
        float: one if the location belongs to the sliding mass and zero if not.
    """

    if resu_stab.get("RAYON") is not None:
        # Surface circulaire
        Rayon = resu_stab["RAYON"]
        x_center = resu_stab["CENTRE_X"]
        y_center = resu_stab["CENTRE_Y"]
        if Rayon**2 - (X - x_center) ** 2 - (Y - y_center) ** 2 >= 1e-6:
            return 1.0
        else:
            return 0.0

    else:
        # Surface non-circulaire
        coor_X = resu_stab["COOR_X"]
        coor_Y = resu_stab["COOR_Y"]
        for n in range(np.size(coor_X, 0) - 1):
            if X >= coor_X[n] and X < coor_X[n + 1]:
                k = (coor_Y[n + 1] - coor_Y[n]) / (coor_X[n + 1] - coor_X[n])
                b = coor_Y[n] - k * coor_X[n]
                if (Y - k * X - b) >= 1e-6:
                    return 1.0
                else:
                    return 0.0
        return 0.0


def outline_slice_circ(X, Y, x0, y0, x1, x2, R, N, sigma):
    """Calculate the indicator of the proximity to the slice outline for circular failure surface.

    Args:
        X (float): Abscissa of the input point.
        Y (float): Ordiante of the input point.
        x0 (float): Abscissa of the sliding circle centre.
        y0 (float): Ordinate of the sliding circle centre.
        x1 (float): Abscissa of the sliding circle endpoint 1.
        x2 (float): Abscissa of the sliding circle endpoint 2.
        R (float): Radius of the sliding circle.
        N (int): Number of slices.
        sigma (float): Standard deviation of gaussian distribution.

    Returns:
        float: Value of proximity indicator.
    """

    d_min = np.inf
    width = np.abs(x2 - x1) / N

    d_cerc = R - np.sqrt((X - x0) ** 2 + (Y - y0) ** 2)
    if d_cerc > 1e-8:
        # maille dans la partie glissante
        for n in range(N - 1):
            d_tran = np.abs(X - np.min([x1, x2]) - width * (n + 1))
            if d_tran < d_min:
                d_min = d_tran
        if d_cerc < d_min:
            d_min = d_cerc
    else:
        d_min = np.abs(d_cerc)

    return 1.0 / np.sqrt(2 * np.pi) / sigma * np.exp(-(d_min**2) / 2 / sigma**2)


def outline_slice_nc(X, Y, para_geom, N, sigma):
    """Calculate the indicator of the proximity to the slice outline for non-circular failure surface.

    Args:
        X (float): Abscissa of the input point.
        Y (float): Ordiante of the input point.
        para_geom (list[float]): Array containing the coordinates of the intermediate points of the failure surface.
        N (int): Number of slices.
        sigma (float): Standard deviation of gaussian distribution.

    Returns:
        float: Value of proximity indicator.
    """

    d_min = np.inf
    # si x_enter < X < x_sorti, on calcule d_min, infini sinon
    for n in range(N):
        if X >= para_geom[n, 0] and X < para_geom[n + 1, 0]:
            d_tran_l = np.abs(X - para_geom[n, 0])
            d_tran_r = np.abs(X - para_geom[n + 1, 0])
            if n == 0:
                d_tran = d_tran_r
            elif n == N - 1:
                d_tran = d_tran_l
            else:
                d_tran = np.min([d_tran_l, d_tran_r])

            k = (para_geom[n + 1, 1] - para_geom[n, 1]) / (para_geom[n + 1, 0] - para_geom[n, 0])
            b = para_geom[n, 1] - k * (para_geom[n, 0])
            d_tran_base = (Y - k * X - b) / np.sqrt(k**2 + 1)
            if d_tran_base > 1e-8:
                d_min = np.min([d_tran, d_tran_base])
            else:
                d_min = np.abs(d_tran_base)

    return 1.0 / np.sqrt(2 * np.pi) / sigma * np.exp(-(d_min**2) / 2 / sigma**2)


class LEM_Solver:
    """Solver object containing the LEM-related methods."""

    def __init__(self, parent, args):
        self.parent = parent

        self.surf_circ = False
        self.bishop = False
        self.spencer = False
        self.chptot = None

        if args["METHODE_LEM"] in ["BISHOP", "FELLENIUS"]:
            self.surf_circ = True
            if args["METHODE_LEM"] == "BISHOP":
                self.bishop = True
        else:
            if args["METHODE_LEM"] == "SPENCER":
                self.spencer = True

        self.cham_defo = args["CHAM_DEFO"]

        self.chmat = args["CHAM_MATER"]
        para_mat = get_para_mat(self.chmat)
        self.affe_mat = para_mat["AFFE_MATER"]
        self.mesh = self.chmat.getMesh()
        self.connex = self.mesh.getMedConnectivity()

        # Récupérer le champ PTOT
        # La pression interstitielle est supposée positive
        para_ext_vari = get_syntax_affevarc(self.chmat)
        for kw_varc in para_ext_vari:
            if kw_varc["NOM_VARC"] == "PTOT":
                self.chptot = kw_varc["CHAM_GD"]
                chptot = kw_varc["CHAM_GD"].toSimpleFieldOnNodes()
                chptotval, mask = chptot.getValues()
                self.chptotval = np.array(chptotval).astype(np.float64)

                if kw_varc.get("TOUT") is None:
                    UTMESS("F", "CALCSTABPENTE_15")
                # On considère seulement le premier PTOT s'il en existe plusieurs
                break

        x1_min = args["X1_MINI"]
        x1_max = args["X1_MAXI"]
        x2_min = args["X2_MINI"]
        x2_max = args["X2_MAXI"]

        # Vérifier la légitimité des x1 et x2
        if x1_min > x1_max or x2_min > x2_max or x1_max > x2_min:
            UTMESS("F", "CALCSTABPENTE_14")

        self.x_lim_1 = [x1_min, x1_max]
        self.x_lim_2 = [x2_min, x2_max]

        self.nb_tran = args["NB_TRANCHE"]
        mc_raff_mail = args["RAFF_MAIL"]
        self.nb_max_adap = mc_raff_mail["NB_RAFF_MAXI"]
        self.resi_FS = mc_raff_mail["RAFF_CRIT_STAB"]
        self.cpt = args["CRITERE"]

        nom_grma_pente = args["GROUP_MA"]
        self.pente = self.get_coord_pente(nom_grma_pente)
        if x1_min < np.min(self.pente[:, 0]) or x2_max > np.max(self.pente[:, 0]):
            UTMESS("F", "CALCSTABPENTE_16")

        self.coord = get_coord(self.mesh)
        self.y_bas = np.min(self.coord[:, 1])
        self.centre = self.calc_cell_centre()

    def get_coord_pente(self, grma_pente):
        """Get the coordinates of the slope surface.

        Args:
            grma_pente (str): Name of the group of 1D elements representing the slope surface.

        Returns:
            list[float]: Array of dimension (nb_nodes, 3) containing the coordinates of the nodes on the slope surface.
        """

        # Coordonnées des neouds sur la pente
        mesh_pente = CREA_MAILLAGE(MAILLAGE=self.mesh, RESTREINT=_F(GROUP_MA=(grma_pente,)))

        coord_pente = mesh_pente.getCoordinates().getValues()
        coord_pente = np.array(coord_pente).reshape((-1, 3))
        indsort = np.argsort(coord_pente[:, 0])
        coord_pente = coord_pente[indsort, :]

        return coord_pente

    def calc_cell_centre(self):
        """Calculate the centre coordinates of mesh cells.

        Returns:
            list[float]: Array of dimension (nb_cells, 2) containing the planar coordinates of centre of the cells.
        """

        connex = self.mesh.getMedConnectivity()
        meshtypes = self.mesh.getMedCellsTypes()

        center = np.zeros((len(connex), 2))
        for ncell in range(len(connex)):
            if meshtypes[ncell] // 100 == 1:
                # éliminer les éléments SEG où le matériau n'est pas défini
                center[ncell, :] = [np.inf, np.inf]
                continue
            x = 0.0
            y = 0.0
            for nnode in connex[ncell]:
                x += self.coord[nnode - 1, 0]
                y += self.coord[nnode - 1, 1]

            center[ncell, :] = np.array([x, y]) / len(connex[ncell])

        return center

    def solve_circle(self, x_enter, x_sorti, R, yb=None):
        """Calculate the centre coordinates of a circle.

        Args:
            x_enter (float): Abscissa of the enter point.
            x_sorti (float): Abscissa of the exit point .
            R (float): Radius of the circle.
            yb (float, optional): Ordinate of the horizontal tangential line. Defaults to None.

        Returns:
            float, float: Coordinates of the circle centre.
        """

        x1 = x_enter
        x2 = x_sorti
        y1 = np.interp(x_enter, self.pente[:, 0], self.pente[:, 1])
        y2 = np.interp(x_sorti, self.pente[:, 0], self.pente[:, 1])

        # x0 = c1 - c2y0
        c1 = (x2**2 - x1**2 + y2**2 - y1**2) / 2 / (x2 - x1)
        c2 = (y2 - y1) / (x2 - x1)

        # Ay0^2 + By0 + C = 0
        if yb is None:
            A = c2**2 + 1
            B = 2 * (x1 - c1) * c2 - 2 * y1
            C = (x1 - c1) ** 2 + y1**2 - R**2
            delta = B**2 - 4 * A * C
            if np.abs(delta) < 1e-6:
                delta = 0.0
            else:
                if delta < 0:
                    raise Exception(
                        "Le rayon doit être supérieur à la distance entre les deux points données."
                    )
            y0 = (-B + np.sqrt(delta)) / 2 / A

        else:
            A = c2**2
            B = 2 * (x1 - c1) * c2 - 2 * y1 + 2 * yb
            C = (x1 - c1) ** 2 + y1**2 - yb**2
            delta = B**2 - 4 * A * C
            if np.abs(delta) < 1e-6:
                delta = 0.0
            else:
                if delta < 0:
                    raise Exception("Echec solve circle.")
            y0 = (-B - np.sqrt(delta)) / 2 / A

        x0 = c1 - c2 * y0

        return x0, y0

    def get_properties_on_surface(self, centre_base):
        """Get the material properties and pore pressure on the failure surface.

        Args:
            centre_base (list[float]): Array of dimension (nb_tran, 2) containing the coordiantes of the centre of slices' base.

        Returns:
            list[float], list[float]: Material properties array (nb_tran, 2) and pore pressure array (nb_tran).
        """

        mat_prop = np.zeros((self.nb_tran, 2))
        ptot = np.zeros(self.nb_tran)

        for ntran in range(np.size(centre_base, 0)):
            dist_min = np.inf
            for ncell in range(np.size(self.centre, 0)):
                dist = np.linalg.norm(centre_base[ntran, :] - self.centre[ncell, :])
                if dist < dist_min:
                    dist_min = dist
                    ncell_base = ncell

            # Material property at the base center
            mat_tran = self.chmat.getMaterialOnCell(self.mesh.getCellName(ncell_base))

            # Vérifier que le comportement existe dans le matériau
            nom_mats = mat_tran.getMaterialNames()
            if self.cpt not in nom_mats:
                UTMESS("F", "CALCSTABPENTE_10")

            # Récupération des propriétés de résistance
            if self.cpt == "MOHR_COULOMB":
                mat_prop[ntran, :] = [
                    mat_tran.getValueReal("MOHR_COULOMB", "PHI") * np.pi / 180.0,
                    mat_tran.getValueReal("MOHR_COULOMB", "COHESION"),
                ]
            else:
                A = mat_tran.getValueReal("DRUCK_PRAGER", "ALPHA")
                SY = mat_tran.getValueReal("DRUCK_PRAGER", "SY")
                mat_prop[ntran, 0] = np.arctan(3 * A / (2 * np.sqrt((2 * A + 1) * (1 - A))))
                mat_prop[ntran, 1] = SY / (2 * np.sqrt((2 * A + 1) * (1 - A)))

            # Pression interstitielle à la base de la tranche
            if self.chptot is not None:
                l_nodes = self.connex[ncell_base][:]
                # Evaluation of field value on centre_base by IDW interpolation
                inv_dist_tot = 0.0
                p = 2
                for nnode in l_nodes:
                    dist = np.linalg.norm(centre_base[ntran, :] - self.coord[nnode - 1, :2])
                    if dist < 1e-10:
                        ptot[ntran] = self.chptotval[nnode - 1, 0]
                        break
                    ptot[ntran] += (1.0 / dist) ** p * self.chptotval[nnode - 1, 0]
                    inv_dist_tot += (1.0 / dist) ** p

                if inv_dist_tot > 0:
                    ptot[ntran] /= inv_dist_tot

        return mat_prop, ptot

    def calc_poids_tranche(self, mesh, centre_base, width_tran, alpha=None):
        """Calculate the deadweight of the soil slices.

        Args:
            mesh (Mesh object): Input mesh.
            centre_base (list[float]): Array of dimension (nb_tran, 2) containing the coordiantes of the centre of slices' base.
            width_tran (float): Uniform width of the slices.
            alpha (list[float], optional): Inclination angles of the slices' base. Defaults to None.

        Returns:
            list[float]: Array containing the deadweight of the soil slices.
        """

        poids = np.zeros(self.nb_tran)
        # Regénérer le cham_mater avec MAILRAF
        kw_mat = []
        for mat in list(self.affe_mat.keys()):
            if self.affe_mat[mat] == "TOUT":
                kw_mat.append(_F(MATER=(mat,), TOUT="OUI"))
            else:
                kw_mat.append(_F(MATER=(mat,), GROUP_MA=self.affe_mat[mat]))

        __CHMRAF = AFFE_MATERIAU(AFFE=kw_mat, MAILLAGE=mesh)

        __MODELE = AFFE_MODELE(
            AFFE=_F(MODELISATION=("D_PLAN",), PHENOMENE="MECANIQUE", TOUT="OUI"), MAILLAGE=mesh
        )

        nom_tran = ["TRAN_" + str(ntran) for ntran in range(self.nb_tran)]
        nom_grma_del = []
        nom_grno_del = []

        for ntran in range(self.nb_tran):
            if self.surf_circ:
                # Surface circulaire - Méthode Fellenius et Bishop
                nom_grma_del += [nom_tran[ntran], nom_tran[ntran] + "_verti"]

                # Si on prend les tranches avec courbure
                if ntran == 0:
                    mesh = DEFI_GROUP(
                        reuse=mesh,
                        CREA_GROUP_MA=(
                            _F(
                                ANGL_NAUT=(0.0,),
                                DIST=width_tran / 2,
                                NOM=nom_tran[ntran] + "_verti",
                                OPTION="BANDE",
                                POINT=centre_base[ntran, :],
                            ),
                            _F(INTERSEC=(nom_tran[ntran] + "_verti", "GLISS"), NOM=nom_tran[ntran]),
                        ),
                        MAILLAGE=mesh,
                    )
                else:
                    mesh = DEFI_GROUP(
                        reuse=mesh,
                        CREA_GROUP_MA=(
                            _F(
                                ANGL_NAUT=(0.0,),
                                DIST=width_tran / 2,
                                NOM=nom_tran[ntran] + "_verti",
                                OPTION="BANDE",
                                POINT=centre_base[ntran, :],
                            ),
                            _F(
                                INTERSEC=(nom_tran[ntran] + "_verti", "GLISS"),
                                NOM=nom_tran[ntran] + "_inter",
                            ),
                            _F(
                                DIFFE=(nom_tran[ntran] + "_inter", nom_tran[ntran - 1]),
                                NOM=nom_tran[ntran],
                            ),
                        ),
                        MAILLAGE=mesh,
                    )

                    nom_grma_del.append(nom_tran[ntran] + "_inter")

            else:
                # Surface non-circulaire - Méthode Morgenstern
                nom_grma_del += [
                    nom_tran[ntran],
                    nom_tran[ntran] + "_droite",
                    nom_tran[ntran] + "_verti",
                ]
                nom_grno_del.append(nom_tran[ntran] + "_N_droite")

                # y = kx + yb pour la droite de la base des tranches
                k = np.tan(alpha[ntran])
                yb = centre_base[ntran, 1] - k * centre_base[ntran, 0]

                nodes_gliss = []
                coord_raf = get_coord(mesh)
                for nnode in range(np.size(coord_raf, 0)):
                    x_node = coord_raf[nnode, 0]
                    y_node = coord_raf[nnode, 1]
                    if y_node - k * x_node - yb > 1e-6:
                        nodes_gliss.append(nnode)

                mesh.setGroupOfNodes(nom_tran[ntran] + "_N_droite", nodes_gliss)

                kw_crea_group_ma = [
                    _F(
                        OPTION="APPUI",
                        GROUP_NO=nom_tran[ntran] + "_N_droite",
                        TYPE_APPUI="TOUT",
                        NOM=nom_tran[ntran] + "_droite",
                    ),
                    _F(
                        ANGL_NAUT=(0.0,),
                        DIST=width_tran / 2,
                        NOM=nom_tran[ntran] + "_verti",
                        OPTION="BANDE",
                        POINT=centre_base[ntran, :],
                    ),
                ]

                if ntran == 0:
                    kw_crea_group_ma.append(
                        _F(
                            INTERSEC=(nom_tran[ntran] + "_verti", nom_tran[ntran] + "_droite"),
                            NOM=nom_tran[ntran],
                        )
                    )
                    mesh = DEFI_GROUP(reuse=mesh, CREA_GROUP_MA=kw_crea_group_ma, MAILLAGE=mesh)
                else:
                    kw_crea_group_ma += [
                        _F(
                            INTERSEC=(nom_tran[ntran] + "_verti", nom_tran[ntran] + "_droite"),
                            NOM=nom_tran[ntran] + "_inter",
                        ),
                        _F(
                            DIFFE=(nom_tran[ntran] + "_inter", nom_tran[ntran - 1]),
                            NOM=nom_tran[ntran],
                        ),
                    ]
                    mesh = DEFI_GROUP(reuse=mesh, CREA_GROUP_MA=kw_crea_group_ma, MAILLAGE=mesh)

                    nom_grma_del.append(nom_tran[ntran] + "_inter")

            __TABM = POST_ELEM(
                CHAM_MATER=__CHMRAF, MODELE=__MODELE, MASS_INER=_F(GROUP_MA=(nom_tran[ntran],))
            )

            poids[ntran] = __TABM["MASSE", 1] * 9.81

            DETRUIRE(CONCEPT=_F(NOM=(__TABM)))

        if not self.surf_circ:
            mesh = DEFI_GROUP(reuse=mesh, MAILLAGE=mesh, DETR_GROUP_NO=_F(NOM=nom_grno_del))

        mesh = DEFI_GROUP(DETR_GROUP_MA=_F(NOM=nom_grma_del), MAILLAGE=mesh, reuse=mesh)

        DETRUIRE(CONCEPT=_F(NOM=(__MODELE, __CHMRAF)))

        return poids

    def calc_bishop(self, x_enter, x_sorti, Rayon, MAILRAF=None):
        """Calculate the factor of security for circular failure surfaces.

        Args:
            x_enter (float): Abscissa of the enter point.
            x_sorti (float): Abscissa of the exit point.
            Rayon (float): Radius of the failure surface.
            MAILRAF (Mesh object, optional): Refined mesh. Defaults to None.

        Returns:
            float: Factor of security.
        """

        if MAILRAF is None:
            MAILRAF = self.mesh

        x_centre, y_centre = self.solve_circle(x_enter, x_sorti, Rayon)

        # Vérification acceptabilité cinématique de la surface critique
        for ind in range(np.size(self.pente, 0)):
            x_p = self.pente[ind, 0]
            y_p = self.pente[ind, 1]
            if (x_p - x_centre) ** 2 + (y_p - y_centre) ** 2 - Rayon**2 >= 0.0:
                if x_p - x_enter > 1e-6 and x_sorti - x_p > 1e-6:
                    # surface critique au-dessus de la pente
                    return None
        if self.y_bas + Rayon - y_centre > 1e-3:
            # surface critique hors du modèle
            return None

        # Division du modèle en tranches
        centre_base = np.zeros((self.nb_tran, 2))
        alpha = np.zeros(self.nb_tran)
        width_tran = np.abs(x_enter - x_sorti) / self.nb_tran

        # ============================================================
        # Calcul des alpha et propriété materiaux (maillage originel)
        # ============================================================
        for ntran in range(self.nb_tran):
            centre_base[ntran, 0] = (ntran + 1 / 2) * width_tran + np.min([x_enter, x_sorti])
            centre_base[ntran, 1] = y_centre - np.sqrt(
                Rayon**2 - (centre_base[ntran, 0] - x_centre) ** 2
            )

            tan_a = (
                1.0
                / np.sqrt(Rayon**2 - (centre_base[ntran, 0] - x_centre) ** 2)
                * (centre_base[ntran, 0] - x_centre)
            )
            alpha[ntran] = np.arctan(tan_a)

        # Identifier les propriétés matériaux sur la surface critique
        mat_prop, ptot = self.get_properties_on_surface(centre_base)

        # =============================================
        # Calcul des poids propres (maillage raffiné)
        # =============================================

        # Création grma de la partie glissante et calcul des poids des tranches
        coord_raf = get_coord(MAILRAF)
        nodes_gliss = []
        for nnode in range(np.size(coord_raf, 0)):
            x_node = coord_raf[nnode, 0]
            y_node = coord_raf[nnode, 1]
            if Rayon**2 - (x_centre - x_node) ** 2 - (y_centre - y_node) ** 2 > 1e-6:
                nodes_gliss.append(nnode)

        MAILRAF.setGroupOfNodes("GLISS_N", nodes_gliss)

        MAILRAF = DEFI_GROUP(
            reuse=MAILRAF,
            CREA_GROUP_MA=_F(
                OPTION="APPUI", GROUP_NO="GLISS_N", TYPE_APPUI="TOUT", NOM="GLISS", TYPE_MAILLE="2D"
            ),
            MAILLAGE=MAILRAF,
        )

        poids = self.calc_poids_tranche(MAILRAF, centre_base, width_tran)

        MAILRAF = DEFI_GROUP(
            DETR_GROUP_MA=_F(NOM=("GLISS")),
            DETR_GROUP_NO=_F(NOM=("GLISS_N")),
            MAILLAGE=MAILRAF,
            reuse=MAILRAF,
        )
        # Calcul du FS
        if not self.bishop:
            FS_resu = np.sum(
                mat_prop[:, 1] * width_tran / np.cos(alpha)
                + (poids * np.cos(alpha) - ptot * width_tran * np.cos(alpha))
                * np.tan(mat_prop[:, 0])
            ) / np.abs(np.sum(poids * np.sin(alpha)))
        else:
            resi_max = 1e-6
            resi = 1.0
            FS = 1.0
            nb_iter_max = 1e3
            nb_iter = 1
            while resi > resi_max:
                FS_new = np.sum(
                    (
                        mat_prop[:, 1] * width_tran
                        + (poids - ptot * width_tran) * np.tan(mat_prop[:, 0])
                    )
                    / (np.cos(alpha) + (np.sin(alpha) * np.tan(mat_prop[:, 0]) / FS))
                ) / np.abs(np.sum(poids * np.sin(alpha)))

                resi = np.abs(FS_new - FS)
                FS = FS_new

                nb_iter += 1
                if nb_iter > nb_iter_max:
                    # Divergence iteration point fixe
                    return None

            FS_resu = FS

        return FS_resu

    def bilan_last_tranche(self, poids, width, alpha, mat_prop, ptot, spencer=True, ang_dir=None):
        """Solve the FS-lambda coupling equations with the condition of null force and moment on the right of the last slice.

        Args:
            poids (list[float]): Array containing the deadweight of the slices.
            width (float): Uniform width of the slices.
            alpha (list[float]): Array containing the inclination angles of the slices' base.
            mat_prop (list[float]): Array of material properties on the base of slices.
            ptot (list[float]): Array of pore pressure on the base of slices.
            spencer (bool, optional): Indicate if the Spencer procedure is employed. Defaults to True.
            ang_dir (list[float], optional): Tangent of the interslice force direction at the ends of the failure surface. Defaults to None.

        Returns:
            float: Factor of security.
        """

        N = np.size(poids, 0)
        max_resi = 1e-6
        max_resi_sub = 1e-6
        lamb = 0.0
        FS = 1.0
        resi_FS = 1.0
        resi_lamb = 1.0
        max_iter = 100
        max_iter_sub = 200
        nb_iter = 0

        # Résolution lambda et FS par itération point fixe
        while resi_FS > max_resi or resi_lamb > max_resi:
            nb_iter += 1

            lamb_old = lamb
            # Calcul des angles tan(beta) = f0(x) + lambda*sin(x) [Chen and Morgenstern, 1983]
            # Cardinal(f) = N
            if spencer:
                # f0 = 0 et f = 1 --> Méthode Spencer
                beta = np.arctan(np.ones(N) * lamb)
                f0 = np.zeros(N)
                f1 = np.ones(N)
            else:
                # f0 = variation linéaire, f = sin(pi/width*(x-x_min)) --> Méthode MP
                f0 = np.linspace(ang_dir[0], ang_dir[1], N + 1)[1:]
                f1 = np.sin(np.linspace(0.0, np.pi, N + 1)[1:])
                beta = np.arctan(f0 + f1 * lamb)

            FS_old = FS
            # Mise à jour du FS
            FS_sub = FS
            resi_sub = 1.0
            nb_iter_sub = 1.0
            while resi_sub > max_resi_sub:
                Phi = FS * np.cos(alpha - beta) + np.sin(alpha - beta) * np.tan(mat_prop[:, 0])
                F_Resist = (poids * np.cos(alpha) - ptot * width / np.cos(alpha)) * np.tan(
                    mat_prop[:, 0]
                ) + mat_prop[:, 1] * width / np.cos(alpha)
                F_Tangent = poids * np.sin(alpha)
                Psi = np.zeros(N - 1)

                numer = F_Resist[-1]
                denom = F_Tangent[-1]
                for n in range(N - 2, -1, -1):
                    Psi[n] = (
                        FS * np.cos(alpha[n + 1] - beta[n])
                        + np.sin(alpha[n + 1] - beta[n]) * np.tan(mat_prop[n + 1, 0])
                    ) / Phi[n]
                    PI_psi = 1.0
                    for j in range(n, N - 1):
                        PI_psi *= Psi[j]
                    numer += F_Resist[n] * PI_psi
                    denom += F_Tangent[n] * PI_psi

                FS = numer / np.abs(denom)
                resi_sub = np.abs(FS_sub - FS)
                FS_sub = FS
                if nb_iter_sub > max_iter_sub:
                    return 1e4
                else:
                    nb_iter_sub += 1

            resi_FS = np.abs(FS - FS_old)

            # Force d'interantion entre tranches, F_n = 0
            F_inter = np.zeros(N)
            F_inter[0] = (F_Resist[0] - FS * F_Tangent[0]) / Phi[0]
            for n in range(1, N - 1):
                F_inter[n] = (
                    Psi[n - 1] * F_inter[n - 1] * Phi[n - 1] + F_Resist[n] - FS * F_Tangent[n]
                ) / Phi[n]
            F_hor = F_inter * np.cos(beta)

            # Mise à jour de lambda
            numer_l = np.tan(alpha[0]) * F_hor[0] - f0[0] * F_hor[0]
            denom_l = f1[0] * F_hor[0]

            for n in range(1, N):
                numer_l += np.tan(alpha[n]) * (F_hor[n] + F_hor[n - 1]) - (
                    f0[n] * F_hor[n] + f0[n - 1] * F_hor[n - 1]
                )
                denom_l += f1[n] * F_hor[n] + f1[n - 1] * F_hor[n - 1]

            lamb = numer_l / denom_l
            resi_lamb = np.abs(lamb_old - lamb)

            if nb_iter == max_iter:
                UTMESS("I", "CALCSTABPENTE_17")
                return 1e4

        return FS

    def calc_morgenstern(self, vari_etat, MAILRAF=None):
        """Calculate the factor of security for non-circular failure surfaces.

        Args:
            vari_etat (list[float]): Array containing the abscissa of the end points and the ordinates of the intermediate points.
            MAILRAF (Mesh object, optional): Refined mesh. Defaults to None.

        Returns:
            float: Factor of security.
        """

        x_enter = vari_etat[0]
        x_sorti = vari_etat[-1]
        l_Y = vari_etat[1:-1]

        N = len(l_Y) + 1
        surf_X = np.linspace(x_enter, x_sorti, N + 1)
        surf_Y = np.array([np.interp(x_enter, self.pente[:, 0], self.pente[:, 1])])
        surf_Y = np.concatenate((surf_Y, np.array(l_Y)))
        surf_Y = np.append(surf_Y, np.interp(x_sorti, self.pente[:, 0], self.pente[:, 1]))

        # Division du modèle en tranches et calculer ALPHA, CENTRE_BASE, C et PHI
        centre_base = np.zeros((N, 2))
        alpha = np.zeros(N)
        width_tran = np.abs(x_enter - x_sorti) / N

        for ntran in range(N):
            centre_base[ntran, 0] = (surf_X[ntran] + surf_X[ntran + 1]) / 2
            centre_base[ntran, 1] = (surf_Y[ntran] + surf_Y[ntran + 1]) / 2
            tan_a = (surf_Y[ntran + 1] - surf_Y[ntran]) / (surf_X[ntran + 1] - surf_X[ntran])
            alpha[ntran] = np.arctan(tan_a)

        mat_prop, ptot = self.get_properties_on_surface(centre_base)

        # Calcul des poids propres des tranches et création du grma de la partie glissante
        if MAILRAF is None:
            MAILRAF = self.mesh

        poids = self.calc_poids_tranche(MAILRAF, centre_base, width_tran, alpha=alpha)

        # Résolution du FS
        if self.spencer:
            FS = self.bilan_last_tranche(poids, width_tran, alpha, mat_prop, ptot)
        else:
            tan_ang = []
            x_extrem = x_enter
            for nnode in range(np.size(self.pente, 0)):
                if len(tan_ang) > 0:
                    x_extrem = x_sorti
                if self.pente[nnode, 0] > x_extrem:
                    tan = (self.pente[nnode, 1] - self.pente[nnode - 1, 1]) / (
                        self.pente[nnode, 0] - self.pente[nnode - 1, 0]
                    )
                    tan_ang.append(tan)
                    if len(tan_ang) > 1:
                        break
            FS = self.bilan_last_tranche(
                poids, width_tran, alpha, mat_prop, ptot, spencer=False, ang_dir=tan_ang
            )

        return FS

    def crea_champ_pilo(self, para_geom, sigma, mesh):
        """Create the field of proximity to which the mesh will be adapted.

        Args:
            para_geom (list[float]): Geometric parameters of the failure surface.
            sigma (float): Standard deviation of gaussien distribution.
            mesh (Mesh object): Input mesh.

        Returns:
            cham_no object: Field of proximity.
        """

        if self.surf_circ:
            x_enter = para_geom["X_ENTER"]
            x_sorti = para_geom["X_SORTI"]
            R = para_geom["RAYON"]
            x_centre, y_centre = self.solve_circle(x_enter, x_sorti, R)
            __FOUTL = FORMULE(
                NOM_PARA=("X", "Y"),
                VALE="outline_slice_circ(X, Y, x_centre, y_centre, x_enter, x_sorti, R, nb_tran, sigma)",
                outline_slice_circ=outline_slice_circ,
                x_centre=x_centre,
                y_centre=y_centre,
                x_enter=x_enter,
                x_sorti=x_sorti,
                R=R,
                nb_tran=self.nb_tran,
                sigma=sigma,
            )
        else:
            __FOUTL = FORMULE(
                NOM_PARA=("X", "Y"),
                VALE="outline_slice_nc(X, Y, para_geom, nb_tran, sigma)",
                outline_slice_nc=outline_slice_nc,
                para_geom=para_geom,
                nb_tran=self.nb_tran,
                sigma=sigma,
            )

        __CHGEO = CREA_CHAMP(
            MAILLAGE=mesh, NOM_CHAM="GEOMETRIE", OPERATION="EXTR", TYPE_CHAM="NOEU_GEOM_R"
        )

        __CHFON = CREA_CHAMP(
            AFFE=_F(NOM_CMP=("X1"), TOUT="OUI", VALE_F=(__FOUTL)),
            MAILLAGE=mesh,
            OPERATION="AFFE",
            TYPE_CHAM="NOEU_NEUT_F",
        )

        __CHPILO = CREA_CHAMP(
            CHAM_F=__CHFON, CHAM_PARA=(__CHGEO,), OPERATION="EVAL", TYPE_CHAM="NOEU_NEUT_R"
        )

        return __CHPILO

    def raff_maillage(self, resu_stab):
        """Refine the mesh around the outline the soil slices.

        Args:
            resu_stab (dict): Dictionnary containing the result of failure surface searching.

        Returns:
            list[float]: Array containing the FS corresponding to the refined meshes.
        """

        l_FS = []
        FS = resu_stab["FS"]
        resi = np.inf
        k_adap = 0
        l_nom_mailraf = [None for i in range(self.nb_max_adap + 1)]
        l_mailraf_detr = []
        l_FS.append(FS)

        if self.surf_circ:
            x_enter_FS = resu_stab["X_ENTRE"]
            x_sorti_FS = resu_stab["X_SORTIE"]
            R_FS = resu_stab["RAYON"]
            para_geom = {"X_ENTER": x_enter_FS, "X_SORTI": x_sorti_FS, "RAYON": R_FS}
        else:
            coord_X = resu_stab["COOR_X"]
            coord_Y = resu_stab["COOR_Y"]
            x_enter_FS = coord_X[0]
            x_sorti_FS = coord_X[-1]
            para_geom = np.vstack((coord_X, coord_Y))
            para_geom = para_geom.transpose()
            # Retrouver l'état des variables optimal
            etat_opti = coord_Y.copy()
            etat_opti[0] = coord_X[0]
            etat_opti[-1] = coord_X[-1]

        sigma = 0.3 * np.abs(x_enter_FS - x_sorti_FS) / self.nb_tran

        MAILRAF = self.mesh

        while resi > self.resi_FS and k_adap < self.nb_max_adap:
            k_adap += 1
            l_nom_mailraf[k_adap] = "MAILRAF_" + str(k_adap)

            __CHPILO = self.crea_champ_pilo(para_geom, sigma, MAILRAF)

            # seuil_raff compte les mailles dont la somme de proximité = 68.2% (Très proches) --> Raffiner
            # seuil_dera compte les mailles dont la somme de proximite = 4.6% (Très loins) --> Déraffiner
            seuil_raff = 1.0 / np.sqrt(2 * np.pi) / sigma * np.exp(-0.5)
            seuil_dera = 1.0 / np.sqrt(2 * np.pi) / sigma * np.exp(-2.0)

            MACR_ADAP_MAIL(
                ADAPTATION="RAFF_DERA",
                MAILLAGE_N=MAILRAF,
                MAILLAGE_NP1=CO(l_nom_mailraf[k_adap]),
                CHAM_GD=__CHPILO,
                CRIT_RAFF_ABS=float(format(seuil_raff, ".4f")),
                CRIT_DERA_ABS=float(format(seuil_dera, ".4f")),
            )

            MAILRAF = globals()[l_nom_mailraf[k_adap]]
            l_mailraf_detr.append(MAILRAF)

            sigma /= 2

            if self.surf_circ:
                FS_raf = self.calc_bishop(x_enter_FS, x_sorti_FS, R_FS, MAILRAF=MAILRAF)
            else:
                FS_raf = self.calc_morgenstern(etat_opti, MAILRAF=MAILRAF)

            resi = np.abs(FS_raf - FS)
            FS = FS_raf
            l_FS.append(FS)

        if self.cham_defo is not None:
            self.regis_cham_defo(resu_stab, MAILRAF)

        DETRUIRE(CONCEPT=_F(NOM=l_mailraf_detr))

        return l_FS

    def regis_cham_defo(self, resu_stab, mesh):
        """Output the field visualizing the failure surface.

        Args:
            resu_stab (list[float]): Array containing the result of failure surface searching.
            mesh (Mesh object): Input mesh.
        """

        __FGLISS = FORMULE(
            NOM_PARA=("X", "Y"),
            VALE="zone_gliss(X, Y, resu_stab)",
            zone_gliss=zone_gliss,
            resu_stab=resu_stab,
        )

        __CHGEO = CREA_CHAMP(
            MAILLAGE=mesh, NOM_CHAM="GEOMETRIE", OPERATION="EXTR", TYPE_CHAM="NOEU_GEOM_R"
        )

        __CHFONC = CREA_CHAMP(
            AFFE=_F(NOM_CMP=("X1"), TOUT="OUI", VALE_F=(__FGLISS)),
            MAILLAGE=mesh,
            OPERATION="AFFE",
            TYPE_CHAM="NOEU_NEUT_F",
        )

        __CHEVAL = CREA_CHAMP(
            CHAM_F=__CHFONC, CHAM_PARA=(__CHGEO,), OPERATION="EVAL", TYPE_CHAM="NOEU_NEUT_R"
        )

        __CHDEPL = CREA_CHAMP(
            ASSE=_F(CHAM_GD=__CHEVAL, NOM_CMP=("X1"), NOM_CMP_RESU=("DX"), TOUT="OUI"),
            MAILLAGE=mesh,
            OPERATION="ASSE",
            TYPE_CHAM="NOEU_DEPL_R",
        )

        __REDEPL = CREA_RESU(
            AFFE=(_F(CHAM_GD=__CHDEPL, INST=(0.0,)),),
            NOM_CHAM="DEPL",
            OPERATION="AFFE",
            TYPE_RESU="EVOL_NOLI",
        )

        self.parent.register_result(__REDEPL, self.cham_defo)


class Surf_Circ_Solver(LEM_Solver):
    """Solver object containing the methods for circular failure surface searching and FS calculation."""

    def __init__(self, parent, args):
        super(Surf_Circ_Solver, self).__init__(parent, args)
        self.x_bande_1 = np.linspace(self.x_lim_1[0], self.x_lim_1[1], num=args["NB_POINT_1"])
        self.x_bande_2 = np.linspace(self.x_lim_2[0], self.x_lim_2[1], num=args["NB_POINT_2"])

        if args.get("Y_MINI") and args.get("Y_MAXI"):
            self.def_y_lim = True
        else:
            self.def_y_lim = False

        if self.def_y_lim:
            self.y_min = args["Y_MINI"]
            self.y_max = args["Y_MAXI"]
            # Vérification du y_min
            y1 = np.interp(self.x_lim_1[1], self.pente[:, 0], self.pente[:, 1])
            y2 = np.interp(self.x_lim_2[0], self.pente[:, 0], self.pente[:, 1])
            if self.y_min > np.min([y1, y2]):
                UTMESS("F", "CALCSTABPENTE_19")
        else:
            if args.get("Y_MINI") or args.get("Y_MAXI"):
                UTMESS("F", "CALCSTABPENTE_18")

    def search_surf_crit(self):
        """Search for the critical surface minimizing the FS.

        Returns:
            dict: Dictionary containing the result of failure surface searching.
        """
        resu_stab_circ = {}
        H_pente = np.max(self.pente[:, 1]) - np.min(self.pente[:, 1])
        incr_R = H_pente * 0.1
        FS_MIN = np.inf

        for x_enter in self.x_bande_1:
            for x_sorti in self.x_bande_2:
                y_enter = np.interp(x_enter, self.pente[:, 0], self.pente[:, 1])
                y_sorti = np.interp(x_sorti, self.pente[:, 0], self.pente[:, 1])
                y_lower = np.min([y_enter, y_sorti])
                if not self.def_y_lim:
                    # On calcule automatiquement les limites du rayon
                    R_test = np.linalg.norm(np.array([x_enter - x_sorti, y_enter - y_sorti])) / 2

                    if (y_enter + y_sorti) / 2 - R_test < self.y_bas:
                        # Commencer par le cercle tangent à la base
                        x_test, y_test = self.solve_circle(x_enter, x_sorti, " ", self.y_bas)
                        R_test = y_test - self.y_bas

                    x_max, y_max = self.solve_circle(x_enter, x_sorti, " ", y_lower)
                    R_max = y_max - y_lower
                    R_sup = 1.5 * R_max

                else:
                    # Les limites du rayon sont définies par l'utilisateur
                    if self.y_min > y_lower:
                        continue
                    x_test, y_test = self.solve_circle(x_enter, x_sorti, " ", self.y_min)
                    R_test = y_test - self.y_min
                    if self.y_max >= y_lower:
                        y_sup = y_lower
                    else:
                        y_sup = self.y_max
                    x_max, y_sup = self.solve_circle(x_enter, x_sorti, " ", y_sup)
                    R_sup = y_sup - y_lower

                # Eviter que l'incrément R dépasse la bande d'essai
                if R_sup - R_test < incr_R:
                    incr_R = (R_sup - R_test) / 2

                FS_0 = self.calc_bishop(x_enter, x_sorti, R_test)

                regis_fs = [FS_0]
                while R_test <= R_sup:
                    R_test += incr_R
                    FS = self.calc_bishop(x_enter, x_sorti, R_test)
                    if FS is None:
                        # surface critique non trouvee
                        break
                    regis_fs.append(FS)

                FS_loc_min = min(regis_fs)

                # Mise à jour du mimima global
                if FS_loc_min < FS_MIN:
                    ind_min = regis_fs.index(FS_loc_min)
                    x_enter_FS = x_enter
                    x_sorti_FS = x_sorti
                    R_FS = R_test - (len(regis_fs) - 1 - ind_min) * incr_R

                    # Trouver la vraie mimina par la méthode de bifurcation
                    resi = np.inf
                    n_div = 1
                    fs_stat = np.zeros(5)
                    for i in range(5):
                        fs_stat[i] = np.inf
                    if ind_min == 0:
                        fs_stat[1:3] = np.array(regis_fs[:2])
                    elif ind_min == len(regis_fs) - 1:
                        fs_stat[:2] = np.array(regis_fs[ind_min - 1 :])
                    else:
                        fs_stat[:3] = np.array(regis_fs[ind_min - 1 : ind_min + 2])

                    while resi > 1e-3:
                        if fs_stat[0] == np.inf:
                            fs_stat[4] = self.calc_bishop(
                                x_enter, x_sorti, R_FS + incr_R / 2**n_div
                            )
                        elif fs_stat[2] == np.inf:
                            fs_stat[3] = self.calc_bishop(
                                x_enter, x_sorti, R_FS - incr_R / 2**n_div
                            )
                        else:
                            fs_stat[3] = self.calc_bishop(
                                x_enter, x_sorti, R_FS - incr_R / 2**n_div
                            )
                            fs_stat[4] = self.calc_bishop(
                                x_enter, x_sorti, R_FS + incr_R / 2**n_div
                            )

                        # Eviter les surfaces illégales et enregistrer les surface testees
                        for ind in [3, 4]:
                            if fs_stat[ind] is None:
                                fs_stat[ind] = np.inf

                        FS = np.min(fs_stat)
                        if FS < FS_loc_min:
                            ind_min = np.argmin(fs_stat)
                            if ind_min == 3:
                                R_FS -= incr_R / 2**n_div
                                fs_stat[:3] = np.array([fs_stat[0], fs_stat[3], fs_stat[1]])
                            else:
                                R_FS += incr_R / 2**n_div
                                fs_stat[:3] = np.array([fs_stat[1], fs_stat[4], fs_stat[2]])

                        resi = np.abs(FS - FS_loc_min)
                        FS_loc_min = FS
                        n_div += 1

                    FS_MIN = FS_loc_min

        resu_stab_circ["RAYON"] = R_FS
        resu_stab_circ["X_ENTRE"] = x_enter_FS
        resu_stab_circ["X_SORTIE"] = x_sorti_FS
        resu_stab_circ["Y_ENTRE"] = np.interp(x_enter_FS, self.pente[:, 0], self.pente[:, 1])
        resu_stab_circ["Y_SORTIE"] = np.interp(x_sorti_FS, self.pente[:, 0], self.pente[:, 1])

        x_centre_FS, y_centre_FS = self.solve_circle(x_enter_FS, x_sorti_FS, R_FS)
        resu_stab_circ["CENTRE_X"] = x_centre_FS
        resu_stab_circ["CENTRE_Y"] = y_centre_FS
        resu_stab_circ["FS"] = FS_MIN

        return resu_stab_circ

    def run(self):
        """Perform the stability analysis using LEM method with circular failure surface.

        Returns:
            table object: Output table of the macro-command.
        """

        resu_stab_circ = self.search_surf_crit()
        l_FS = self.raff_maillage(resu_stab_circ)

        TABFS = CREA_TABLE(
            LISTE=(
                _F(LISTE_I=[i for i in range(len(l_FS))], PARA="NUME_RAFF"),
                _F(LISTE_R=l_FS, PARA="FS"),
                _F(LISTE_R=[resu_stab_circ["X_ENTRE"]], PARA="X_1"),
                _F(LISTE_R=[resu_stab_circ["Y_ENTRE"]], PARA="Y_1"),
                _F(LISTE_R=[resu_stab_circ["X_SORTIE"]], PARA="X_2"),
                _F(LISTE_R=[resu_stab_circ["Y_SORTIE"]], PARA="Y_2"),
                _F(LISTE_R=[resu_stab_circ["CENTRE_X"]], PARA="CENTRE_X"),
                _F(LISTE_R=[resu_stab_circ["CENTRE_Y"]], PARA="CENTRE_Y"),
                _F(LISTE_R=[resu_stab_circ["RAYON"]], PARA="RAYON"),
            )
        )

        return TABFS


class Surf_Non_Circ_Solver(LEM_Solver):
    """Solver object containing the methods for non-circular failure surface searching and FS calculation."""

    def __init__(self, parent, args):
        super(Surf_Non_Circ_Solver, self).__init__(parent, args)
        self.kw_efwa = args["ALGO_EFWA"]

    def search_surf_crit(self):
        """Search for the critical surface minimizing the FS using the EFWA algorithm.

        Returns:
            list[float], float: Geometric parameters of the critical surface and the associated FS.
        """

        efwa_optimizer = Efwa_Optimizer(self)
        etat_opti, FS_opti = efwa_optimizer.optimize()

        return etat_opti, FS_opti

    def run(self):
        """Perform the stability analysis using LEM method with circular failure surface.

        Returns:
            table object: Output table of the macro-command.
        """

        etat_opti, FS_opti = self.search_surf_crit()
        resu_stab_nc = {}

        # Restaurer les coordonnées des points sur la surface critique
        coord_Y = etat_opti.copy()
        coord_Y[0] = np.interp(etat_opti[0], self.pente[:, 0], self.pente[:, 1])
        coord_Y[-1] = np.interp(etat_opti[-1], self.pente[:, 0], self.pente[:, 1])
        coord_X = np.linspace(etat_opti[0], etat_opti[-1], self.nb_tran + 1)

        resu_stab_nc["COOR_X"] = coord_X
        resu_stab_nc["COOR_Y"] = coord_Y
        resu_stab_nc["FS"] = FS_opti

        l_FS = self.raff_maillage(resu_stab_nc)

        TABFS = CREA_TABLE(
            LISTE=(
                _F(LISTE_I=[i for i in range(len(l_FS))], PARA="NUME_RAFF"),
                _F(LISTE_R=l_FS, PARA="FS"),
                _F(LISTE_I=[i for i in range(self.nb_tran + 1)], PARA="NUME_POINT"),
                _F(LISTE_R=resu_stab_nc["COOR_X"], PARA="COOR_X"),
                _F(LISTE_R=resu_stab_nc["COOR_Y"], PARA="COOR_Y"),
            )
        )

        return TABFS


class Efwa_Optimizer:
    """Optimizer object containing the methods related to the EFWA algorithm."""

    def __init__(self, lem_solver):
        kw_efwa = lem_solver.kw_efwa
        self.lem_solver = lem_solver
        self.etat_init = kw_efwa.get("ETAT_INIT")
        self.iter_maxi = kw_efwa["ITER_MAXI"]
        self.A = kw_efwa["A"]
        self.N = kw_efwa["N"]
        self.M = kw_efwa["M"]
        self.MG = kw_efwa["MG"]
        self.SA = kw_efwa["SA"]
        self.SB = kw_efwa["SB"]
        self.nb_stab_maxi = kw_efwa["NB_STAB_MAXI"]
        self.resi_maxi = kw_efwa["CRIT_STAB"]
        self.edge_bc_ctl = kw_efwa["MARGE_PENTE"]

        self.bbas = lem_solver.x_lim_1
        self.bhaut = lem_solver.x_lim_2
        self.nb_tran = lem_solver.nb_tran
        self.nb_vari = lem_solver.nb_tran + 1

        if self.N < 1:
            UTMESS("F", "CALCSTABPENTE_20")

    def calc_bound(self, vari_etat, ivar):
        """Calculate the lower and upper limits of a stat variable.

        Args:
            vari_etat (list[float]): Array of the current stat variables.
            ivar (int): Index of the stat variable.

        Returns:
            float, float: Lower and upper limits.
        """

        if ivar == 0:
            return self.bbas[0], self.bbas[1]
        if ivar == np.size(vari_etat, 0) - 1:
            return self.bhaut[0], self.bhaut[1]

        x_enter = vari_etat[0]
        x_sorti = vari_etat[-1]
        y_enter = np.interp(x_enter, self.lem_solver.pente[:, 0], self.lem_solver.pente[:, 1])
        y_sorti = np.interp(x_sorti, self.lem_solver.pente[:, 0], self.lem_solver.pente[:, 1])

        width_tran = (x_sorti - x_enter) / self.nb_tran
        x_i = x_enter + ivar * width_tran

        # Retrouver les coordonées des points précédents
        if ivar == 1:
            x2 = x_enter
            y2 = y_enter
        else:
            x1 = x_i - 2 * width_tran
            x2 = x_i - width_tran
            if ivar == 2:
                y1 = y_enter
            else:
                y1 = vari_etat[ivar - 2]
            y2 = vari_etat[ivar - 1]

        # Calcul de la limite supérieure
        upper = ((x2 - x_i) * y_sorti + (x_i - x_sorti) * y2) / (x2 - x_sorti)

        k = (y2 - upper) / (x2 - x_i)
        for nnode in range(np.size(self.lem_solver.pente, 0)):
            x_n = self.lem_solver.pente[nnode, 0]
            y_n = self.lem_solver.pente[nnode, 1]

            if x_n < x2:
                continue
            if x_n > x_sorti:
                break
            if k * (x_n - x2) + y2 - y_n > 1e-6:
                new_upper = ((x2 - x_i) * y_n + (x_i - x_n) * y2) / (x2 - x_n)
                if new_upper < upper:
                    upper = new_upper

        # Calcul de la limite inférieure
        if ivar == 1:
            # L'angle entre la surface et le profil de pente < 45 deg afin d'avoir assez de profondeur
            for nnode in range(np.size(self.lem_solver.pente, 0)):
                if (
                    self.lem_solver.pente[nnode, 0] <= x_enter
                    and self.lem_solver.pente[nnode + 1, 0] > x_enter
                ):
                    angle_enter = np.arctan(
                        (self.lem_solver.pente[nnode + 1, 1] - self.lem_solver.pente[nnode, 1])
                        / (self.lem_solver.pente[nnode + 1, 0] - self.lem_solver.pente[nnode, 0])
                    )

            lower = y_enter + np.tan(angle_enter - np.pi / 4) * width_tran
        else:
            lower = ((x1 - x_i) * y2 + (x_i - x2) * y1) / (x1 - x2)

        # Eviter que la surface croisse avec le fond du modèle
        if lower < self.lem_solver.y_bas:
            lower = self.lem_solver.y_bas

        if lower > upper:
            UTMESS("F", "CALCSTABPENTE_21", vali=ivar)

        return lower, upper

    def check_vari_etat(self, vari_etat, ivar):
        """Verify the stat variables on the right of ivar.

        Args:
            vari_etat (list[float]): Array of the current stat variables.
            ivar (int): Index of the stat variable.

        Returns:
            list[float]: Corrected array of stat variables.
        """

        # Les bornes des variabels dépendant des points précédents,
        # on vérifie et ajuste les coor_Y des points à droite de ivar.

        for iivar in range(ivar, self.nb_vari):
            lower, upper = self.calc_bound(vari_etat, iivar)
            if vari_etat[iivar] < lower or vari_etat[iivar] >= upper:
                vari_etat[iivar] = lower + np.random.rand() * (upper - lower)
            if iivar > 0 and iivar < self.nb_vari - 1:
                # Eviter que le point est trop proche du bord
                xx = vari_etat[0] + (vari_etat[-1] - vari_etat[0]) / self.nb_tran * iivar
                yy = np.interp(xx, self.lem_solver.pente[:, 0], self.lem_solver.pente[:, 1])
                if yy - vari_etat[iivar] < self.edge_bc_ctl:
                    vari_etat[iivar] -= self.edge_bc_ctl

        return vari_etat

    def optimize(self):
        """EFWA optimisation

        Returns:
            list[float], float: Geometric parameters of the critical surface and the associated FS.
        """

        # ---------------- INITIALIZATION ------------------

        fireworks = np.zeros((self.N, self.nb_vari))

        fireworks[:, 0] = np.random.uniform(self.bbas[0], self.bbas[1], size=self.N)
        fireworks[:, -1] = np.random.uniform(self.bhaut[0], self.bhaut[1], size=self.N)

        # INITIALIZER : Générer les Y des points intermédiaires
        for ifw in range(self.N):
            for ivar in range(1, self.nb_vari - 1):
                lower, upper = self.calc_bound(fireworks[ifw, :], ivar)
                fireworks[ifw, ivar] = np.random.uniform(lower, upper)
                # Eviter que le point est trop proche du bord
                xx = (
                    fireworks[ifw, 0]
                    + (fireworks[ifw, -1] - fireworks[ifw, 0]) / self.nb_tran * ivar
                )
                yy = np.interp(xx, self.lem_solver.pente[:, 0], self.lem_solver.pente[:, 1])
                if yy - fireworks[ifw, ivar] < self.edge_bc_ctl:
                    fireworks[ifw, ivar] -= self.edge_bc_ctl

        if self.etat_init is not None:
            tab_init = self.etat_init.EXTR_TABLE().values()
            etat_init = np.zeros(self.nb_vari)
            etat_init[0] = tab_init["COOR_X"][0]
            etat_init[-1] = tab_init["COOR_X"][-1]
            etat_init[1:-1] = tab_init["COOR_Y"][1:-1]
            fireworks[0, :] = etat_init

        fit = np.zeros(self.N)
        for ifw in range(self.N):
            fit[ifw] = self.lem_solver.calc_morgenstern(fireworks[ifw, :])

        FS_opti = np.min(fit)
        id_opti = np.argmin(fit)

        # ----------------- OPTIMISATION -------------------
        nb_iter = 1
        nb_stab = 0
        epsi = 1e-8
        hist_opti = [FS_opti]
        nb_etinc = np.zeros(self.N)
        ampli = np.zeros(self.N)

        while nb_iter <= self.iter_maxi and nb_stab < self.nb_stab_maxi:
            fs_max = np.max(fit)
            denom_etin = np.sum(fs_max - fit) + epsi
            denom_ampli = np.sum(fit - FS_opti) + epsi

            for ifw in range(self.N):
                # CALC NB ETINCELLES
                Si = self.M * (fs_max - fit[ifw] + epsi) / denom_etin
                if Si < self.SA * self.M:
                    Si = self.SA * self.M
                if Si > self.SB * self.M:
                    Si = self.SB * self.M
                nb_etinc[ifw] = round(Si)

                # CALC AMPLITUDE
                denom_ampli = np.sum(fit - FS_opti) + epsi
                Ai = self.A * (fit[ifw] - FS_opti + epsi) / denom_ampli
                ampli[ifw] = Ai

            nb_etinc = nb_etinc.astype("int64")

            # EXPLOSION
            etinc = np.zeros((np.sum(nb_etinc), self.nb_vari))
            etinc_gauss = np.zeros((self.MG, self.nb_vari))

            # Générer les étincelles ordinaires et évaluer FS
            fit_etinc = np.zeros(np.sum(nb_etinc))
            for ifw in range(self.N):
                for ie in range(nb_etinc[ifw]):
                    ind = np.sum(nb_etinc[:ifw]) + ie
                    etinc[ind, :] = fireworks[ifw, :].copy()

                    for ivar in range(self.nb_vari):
                        if round(np.random.rand()) == 1:
                            dx = ampli[ifw] * np.random.uniform(-1, 1)
                            etinc[ind, ivar] += dx
                            etinc[ind, :] = self.check_vari_etat(etinc[ind, :], ivar)

                    fit_etinc[ind] = self.lem_solver.calc_morgenstern(etinc[ind, :])

            # Générer les étincelles gaussiennes et évaluer FS
            fit_gauss = np.zeros(self.MG)
            for igauss in range(self.MG):
                id_luck = int(np.random.rand() * self.N) - 1
                etinc_gauss[igauss, :] = fireworks[id_luck, :].copy()

                for ivar in range(self.nb_vari):
                    if round(np.random.rand()) == 1:
                        coef_gauss = np.random.normal(1, 1)
                        etinc_gauss[igauss, ivar] += (
                            fireworks[id_opti, ivar] - etinc_gauss[igauss, ivar]
                        ) * coef_gauss

                        etinc_gauss[igauss, :] = self.check_vari_etat(etinc_gauss[igauss, :], ivar)

                fit_gauss[igauss] = self.lem_solver.calc_morgenstern(etinc_gauss[igauss, :])

            # SELECTION PAR L'ALGORITHME DE ROULETTE
            etinc_tot = np.vstack((fireworks, etinc, etinc_gauss))
            fit_tot = np.concatenate((fit, fit_etinc, fit_gauss))

            resi = np.abs(FS_opti - np.min(fit_tot))
            if resi < self.resi_maxi:
                nb_stab += 1
            else:
                nb_stab = 0

            FS_opti = np.min(fit_tot)
            idmin = np.argmin(fit_tot)
            etat_opti = etinc_tot[idmin, :]
            fireworks[0, :] = etat_opti
            fit[0] = FS_opti

            p = (np.max(fit_tot) - fit_tot) / (np.max(fit_tot) - FS_opti)
            p_cumul = np.zeros_like(p)
            for ifw in range(np.size(p, 0)):
                p_cumul[ifw] = np.sum(p[: ifw + 1])

            shotgun = np.random.uniform(0, np.sum(p), self.N - 1)
            for ifw, val in enumerate(shotgun):
                for ie, prob in enumerate(p_cumul):
                    if prob >= val:
                        fireworks[ifw + 1, :] = etinc_tot[ie, :]
                        fit[ifw + 1] = fit_tot[ie]
                        break

            # MAJ les variables d'itération
            hist_opti.append(FS_opti)
            nb_iter += 1
            id_opti = 0

        return etat_opti, FS_opti
