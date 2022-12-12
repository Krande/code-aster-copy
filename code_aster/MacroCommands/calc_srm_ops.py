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
from ..Objects import Material, ExternalVariableTraits
from ..Commands import (
    DEFI_GROUP,
    AFFE_MATERIAU,
    STAT_NON_LINE,
    CALC_CHAMP,
    CREA_CHAMP,
    POST_RELEVE_T,
    CREA_TABLE,
    FORMULE,
    DETRUIRE,
)


def calc_srm_ops(self, **args):

    chmat = args["CHAM_MATER"]
    modele = args["MODELE"]

    # Vérification du maillage en entrée
    if not chmat.getMesh() == modele.getMesh():
        UTMESS("F", "CALCSRM_1")
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
        UTMESS("F", "CALCSRM_13")

    if prec_fin > prec:
        UTMESS("A", "CALCSRM_2", valr=prec)
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
            UTMESS("F", "CALCSRM_3")
            break

        KW_Mat_Deg = gene_deg_kwmat(CPT_UTIL, Para_Mat, dict_zone_deg, fact, args)

        __CM_DEG = AFFE_MATERIAU(AFFE=KW_Mat_Deg, **KW_Chmat_Fix)

        try:
            __RESU = STAT_NON_LINE(CHAM_MATER=__CM_DEG, **KW_Solveur_Fix)
        except:
            if num_iter == 0:
                # FS_INIT trop grand.
                UTMESS("F", "CALCSRM_4")

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
        args (_type_): Input arguments of the macro-command.

    Returns:
        str: Name of the material behaviour used in the SRM algorithm.
    """

    cpt_auto = ["MOHR_COULOMB", "DRUCK_PRAGER"]
    cpt_exis = []
    if args.get("TOUT") is not None:
        for mcf_cpt in mc_cpt:
            if mcf_cpt["RELATION"] not in cpt_auto:
                UTMESS("F", "CALCSRM_5", valk=mcf_cpt["RELATION"])
            cpt_exis.append(mcf_cpt["RELATION"])
    else:
        lgrma_deg = args["GROUP_MA"]
        for mcf_cpt in mc_cpt:
            if mcf_cpt.get("MAILLE"):
                UTMESS("F", "CALCSRM_6")

            isGPMA = mcf_cpt.get("GROUP_MA")

            cpt = mcf_cpt["RELATION"]
            if isGPMA is None:
                if cpt not in cpt_auto:
                    UTMESS("F", "CALCSRM_5", valk=cpt)
                cpt_exis.append(cpt)
            else:
                lgrma_cpt = mcf_cpt["GROUP_MA"]
                for grma in lgrma_deg:
                    lgrma_inter = check_ma_intersec(mesh, grma, lgrma_cpt)
                    if len(lgrma_inter) > 1:
                        if cpt not in cpt_auto:
                            UTMESS("F", "CALCSRM_5", valk=cpt)
                        cpt_exis.append(cpt)

    cpt_exis = list(set(cpt_exis))

    if len(cpt_exis) != 1:
        UTMESS("F", "CALCSRM_7")

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
                            UTMESS("F", "CALCSRM_8", valk=grma_deg)
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
        UTMESS("F", "CALCSRM_9", valk=grma_obj)

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
        args (dict): Arguments of CALC_SRM.

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
                UTMESS("F", "CALCSRM_10")

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
                UTMESS("F", "CALCSRM_10")

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
            UTMESS("A", "CALCSRM_11")
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
            UTMESS("F", "CALCSRM_12")

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
