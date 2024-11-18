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
"""
Compute:
  - Weibull stress
  - Weibull stress**m
  - Probability of failure

  - Strain model -> Stress to consider (6 components: XX, YY, ZZ, XY, XZ, YZ)
      - "PETIT"      -> "SIEF"
      - "PETIT_REAC" -> "SIEF"
      - "GDEF_LOG"   -> "VARI_ELGA"
      - "SIMO_MIEHE" -> To retrieve Kirchhof stresses ? Not obvious ...
                        Normally, we can do without
                        Use of HHO models instead ?

Integral must be computed on the initial configuration

_biblio : CR-T66-2017-132, Lorentz
"""
import numpy as np

from ..Cata.Syntax import _F
from ..CodeCommands import (
    CALC_TABLE,
    CALC_CHAMP,
    CREA_CHAMP,
    CALC_CHAM_ELEM,
    CREA_TABLE,
    DEFI_GROUP,
)
from ..Messages import UTMESS
from .Fracture.post_beremin_utils import get_beremin_properties, get_resu_from_deftype
from ..Objects import NonLinearResult, ExternalVariableTraits, FieldOnCellsReal, Formula


def post_beremin_ops(self, **args):
    """
    Beremin post-processor

    Arguments:
        resultat (NonLinearResult): Mechanical code_aster result
        group_ma (str): Mesh cells group on which Beremin post-treatment
            is carried out
            Only 1 mesh cells group authorized
        deformation (str): "PETIT", "PETIT_REAC", "GDEF_LOG"
        filtre_sigm (str): "SIGM_ELGA" or "SIGM_ELMOY".
        coef_mult (float): As explained in doc aster u4.81.22 (weibull) p.15

    Returns:
        Table:

            - Weibull stress
            - Weibull stress**m
            - Failure probability

    if SIGM_MAXI in command file, MED file containing:

        - where the plasticity is active : maximum of major principal
          stress
        - where the plasticity is not active : 0
          Result of type ELGA (apply ElgaFieldToSurface filter in
          Paraview to see it)
    """
    resupb = args.get("RESULTAT")
    grmapb = args.get("GROUP_MA")
    deformation = args.get("DEFORMATION")
    numvs = {"vari": args.get("LIST_NUME_VARI")}

    if deformation == "GDEF_LOG":
        numvs.update({"sief": args.get("LIST_NUME_SIEF")})

    data_resu = get_resu_from_deftype(resupb, grmapb, deformation, numvs)

    if len(tuple(elt[1] for elt in data_resu[1])) == 0:
        lentable = len(resupb.getAccessParameters()["INST"])
        liste_r = [0] * lentable
        table = CREA_TABLE(
            LISTE=(
                _F(PARA="NUME_ORDRE", LISTE_I=resupb.getAccessParameters()["NUME_ORDRE"]),
                _F(PARA="INST", LISTE_R=resupb.getAccessParameters()["INST"]),
                _F(PARA="GROUP_MA", LISTE_K=[grmapb] * lentable),
                _F(PARA="SIGMA_WEIBULL", LISTE_R=liste_r),
                _F(PARA="SIGMA_WEIBULL**M", LISTE_R=liste_r),
                _F(PARA="PROBA_WEIBULL", LISTE_R=liste_r),
            )
        )

    else:
        table = compute_pb(self, resupb, grmapb, numvs["vari"], data_resu, **args)

    return table


def compute_pb(self, resupb, grmapb, numvi, data_resu, **args):
    """
    Core function of POST_BEREMIN

    Arguments:
        resupb (NonLinearResult): Resultat input of POST_BEREMIN
        grmapb (str): Mesh cells group given in POST_BEREMIN
        numvi (tuple): Indices of EPSPEQ and INDIPLAS
        data_resu (tuple): (result to consider for Weibull stress,
                            list of time steps (nume_ordre, time step) in plasticity only,
                            list of values of maximum of plasticity, if strictly greater than 0)

    Returns:
        Table: Final table with Weibull stress and failure probability
    """
    resanpb = args.get("SIGM_MAXI")

    dwb = get_beremin_properties(resupb, grmapb)

    mawbrest = make_plasticity_groups(data_resu, grmapb, numvi[0], dwb[grmapb]["SEUIL_EPSP_CUMU"])

    if args.get("FILTRE_SIGM") == "SIGM_ELGA":

        rsieq = CALC_CHAMP(
            RESULTAT=data_resu[0],
            GROUP_MA="mgrplasfull",
            INST=tuple(elt[1] for elt in data_resu[1]),
            CRITERES="SIEQ_ELGA",
        )

        (sigw, resimpr) = tps_maxsigm(
            rsieq,
            data_resu[1],
            sig1plasac(data_resu[0], rsieq, numvi, dwb, resupb, grmapb, data_resu[1]),
            resanpb,
            dwb[grmapb]["M"],
            args.get("MAXI_TPS"),
        )

    elif args.get("FILTRE_SIGM") == "SIGM_ELMOY":

        (sigw, resimpr) = compute_sigm_elmoy(
            data_resu, numvi, dwb, resupb, grmapb, resanpb, args.get("MAXI_TPS")
        )

    else:
        assert False

    if resanpb is not None:
        self.register_result(resimpr, resanpb)

    table = compute_beremin_integral(args.get("COEF_MULT"), sigw, dwb, grmapb, resupb)

    itlist = [elt[0] for elt in data_resu[1]]

    mawbrest = DEFI_GROUP(
        reuse=mawbrest,
        MAILLAGE=mawbrest,
        DETR_GROUP_NO=_F(
            NOM=tuple(
                ["ngrmapb"]
                + [f"vale_{iteration}" for iteration in itlist]
                + [f"ngrplas_{iteration}" for iteration in itlist]
            )
        ),
    )

    mawbrest = DEFI_GROUP(
        reuse=mawbrest,
        MAILLAGE=mawbrest,
        DETR_GROUP_MA=_F(
            NOM=tuple(
                [f"mgrplas_{iteration}" for iteration in itlist]
                + [f"ngrplas_{iteration}" for iteration in itlist]
                + ["mgrplasfull"]
            )
        ),
    )

    mawbrest.resetDependencies()

    return table


def compute_sigm_elmoy(data_resu, numvi, dwb, resupb, grmapb, resanpb, mtpb):
    """
    Compute POST_BEREMIN/SIGM_ELMOY

    Arguments:
        data_resu (tuple): (result to consider for Weibull stress,
                            list of time steps (nume_ordre, time step) in plasticity only,
                            list of values of maximum of plasticity, if strictly greater than 0)
        numvi (tuple): Indices of EPSPEQ and INDIPLAS
        dwb (dict): Beremin parameters
        resupb (NonLinearResult): Resultat input of POST_BEREMIN
        grmapb (str): Mesh cells group given in POST_BEREMIN
        resanpb (NonLinearResult): Name of auxiliary result
        mtpb (str): {"OUI", "NON"} Value of keyword MAXI_TPS

    Returns:
        (NonLinearResult, NonLinearResult): (
            ELGA_DEPL_R filled by PRIN_3,
            informative MED file: where the plasticity is active: maximum of major principal stress,
                                            )
    """
    rmelmoy = CALC_CHAMP(
        RESULTAT=data_resu[0],
        GROUP_MA="mgrplasfull",
        INST=tuple(elt[1] for elt in data_resu[1]),
        CONTRAINTE="SIMY_ELGA",
    )

    relmoysief = NonLinearResult()
    relmoysief.allocate(rmelmoy.getNumberOfIndexes())

    for nume_inst in rmelmoy.getIndexes():
        relmoysief.setField(rmelmoy.getField("SIMY_ELGA", nume_inst), "SIEF_ELGA", nume_inst)
        for field in ("COMPORTEMENT", "VARI_ELGA"):
            relmoysief.setField(data_resu[0].getField(field, nume_inst), field, nume_inst)
            relmoysief.setModel(rmelmoy.getModel(nume_inst), nume_inst)
            relmoysief.setTime(rmelmoy.getTime(nume_inst), nume_inst)

    rsieq = CALC_CHAMP(
        RESULTAT=relmoysief,
        GROUP_MA="mgrplasfull",
        INST=tuple(elt[1] for elt in data_resu[1]),
        CRITERES="SIEQ_ELGA",
    )

    (sigw, resimpr) = tps_maxsigm(
        rsieq,
        data_resu[1],
        sig1plasac(relmoysief, rsieq, numvi, dwb, resupb, grmapb, data_resu[1]),
        resanpb,
        dwb[grmapb]["M"],
        mtpb,
    )

    return (sigw, resimpr)


def sigma1(rsieq, nume_inst, dwb, reswbrest, grwb):
    """
    Major principal stress

    Arguments:
        rsieq (NonLinearResult): SIEQ_ELGA field
        nume_inst (int): Rank of instant
        stress depends of temperature
        dwb (dict): Beremin parameters
        reswbrest (NonLinearResult): result restricted where plasticity is
            greater than 0
        grwb (str): name of mesh cells group given in POST_BEREMIN

    Returns:
        FieldOnCells: ELGA_DEPL_R filled by PRIN_3
    """
    modele = reswbrest.getModel()
    if modele.getMesh().hasGroupOfCells(f"mgrplas_{nume_inst}"):

        if "SIGM_CNV" not in dwb[grwb]:

            sg1 = FieldOnCellsReal(modele, "ELGA", "DEPL_R")
            sg1.setValues(
                (
                    rsieq.getField("SIEQ_ELGA", nume_inst)
                    .asPhysicalQuantity("DEPL_R", {"PRIN_3": "DX"})
                    .toSimpleFieldOnCells()
                )
                .restrict(["DX"], [f"mgrplas_{nume_inst}"])
                .toFieldOnCells(modele.getFiniteElementDescriptor(), "TOU_INI_ELGA", "")
                .getValues()
            )

        else:

            sg1 = sigma1_f(rsieq, nume_inst, dwb, reswbrest, grwb)

    else:
        sg1 = FieldOnCellsReal(modele, "ELGA", "DEPL_R")
        sg1.setValues(0)

    if "CRIT_SIGM" in dwb[grwb]:

        def crit_sigm(sigma):
            return max(sigma - dwb[grwb]["CRIT_SIGM"], 0.0)

        sg1 = sg1.transform(crit_sigm)

    return sg1


def sigma1_f(rsieq, nume_inst, dwb, reswbrest, grwb):
    """
    Major principal stress (WEIBULL_FO)

    Arguments:
        rsieq (NonLinearResult): SIEQ_ELGA field
        nume_inst (int): Rank of instant
        stress depends of temperature
        dwb (dict): Beremin parameters
        reswbrest (NonLinearResult): result restricted where plasticity is
            greater than 0
        grwb (str): name of mesh cells group given in POST_BEREMIN

    Returns:
        FieldOnCells: ELGA_DEPL_R filled by PRIN_3
    """
    grmacalc = f"mgrplas_{nume_inst}"
    modele = reswbrest.getModel()

    chf = CREA_CHAMP(
        AFFE=_F(NOM_CMP="X1", GROUP_MA=grmacalc, VALE_F=dwb[grwb]["SIGM_REFE"]),
        TYPE_CHAM="NOEU_NEUT_F",
        MODELE=modele,
        OPERATION="AFFE",
    )

    for curiter in reswbrest.getMaterialField().getExtStateVariablesOnMeshEntities():
        extvari = curiter[0]
        nom_varc = ExternalVariableTraits.getExternVarTypeStr(extvari.getType())
        if nom_varc == "TEMP":
            inputfield = extvari.getField()
            evolparam = extvari.getEvolutionParameter()

            if inputfield:
                chamchpar = inputfield
            if evolparam:
                chamchpar = evolparam.getTransientResult().getField("TEMP", nume_inst)

    sgrefeno = CREA_CHAMP(
        OPERATION="EVAL",
        TYPE_CHAM="NOEU_NEUT_R",
        CHAM_F=chf,
        CHAM_PARA=CREA_CHAMP(
            OPERATION="ASSE",
            TYPE_CHAM="NOEU_TEMP_R",
            MODELE=modele,
            ASSE=_F(GROUP_MA=grmacalc, CHAM_GD=chamchpar, NOM_CMP="TEMP", NOM_CMP_RESU="TEMP"),
        ),
    )

    chno = CREA_CHAMP(
        OPERATION="ASSE",
        TYPE_CHAM="NOEU_DEPL_R",
        MODELE=modele,
        ASSE=_F(GROUP_MA=grmacalc, CHAM_GD=sgrefeno, NOM_CMP="X1", NOM_CMP_RESU="DX"),
    )

    sgrefega = CREA_CHAMP(
        TYPE_CHAM="ELGA_DEPL_R", OPERATION="DISC", MODELE=modele, PROL_ZERO="OUI", CHAM_GD=chno
    )

    sqsursgrefe = CREA_CHAMP(
        OPERATION="ASSE",
        TYPE_CHAM="ELGA_DEPL_R",
        MODELE=modele,
        PROL_ZERO="OUI",
        ASSE=(
            _F(
                GROUP_MA=grmacalc,
                CHAM_GD=rsieq.getField("SIEQ_ELGA", nume_inst),
                NOM_CMP="PRIN_3",
                NOM_CMP_RESU="DX",
            ),
            _F(GROUP_MA=grmacalc, CHAM_GD=sgrefega, NOM_CMP="DX", NOM_CMP_RESU="DY"),
        ),
    )

    rdivaux = NonLinearResult()
    rdivaux.allocate(1)
    rdivaux.setField(sqsursgrefe, "DEPL_ELGA", 0)
    rdivaux.setModel(modele, 0)

    formule = Formula()
    formule.setExpression("DX/DY*sigm_cnv")
    formule.setVariables(["DX", "DY"])
    formule.setContext({"sigm_cnv": dwb[grwb]["SIGM_CNV"]})

    rdiv1 = CALC_CHAMP(
        RESULTAT=rdivaux,
        GROUP_MA=grmacalc,
        CHAM_UTIL=_F(NOM_CHAM="DEPL_ELGA", FORMULE=formule, NUME_CHAM_RESU=1),
    )

    return rdiv1.getField("UT01_ELGA", 0).asPhysicalQuantity("DEPL_R", {"X1": "DX"})


def sig1plasac(resultat, rsieq, numvi, dwb, reswbrest, grmapb, l_instplas):
    """
    Major principal stress where the plasticity is active

    Arguments:
        resultat (NonLinearResult): Result to consider to compute Weibull
            stress
        rsieq (NonLinearResult): SIEQ_ELGA field
        numvi (tuple): Indices of internal variables EPSPEQ and INDIPLAS
        dwb (dict): Weibull parameters
        reswbrest (NonLinearResult): Result where plasticity is strictly
            positive
        grmapb (str): Mesh cells group given in POST_BEREMIN
        l_instplas (list): List of time steps (order, time step) where there is plasticity

    Returns:
        NonLinearResult: ELGA_DEPL_R filled by PRIN_3
    """
    modele = resultat.getModel()
    if grmapb not in dwb:
        UTMESS("F", "RUPTURE1_88", valk=(grmapb))

    maxsig = NonLinearResult()
    maxsig.allocate(rsieq.getNumberOfIndexes())
    seuil = dwb[grmapb]["SEUIL_EPSP_CUMU"]
    indice = 0
    for nume_inst in rsieq.getAccessParameters()["NUME_ORDRE"]:
        inst = rsieq.getTime(nume_inst)

        if inst in [elt[1] for elt in l_instplas]:

            tronque = FieldOnCellsReal(modele, "ELGA", "DEPL_R")
            tronque.setValues(
                [
                    np.sign(max(_v1 - seuil, 0) * _v2)
                    for (_v1, _v2) in zip(
                        (
                            resultat.getField("VARI_ELGA", nume_inst)
                            .asPhysicalQuantity("DEPL_R", {f"V{numvi[0]}": "DX"})
                            .toSimpleFieldOnCells()
                        )
                        .restrict(["DX"])
                        .toFieldOnCells(modele.getFiniteElementDescriptor(), "TOU_INI_ELGA", "")
                        .getValues(),
                        (
                            resultat.getField("VARI_ELGA", nume_inst)
                            .asPhysicalQuantity("DEPL_R", {f"V{numvi[1]}": "DX"})
                            .toSimpleFieldOnCells()
                        )
                        .restrict(["DX"])
                        .toFieldOnCells(modele.getFiniteElementDescriptor(), "TOU_INI_ELGA", "")
                        .getValues(),
                    )
                ]
            )

            sigtyp = FieldOnCellsReal(modele, "ELGA", "DEPL_R")

            sigtyp.setValues(
                [
                    vxx * vyy
                    for (vxx, vyy) in zip(
                        (sigma1(rsieq, nume_inst, dwb, reswbrest, grmapb).toSimpleFieldOnCells())
                        .restrict(["DX"])
                        .toFieldOnCells(modele.getFiniteElementDescriptor(), "TOU_INI_ELGA", "")
                        .getValues(),
                        (tronque.toSimpleFieldOnCells())
                        .restrict(["DX"])
                        .toFieldOnCells(modele.getFiniteElementDescriptor(), "TOU_INI_ELGA", "")
                        .getValues(),
                    )
                ]
            )

            maxsig.setField(sigtyp, "DEPL_ELGA", indice)
            maxsig.setTime(resultat.getTime(nume_inst), indice)
            indice = indice + 1

    maxsig.setModel(modele)

    return maxsig


def tps_maxsigm(rsieq, l_instplas, maxsig, resanpb, bere_m, mtpb):
    """
    Compute temporal maximum for each Gauss point and elevation to power m

    Arguments:
        rsieq (NonLinearResult): SIEQ_ELGA field
        l_instplas (list): List of time steps (order, time step) where there is plasticity
        maxsig (NonLinearResult): ELGA_DEPL_R filled by PRIN_3
        resanpb (NonLinearResult): Name of auxiliary result
        bere_m (float): Value of Beremin parameter M
        mtpb (str): {"OUI", "NON"} Value of keyword MAXI_TPS

    Returns:
        (NonLinearResult, NonLinearResult): (
            ELGA_DEPL_R filled by PRIN_3,
            informative MED file: where the plasticity is active: maximum of major principal stress,
                                            )
    """
    if mtpb == "OUI":

        (sigw, resimpr) = maxi_tps_oui(rsieq, l_instplas, maxsig, resanpb, bere_m)

    elif mtpb == "NON":

        (sigw, resimpr) = maxi_tps_non(rsieq, l_instplas, maxsig, resanpb, bere_m)

    sigw.setModel(maxsig.getModel())

    return (sigw, resimpr)


def maxi_tps_oui(rsieq, l_instplas, maxsig, resanpb, bere_m):
    """
    Compute temporal maximum for each Gauss point and elevation to power m

    Arguments:
        rsieq (NonLinearResult): SIEQ_ELGA field
        l_instplas (list): List of time steps (order, time step) where there is plasticity
        maxsig (NonLinearResult): ELGA_DEPL_R filled by PRIN_3
        resanpb (NonLinearResult): Name of auxiliary result
        bere_m (float): Value of Beremin parameter M

    Returns:
        (NonLinearResult, NonLinearResult): (
            ELGA_DEPL_R filled by PRIN_3,
            informative MED file: where the plasticity is active: maximum of major principal stress,
                                            )

    if SIGM_MAXI in command file, MED file containing:

        - where the plasticity is active: maximum of major principal stress
        - where the plasticity is not active: 0
    """
    resimpr = None
    if resanpb is not None:
        if resanpb.is_typco():
            resimpr = NonLinearResult()
            resimpr.allocate(rsieq.getNumberOfIndexes())

    sigw = NonLinearResult()
    sigw.allocate(rsieq.getNumberOfIndexes())

    def puiss_m(valsixx):
        return valsixx**bere_m

    indice = 0
    for nume_inst, inst in enumerate(rsieq.getAccessParameters()["INST"]):

        if inst in [elt[1] for elt in l_instplas]:

            if inst == maxsig.getTime(0):
                chmaxsig = maxsig.getField("DEPL_ELGA", 0)
                maxsig_r = NonLinearResult()
                maxsig_r.allocate(2)
                maxsig_r.setField(chmaxsig, "DEPL_ELGA", 0)
            elif inst > maxsig.getTime(0):

                maxsig_r.setField(maxsig.getField("DEPL_ELGA", indice), "DEPL_ELGA", 1)

                chmaxsig = CREA_CHAMP(
                    TYPE_CHAM="ELGA_DEPL_R",
                    OPERATION="EXTR",
                    RESULTAT=maxsig_r,
                    NOM_CHAM="DEPL_ELGA",
                    TYPE_MAXI="MAXI",
                    TYPE_RESU="VALE",
                    NUME_ORDRE=(0, 1),
                )

                maxsig_r.setField(chmaxsig, "DEPL_ELGA", 0)

            chsixxm = chmaxsig.transform(puiss_m)

            if resanpb is not None:
                resimpr.setField(chsixxm, "DEPL_ELGA", nume_inst)
                resimpr.setTime(inst, nume_inst)

            indice = indice + 1

        else:

            if resanpb is not None:

                chsixxm = FieldOnCellsReal(chmaxsig.getModel(), "ELGA", "DEPL_R")
                chsixxm.setValues(0)
                resimpr.setField(chsixxm, "DEPL_ELGA", nume_inst)
                resimpr.setTime(inst, nume_inst)

        sigw.setField(chsixxm, "DEPL_ELGA", nume_inst)
        sigw.setTime(inst, nume_inst)

    return (sigw, resimpr)


def maxi_tps_non(rsieq, l_instplas, maxsig, resanpb, bere_m):
    """
    Compute instantaneous maximum for each Gauss point and elevation to power m

    Arguments:
        rsieq (NonLinearResult): SIEQ_ELGA field
        l_instplas (list): List of time steps (order, time step) where there is plasticity
        maxsig (NonLinearResult): ELGA_DEPL_R filled by PRIN_3
        resanpb (NonLinearResult): Name of auxiliary result
        bere_m (float): Value of Beremin parameter M

    Returns:
        (NonLinearResult, NonLinearResult): (
            ELGA_DEPL_R filled by PRIN_3,
            informative MED file: where the plasticity is active: maximum of major principal stress,
                                            )

    if SIGM_MAXI in command file, MED file containing:

        - where the plasticity is active: maximum of major principal stress
        - where the plasticity is not active: 0
    """
    resimpr = None
    if resanpb is not None:
        if resanpb.is_typco():
            resimpr = NonLinearResult()
            resimpr.allocate(rsieq.getNumberOfIndexes())

    sigw = NonLinearResult()
    sigw.allocate(rsieq.getNumberOfIndexes())

    def puiss_m(valsixx):
        return valsixx**bere_m

    indice = 0
    for nume_inst, inst in enumerate(rsieq.getAccessParameters()["INST"]):

        if inst in [elt[1] for elt in l_instplas]:

            if inst == maxsig.getTime(0):

                chmaxsig = maxsig.getField("DEPL_ELGA", 0)

            elif inst > maxsig.getTime(0):

                chmaxsig = maxsig.getField("DEPL_ELGA", indice)

            chsixxm = chmaxsig.transform(puiss_m)

            if resanpb is not None:
                resimpr.setField(chsixxm, "DEPL_ELGA", nume_inst)
                resimpr.setTime(inst, nume_inst)

            indice = indice + 1

        else:

            if resanpb is not None:

                chsixxm = FieldOnCellsReal(chmaxsig.getModel(), "ELGA", "DEPL_R")
                chsixxm.setValues(0)
                resimpr.setField(chsixxm, "DEPL_ELGA", nume_inst)
                resimpr.setTime(inst, nume_inst)

        sigw.setField(chsixxm, "DEPL_ELGA", nume_inst)
        sigw.setTime(inst, nume_inst)

    return (sigw, resimpr)


def compute_beremin_integral(coefmultpb, sigw, dwb, grmapb, resupb):
    """
    Integral computation for Weibull stress
    Produces also sigw**m and probability of failure

    Arguments:
        coefmultpb (float): Value of COEF_MULT in POST_BEREMIN
        sigw (NonLinearResult): ELGA_DEPL_R filled by PRIN_3
        dwb (dict): POST_BEREMIN keywords
        grmapb (str): Mesh cells group given in POST_BEREMIN
        resupb (NonLinearResult): Resultat input of POST_BEREMIN

    Returns:
        Table: POST_BEREMIN final table
    """
    sigwinst = sigw.getAccessParameters()["NUME_ORDRE"]
    poids = CALC_CHAM_ELEM(MODELE=sigw.getModel(), OPTION="COOR_ELGA").getValuesWithDescription(
        "W", ["mgrplasfull"]
    )[0]
    valinte_sixx = [
        np.sum(
            np.array(
                (
                    [
                        elt_poids * elt_sigw
                        for (elt_poids, elt_sigw) in zip(
                            poids,
                            sigw.getField("DEPL_ELGA", nume_inst).getValuesWithDescription(
                                "DX", ["mgrplasfull"]
                            )[0],
                        )
                    ]
                )
            )
        )
        for nume_inst in sigwinst
    ]

    sigwaux = CREA_TABLE(
        LISTE=(_F(PARA="NUME_ORDRE", LISTE_I=sigwinst), _F(PARA="INTE_SIXX", LISTE_R=valinte_sixx))
    )

    if "SIGM_CNV" in dwb[grmapb]:
        sigma_u = dwb[grmapb]["SIGM_CNV"]
    else:
        sigma_u = dwb[grmapb]["SIGM_REFE"]

    d_form = defi_formule(dwb, grmapb, coefmultpb, sigma_u)

    sigwaux = CALC_TABLE(
        TABLE=sigwaux,
        reuse=sigwaux,
        ACTION=(
            _F(OPERATION="OPER", FORMULE=d_form["intfin"], NOM_PARA="SIGMA_WEIBULL"),
            _F(OPERATION="OPER", FORMULE=d_form["sigwm"], NOM_PARA="SIGMA_WEIBULL**M"),
            _F(OPERATION="OPER", FORMULE=d_form["probaw"], NOM_PARA="PROBA_WEIBULL"),
        ),
    )

    tab_aux = {
        clef: (sigwaux.EXTR_TABLE()).values()[clef]
        for clef in ("SIGMA_WEIBULL", "SIGMA_WEIBULL**M", "PROBA_WEIBULL")
    }

    lentable = len(resupb.getAccessParameters()["INST"])
    liste_zero = [0] * (lentable - len(tab_aux["SIGMA_WEIBULL"]))
    tab_final = CREA_TABLE(
        LISTE=(
            _F(PARA="NUME_ORDRE", LISTE_I=resupb.getAccessParameters()["NUME_ORDRE"]),
            _F(PARA="INST", LISTE_R=resupb.getAccessParameters()["INST"]),
            _F(PARA="GROUP_MA", LISTE_K=[grmapb] * lentable),
            _F(PARA="SIGMA_WEIBULL", LISTE_R=liste_zero + tab_aux["SIGMA_WEIBULL"]),
            _F(PARA="SIGMA_WEIBULL**M", LISTE_R=liste_zero + tab_aux["SIGMA_WEIBULL**M"]),
            _F(PARA="PROBA_WEIBULL", LISTE_R=liste_zero + tab_aux["PROBA_WEIBULL"]),
        )
    )

    return tab_final


def defi_formule(dwb, grmapb, coefmultpb, sigma_u):
    """
    Define final formulae.

    Arguments:
        dwb (dict): POST_BEREMIN keywords
        grmapb (str): Mesh cells group given in POST_BEREMIN
        coefmultpb (float): Value of COEF_MULT in POST_BEREMIN
        sigma_u (float): fracture stress

    Returns:
        Table: POST_BEREMIN final table

    """
    t_clef = ("intfin", "sigwm", "probaw")
    t_expr = (
        "(COEF_MULT/v_0*INTE_SIXX)**(1/bere_m)",
        "SIGMA_WEIBULL**bere_m",
        "1-exp(-SIGMA_WEIBULL**bere_m/sigma_u**bere_m)",
    )
    t_var = (["INTE_SIXX"], ["SIGMA_WEIBULL"], ["SIGMA_WEIBULL"])
    t_context = (
        {"COEF_MULT": coefmultpb, "v_0": dwb[grmapb]["VOLU_REFE"], "bere_m": dwb[grmapb]["M"]},
        {"bere_m": dwb[grmapb]["M"]},
        {"sigma_u": sigma_u, "bere_m": dwb[grmapb]["M"], "exp": np.exp},
    )

    d_form = {clef: Formula() for clef in t_clef}

    for (clef, expression, variable, context) in zip(t_clef, t_expr, t_var, t_context):
        d_form[clef].setExpression(expression)
        d_form[clef].setVariables(variable)
        d_form[clef].setContext(context)

    return d_form


def make_plasticity_groups(data_resu, grmapb, numv1, seuil):
    """
    Construct groups of plasticity for each time step when there is plasticity

    Arguments:
        data_resu (tuple): (result to consider for Weibull stress,
                            list of time steps (nume_ordre, time step) in plasticity only,
                            list of values of maximum of plasticity, if strictly greater than 0)
        numv1 (int): Indice of internal variable EPSPEQ
        seuil (float): value of SEUIL_EPSP_CUMU

    Returns:
        mawbrest (Mesh): Initial mesh with groups useful for BEREMIN computation
    """
    mawbrest = data_resu[0].getModel().getMesh()
    mawbrest.setGroupOfNodes("ngrmapb", mawbrest.getNodesFromCells(grmapb))

    liter = [elt[0] for elt in data_resu[1]]
    dval = {}
    l_inte_vale = []

    for indice, iteration in enumerate(liter):

        if data_resu[2][indice] < seuil:
            dval[f"min{indice}"] = 0
            dval[f"max{indice}"] = data_resu[2][indice]
        else:
            dval[f"min{indice}"] = seuil
            dval[f"max{indice}"] = data_resu[2][indice]

    for (indice, iteration) in enumerate(liter):
        cham_gd = data_resu[0].getField("VARI_ELGA", iteration).toFieldOnNodes()
        np_cham_gd = np.array(cham_gd.getValuesWithDescription(f"V{numv1}", [])[0])

        if np.any((np_cham_gd > dval[f"min{indice}"]) & (np_cham_gd < dval[f"max{indice}"])):
            l_inte_vale.append(
                _F(
                    NOM=f"vale_{iteration}",
                    OPTION="INTERVALLE_VALE",
                    CHAM_GD=cham_gd,
                    NOM_CMP=f"V{numv1}",
                    VALE=(dval[f"min{indice}"], dval[f"max{indice}"]),
                )
            )

    mawbrest = DEFI_GROUP(reuse=mawbrest, MAILLAGE=mawbrest, CREA_GROUP_NO=tuple(l_inte_vale))

    d_gr = {"l_ngrplas": [], "l_mgrplas": [], "l_union": []}
    for iteration in liter:
        if mawbrest.hasGroupOfNodes(f"vale_{iteration}"):
            d_gr["l_ngrplas"].append(
                _F(NOM=f"ngrplas_{iteration}", INTERSEC=(f"vale_{iteration}", "ngrmapb"))
            )
            d_gr["l_mgrplas"].append(
                _F(
                    NOM=f"mgrplas_{iteration}",
                    OPTION="APPUI",
                    TYPE_MAILLE="{}D".format(mawbrest.getDimension()),
                    GROUP_NO=f"ngrplas_{iteration}",
                    TYPE_APPUI="AU_MOINS_UN",
                )
            )
            d_gr["l_union"].append(f"mgrplas_{iteration}")

    mawbrest = DEFI_GROUP(reuse=mawbrest, MAILLAGE=mawbrest, CREA_GROUP_NO=tuple(d_gr["l_ngrplas"]))

    mawbrest = DEFI_GROUP(
        reuse=mawbrest,
        MAILLAGE=mawbrest,
        CREA_GROUP_MA=tuple(d_gr["l_mgrplas"])
        + tuple(
            [
                _F(
                    NOM="mgrplasfull",
                    TYPE_MAILLE="{}D".format(mawbrest.getDimension()),
                    UNION=tuple(d_gr["l_union"]),
                )
            ]
        ),
    )

    return mawbrest
