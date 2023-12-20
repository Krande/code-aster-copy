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
"""
Compute:
  - Weibull stress
  - Weibull stress**m
  - Probability of failure

  - Strain model -> Stress to consider (6 components: XX, YY, ZZ, XY, XZ, YZ)
      - "PETIT"      -> "SIEF"
      - "PETIT_REAC" -> "SIEF"
      - "GDEF_LOG"   -> "VARI_ELGA_NOMME"
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
    FORMULE,
    POST_ELEM,
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
    resanpb = args.get("SIGM_MAXI")
    fspb = args.get("FILTRE_SIGM")

    (reswbrest, numv1v2, mclinst, l_epspmax) = get_resu_from_deftype(
        resupb, grmapb, args.get("DEFORMATION")
    )

    linstplasac = tuple(elt[2] for elt in mclinst)

    if len(linstplasac) == 0:
        lentable = len(resupb.getAccessParameters()["INST"])
        liste_r = [0] * lentable
        table = CREA_TABLE(
            LISTE=(
                _F(PARA="NUME_ORDRE", LISTE_I=resupb.getAccessParameters()["NUME_ORDRE"]),
                _F(PARA="INST", LISTE_R=resupb.getAccessParameters()["INST"]),
                _F(PARA="GROUP_MA", LISTE_K=[grmapb[0]] * lentable),
                _F(PARA="SIGMA_WEIBULL", LISTE_R=liste_r),
                _F(PARA="SIGMA_WEIBULL**M", LISTE_R=liste_r),
                _F(PARA="PROBA_WEIBULL", LISTE_R=liste_r),
            )
        )

    else:

        dwb = get_beremin_properties(resupb, grmapb)

        mawbrest = reswbrest.getModel().getMesh()

        mawbrest = DEFI_GROUP(
            reuse=mawbrest, MAILLAGE=mawbrest, CREA_GROUP_NO=_F(GROUP_MA=grmapb[0], NOM="ngrmapb")
        )

        make_plasticity_groups(
            reswbrest, numv1v2, mclinst, dwb[grmapb[0]]["SEUIL_EPSP_CUMU"], l_epspmax
        )

        if fspb == "SIGM_ELGA":

            rsieq = CALC_CHAMP(
                RESULTAT=reswbrest, GROUP_MA="mgrplasfull", INST=linstplasac, CRITERES="SIEQ_ELGA"
            )

            (sigw, resimpr) = tps_maxsigm(
                rsieq,
                mclinst,
                sig1plasac(reswbrest, rsieq, numv1v2, dwb, resupb, grmapb, mclinst),
                resanpb,
                dwb[grmapb[0]]["M"],
            )

        elif fspb == "SIGM_ELMOY":

            rmelmoy = CALC_CHAMP(
                RESULTAT=reswbrest, GROUP_MA="mgrplasfull", INST=linstplasac, CONTRAINTE="SIMY_ELGA"
            )

            relmoysief = NonLinearResult()
            relmoysief.allocate(rmelmoy.getNumberOfIndexes())

            for nume_inst in rmelmoy.getIndexes():
                relmoysief.setField(
                    rmelmoy.getField("SIMY_ELGA", nume_inst), "SIEF_ELGA", nume_inst
                )
                for field in ("COMPORTEMENT", "VARI_ELGA"):
                    relmoysief.setField(reswbrest.getField(field, nume_inst), field, nume_inst)
                relmoysief.setModel(rmelmoy.getModel(nume_inst), nume_inst)
                relmoysief.setTime(rmelmoy.getTime(nume_inst), nume_inst)

            rsieq = CALC_CHAMP(
                RESULTAT=relmoysief, GROUP_MA="mgrplasfull", INST=linstplasac, CRITERES="SIEQ_ELGA"
            )

            (sigw, resimpr) = tps_maxsigm(
                rsieq,
                mclinst,
                sig1plasac(relmoysief, rsieq, numv1v2, dwb, resupb, grmapb, mclinst),
                resanpb,
                dwb[grmapb[0]]["M"],
            )

        else:
            assert False

        if resanpb is not None:
            self.register_result(resimpr, resanpb)

        table = compute_beremin_integral(
            reswbrest.getModel(), args.get("COEF_MULT"), sigw, dwb, grmapb, resupb
        )

        itlist = [elt[0] for elt in mclinst]

        mawbrest = reswbrest.getModel().getMesh()
        mawbrest = DEFI_GROUP(
            reuse=mawbrest,
            MAILLAGE=mawbrest,
            DETR_GROUP_NO=_F(
                NOM=tuple(
                    ["ngrmapb"]
                    + ["vale_{}".format(iteration) for iteration in itlist]
                    + ["ngrplas_{}".format(iteration) for iteration in itlist]
                )
            ),
        )

        mawbrest = DEFI_GROUP(
            reuse=mawbrest,
            MAILLAGE=mawbrest,
            DETR_GROUP_MA=_F(
                NOM=tuple(
                    ["mgrplas_{}".format(iteration) for iteration in itlist]
                    + ["ngrplas_{}".format(iteration) for iteration in itlist]
                    + ["mgrplasfull"]
                )
            ),
        )

    mawbrest.resetDependencies()

    return table


def indic_plasac(_v1, _v2, seuil):
    """
    Characteristic function of active plasticity

    Arguments:
        _v1 (float): Cumulated plastic strain
        _v2 (float): Plasticity indicator
        seuil (float): Value of SEUIL_EPSP_CUMU

    Returns:
        Float: Characteristic function of active plasticity
    """
    vale = 0
    if (_v1 >= seuil) and (_v2 > 0.99):
        vale = 1
    return vale


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
        FieldOnCells: ELGA_NEUT_R filled by PRIN_3
    """
    if not "SIGM_CNV" in dwb[grwb]:

        grmacalc = "mgrplas_{}".format(nume_inst)
        modele = reswbrest.getModel()

        __sg1neut = CREA_CHAMP(
            OPERATION="ASSE",
            TYPE_CHAM="ELGA_NEUT_R",
            MODELE=modele,
            PROL_ZERO="OUI",
            ASSE=_F(
                GROUP_MA=grmacalc,
                CHAM_GD=rsieq.getField("SIEQ_ELGA", nume_inst),
                NOM_CMP="PRIN_3",
                NOM_CMP_RESU="X1",
            ),
        )

    else:

        __sg1neut = sigma1_f(rsieq, nume_inst, dwb, reswbrest, grwb)

    return __sg1neut


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
        FieldOnCells: ELGA_NEUT_R filled by PRIN_3
    """
    grmacalc = "mgrplas_{}".format(nume_inst)
    modele = reswbrest.getModel()

    chf = CREA_CHAMP(
        AFFE=_F(NOM_CMP="X1", GROUP_MA=grmacalc, VALE_F=dwb[grwb]["SIGM_REFE"]),
        TYPE_CHAM="NOEU_NEUT_F",
        MODELE=modele,
        OPERATION="AFFE",
    )

    chmat = reswbrest.getMaterialField()
    extvariaffe = chmat.getExtStateVariablesOnMeshEntities()
    for curiter in extvariaffe:
        extvari = curiter[0]
        nom_varc = ExternalVariableTraits.getExternVarTypeStr(extvari.getType())
        if nom_varc == "TEMP":
            inputfield = extvari.getField()
            evolparam = extvari.getEvolutionParameter()

            if inputfield:
                chamchpar = inputfield
            if evolparam:
                chamchpar = evolparam.getTransientResult().getField("TEMP", nume_inst)

    chpar = CREA_CHAMP(
        OPERATION="ASSE",
        TYPE_CHAM="NOEU_TEMP_R",
        MODELE=modele,
        ASSE=_F(GROUP_MA=grmacalc, CHAM_GD=chamchpar, NOM_CMP="TEMP", NOM_CMP_RESU="TEMP"),
    )

    sgrefeno = CREA_CHAMP(OPERATION="EVAL", TYPE_CHAM="NOEU_NEUT_R", CHAM_F=chf, CHAM_PARA=chpar)

    sgrefega = CREA_CHAMP(
        TYPE_CHAM="ELGA_NEUT_R", OPERATION="DISC", MODELE=modele, PROL_ZERO="OUI", CHAM_GD=sgrefeno
    )

    sqsursgrefe = CREA_CHAMP(
        OPERATION="ASSE",
        TYPE_CHAM="ELGA_NEUT_R",
        MODELE=modele,
        PROL_ZERO="OUI",
        ASSE=(
            _F(
                GROUP_MA=grmacalc,
                CHAM_GD=rsieq.getField("SIEQ_ELGA", nume_inst),
                NOM_CMP="PRIN_3",
                NOM_CMP_RESU="X1",
            ),
            _F(GROUP_MA=grmacalc, CHAM_GD=sgrefega, NOM_CMP="X1", NOM_CMP_RESU="X2"),
        ),
    )

    sq2 = CREA_CHAMP(
        OPERATION="ASSE",
        TYPE_CHAM="ELGA_SIEF_R",
        MODELE=modele,
        PROL_ZERO="OUI",
        ASSE=(
            _F(GROUP_MA=grmacalc, CHAM_GD=sqsursgrefe, NOM_CMP="X1", NOM_CMP_RESU="SIXX"),
            _F(GROUP_MA=grmacalc, CHAM_GD=sqsursgrefe, NOM_CMP="X2", NOM_CMP_RESU="SIYY"),
        ),
    )

    rdivaux = NonLinearResult()
    rdivaux.allocate(1)
    rdivaux.setField(sq2, "SIEF_ELGA", 0)
    rdivaux.setModel(modele, 0)

    formule = Formula()
    formule.setExpression("SIXX/SIYY*sigm_cnv")
    formule.setVariables(["SIXX", "SIYY"])
    formule.setContext({"sigm_cnv": dwb[grwb]["SIGM_CNV"]})

    rdiv1 = CALC_CHAMP(
        RESULTAT=rdivaux,
        GROUP_MA=grmacalc,
        CHAM_UTIL=_F(NOM_CHAM="SIEF_ELGA", FORMULE=formule, NUME_CHAM_RESU=1),
    )

    return CREA_CHAMP(
        OPERATION="ASSE",
        TYPE_CHAM="ELGA_NEUT_R",
        MODELE=modele,
        PROL_ZERO="OUI",
        ASSE=_F(
            GROUP_MA=grmacalc,
            CHAM_GD=rdiv1.getField("UT01_ELGA", 0),
            NOM_CMP="X1",
            NOM_CMP_RESU="X1",
        ),
    )


def sig1plasac(resultat, rsieq, numv1v2, dwb, reswbrest, grmapb, mclinst):
    """
    Major principal stress where the plasticity is active

    Arguments:
        resultat (NonLinearResult): Result to consider to compute Weibull
            stress
        rsieq (NonLinearResult): SIEQ_ELGA field
        numv1v2 (tuple): Indices of internal variables EPSPEQ and INDIPLAS
        dwb (dict): Weibull parameters
        reswbrest (NonLinearResult): Result where plasticity is strictly
            positive
        grmapb (str): Mesh cells group given in POST_BEREMIN
        mclinst (list): List of medcoupling time steps
            (iteration, order, time step) where there is plasticity

    Returns:
        FieldOnCells: ELGA_NEUT_R filled by PRIN_3
    """
    modele = resultat.getModel()
    if not grmapb[0] in dwb.keys():
        UTMESS("F", "RUPTURE1_88", valk=(grmapb[0]))

    maxsig = NonLinearResult()
    maxsig.allocate(rsieq.getNumberOfIndexes())

    indice = 0
    for nume_inst in rsieq.getAccessParameters()["NUME_ORDRE"]:
        grcalc = "mgrplas_{}".format(nume_inst)
        inst = rsieq.getTime(nume_inst)

        if inst in [elt[2] for elt in mclinst]:
            tronque = CALC_CHAMP(
                RESULTAT=resultat,
                INST=inst,
                GROUP_MA=grmapb[0],
                CHAM_UTIL=_F(
                    NOM_CHAM="VARI_ELGA",
                    FORMULE=FORMULE(
                        NOM_PARA=("V{}".format(numv1v2[0]), "V{}".format(numv1v2[1])),
                        VALE="indic_plasac(V{}, V{}, seuil)".format(numv1v2[0], numv1v2[1]),
                        indic_plasac=indic_plasac,
                        seuil=dwb[grmapb[0]]["SEUIL_EPSP_CUMU"],
                    ),
                    NUME_CHAM_RESU=1,
                ),
            )

            sign = CREA_CHAMP(
                OPERATION="ASSE",
                TYPE_CHAM="ELGA_SIEF_R",
                MODELE=modele,
                PROL_ZERO="OUI",
                ASSE=(
                    _F(
                        GROUP_MA=grcalc,
                        CHAM_GD=sigma1(rsieq, nume_inst, dwb, reswbrest, grmapb[0]),
                        NOM_CMP="X1",
                        NOM_CMP_RESU="SIXX",
                    ),
                    _F(
                        GROUP_MA=grcalc,
                        CHAM_GD=tronque.getField("UT01_ELGA", nume_inst),
                        NOM_CMP="X1",
                        NOM_CMP_RESU="SIYY",
                    ),
                ),
            )

            rsig1aux = NonLinearResult()
            rsig1aux.allocate(1)
            rsig1aux.setField(sign, "SIEF_ELGA", 0)
            rsig1aux.setModel(modele, 0)

            fmult = Formula()
            fmult.setExpression("SIXX*SIYY")
            fmult.setVariables(["SIXX", "SIYY"])
            fmult.setContext(dict())

            rsig1 = CALC_CHAMP(
                RESULTAT=rsig1aux,
                GROUP_MA=grcalc,
                CHAM_UTIL=_F(NOM_CHAM="SIEF_ELGA", FORMULE=fmult, NUME_CHAM_RESU=1),
            )

            sigtyp = CREA_CHAMP(
                OPERATION="ASSE",
                TYPE_CHAM="ELGA_SIEF_R",
                MODELE=modele,
                PROL_ZERO="OUI",
                ASSE=_F(
                    GROUP_MA=grcalc,
                    CHAM_GD=rsig1.getField("UT01_ELGA", 0),
                    NOM_CMP="X1",
                    NOM_CMP_RESU="SIXX",
                ),
            )

            maxsig.setField(sigtyp, "SIEF_ELGA", indice)
            maxsig.setTime(resultat.getTime(nume_inst), indice)
            indice = indice + 1

    return maxsig


def tps_maxsigm(rsieq, mclinst, maxsig, resanpb, bere_m):
    """
    Compute temporal maximum for each Gauss point and elevation to power m

    Arguments:
        rsieq (NonLinearResult): SIEQ_ELGA field
        mclinst (list): List of medcoupling time steps
            (iteration, order, time step) where there is plasticity
        maxsig (NonLinearResult): ELGA_NEUT_R filled by PRIN_3
        resanpb (NonLinearResult): Name of auxiliary result
        bere_m (float): Value of Beremin parameter M

    Returns:
        FieldOnCells:

            ELGA_NEUT_R filled by PRIN_3

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
    linstants = rsieq.getAccessParameters()["INST"]

    def puiss_m(valsixx):
        return valsixx ** bere_m

    indice = 0
    for nume_inst, inst in enumerate(linstants):

        if inst in [elt[2] for elt in mclinst]:

            if inst == maxsig.getTime(0):
                chmaxsig = maxsig.getField("SIEF_ELGA", 0)
            elif inst > maxsig.getTime(0):
                maxsig_r = NonLinearResult()
                maxsig_r.allocate(2)
                maxsig_r.setField(maxsig.getField("SIEF_ELGA", indice - 1), "SIEF_ELGA", 0)
                maxsig_r.setField(maxsig.getField("SIEF_ELGA", indice), "SIEF_ELGA", 1)

                chmaxsig = CREA_CHAMP(
                    TYPE_CHAM="ELGA_SIEF_R",
                    OPERATION="EXTR",
                    RESULTAT=maxsig_r,
                    NOM_CHAM="SIEF_ELGA",
                    TYPE_MAXI="MAXI",
                    TYPE_RESU="VALE",
                    NUME_ORDRE=(0, 1),
                )

            chsixxm = chmaxsig.transform(puiss_m)

            if resanpb is not None:
                resimpr.setField(chsixxm, "SIEF_ELGA", nume_inst)
                resimpr.setTime(inst, nume_inst)

            indice = indice + 1

        else:

            if resanpb is not None:

                chsixxm = FieldOnCellsReal(chmaxsig.getModel(), "ELGA", "SIEF_R")
                chsixxm.setValues(0)
                resimpr.setField(chsixxm, "SIEF_ELGA", nume_inst)
                resimpr.setTime(inst, nume_inst)

        sigw.setField(chsixxm, "SIEF_ELGA", nume_inst)
        sigw.setTime(inst, nume_inst)

    return (sigw, resimpr)


def compute_beremin_integral(model, coefmultpb, sigw, dwb, grmapb, resupb):
    """
    Integral computation for Weibull stress
    Produces also sigw**m and probability of failure

    Arguments:
        model (Model): Model used for resultat
        coefmultpb (float): Value of COEF_MULT in POST_BEREMIN
        sigw (NonLinearResult): ELGA_NEUT_R filled by PRIN_3
        dwb (dict): POST_BEREMIN keywords
        grmapb (str): Mesh cells group given in POST_BEREMIN
        resupb (NonLinearResult): Resultat input of POST_BEREMIN

    Returns:
        Table: POST_BEREMIN final table
    """
    poids = CALC_CHAM_ELEM(MODELE=model, OPTION="COOR_ELGA").getValuesWithDescription(
        "W", ["mgrplasfull"]
    )[0]
    sigwinst = sigw.getAccessParameters()["NUME_ORDRE"]
    valinte_sixx = []
    for nume_inst in sigwinst:
        nzval = np.array(
            (
                [
                    elt_poids * elt_sigw
                    for (elt_poids, elt_sigw) in zip(
                        poids,
                        sigw.getField("SIEF_ELGA", nume_inst).getValuesWithDescription(
                            "SIXX", ["mgrplasfull"]
                        )[0],
                    )
                ]
            )
        )
        valinte_sixx.append(np.sum(nzval))

    sigwaux = CREA_TABLE(
        LISTE=(_F(PARA="NUME_ORDRE", LISTE_I=sigwinst), _F(PARA="INTE_SIXX", LISTE_R=valinte_sixx))
    )

    if "SIGM_CNV" in dwb[grmapb[0]]:
        sigma_u = dwb[grmapb[0]]["SIGM_CNV"]
    else:
        sigma_u = dwb[grmapb[0]]["SIGM_REFE"]

    d_form = {}
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

    for clef in t_clef:
        d_form[clef] = Formula()
    for (clef, expression, variable, context) in zip(t_clef, t_expr, t_var, t_context):
        d_form[clef].setExpression(expression)
        d_form[clef].setVariables(variable)
        d_form[clef].setContext(context)

    sigwaux = CALC_TABLE(
        TABLE=sigwaux,
        reuse=sigwaux,
        ACTION=(
            _F(OPERATION="OPER", FORMULE=d_form["intfin"], NOM_PARA="SIGMA_WEIBULL"),
            _F(OPERATION="OPER", FORMULE=d_form["sigwm"], NOM_PARA="SIGMA_WEIBULL**M"),
            _F(OPERATION="OPER", FORMULE=d_form["probaw"], NOM_PARA="PROBA_WEIBULL"),
        ),
    )

    tab_aux = dict()
    for clef in ("SIGMA_WEIBULL", "SIGMA_WEIBULL**M", "PROBA_WEIBULL"):
        tab_aux[clef] = (sigwaux.EXTR_TABLE()).values()[clef]

    lentable = len(resupb.getAccessParameters()["INST"])
    liste_zero = [0] * (lentable - len(tab_aux["SIGMA_WEIBULL"]))
    tab_final = CREA_TABLE(
        LISTE=(
            _F(PARA="NUME_ORDRE", LISTE_I=resupb.getAccessParameters()["NUME_ORDRE"]),
            _F(PARA="INST", LISTE_R=resupb.getAccessParameters()["INST"]),
            _F(PARA="GROUP_MA", LISTE_K=[grmapb[0]] * lentable),
            _F(PARA="SIGMA_WEIBULL", LISTE_R=liste_zero + tab_aux["SIGMA_WEIBULL"]),
            _F(PARA="SIGMA_WEIBULL**M", LISTE_R=liste_zero + tab_aux["SIGMA_WEIBULL**M"]),
            _F(PARA="PROBA_WEIBULL", LISTE_R=liste_zero + tab_aux["PROBA_WEIBULL"]),
        )
    )

    return tab_final


def make_plasticity_groups(reswbrest, numv1v2, mclinst, seuil, l_epspmax):
    """
    Construct groups of plasticity for each time step when there is plasticity

    Arguments:
        reswbrest (NonLinearResult): result to consider to compute Weibull
            stress
        numv1v2 (tuple): Indices of internal variables EPSPEQ and INDIPLAS
        mclinst (list): List of medcoupling time steps
            (iteration, order, time step) where there is plasticity
        seuil (float): value of SEUIL_EPSP_CUMU
        l_epspmax (list): list of values of maximum of plasticity, if
            strictly greater than 0
    """
    mawbrest = reswbrest.getModel().getMesh()
    liter = [elt[0] for elt in mclinst]
    dval = {}
    for indice, iteration in enumerate(liter):

    for (indice, iteration) in enumerate(liter):

        if l_epspmax[indice] < seuil:
            dval["min{}".format(indice)] = 0
            dval["max{}".format(indice)] = l_epspmax[indice]
        else:
            dval["min{}".format(indice)] = seuil
            dval["max{}".format(indice)] = l_epspmax[indice]

    mawbrest = DEFI_GROUP(
        reuse=mawbrest,
        MAILLAGE=mawbrest,
        CREA_GROUP_NO=tuple(
            [
                _F(
                    NOM="vale_{}".format(iteration),
                    OPTION="INTERVALLE_VALE",
                    CHAM_GD=reswbrest.getField("VARI_ELGA", iteration).toFieldOnNodes(),
                    NOM_CMP="V{}".format(numv1v2[0]),
                    VALE=(dval["min{}".format(indice)], dval["max{}".format(indice)]),
                )
                for (indice, iteration) in enumerate(liter)
            ]
            + [
                _F(
                    NOM="ngrplas_{}".format(iteration),
                    INTERSEC=("vale_{}".format(iteration), "ngrmapb"),
                )
                for iteration in liter
            ]
        ),
    )

    mawbrest = DEFI_GROUP(
        reuse=mawbrest,
        MAILLAGE=mawbrest,
        CREA_GROUP_MA=tuple(
            _F(
                NOM="mgrplas_{}".format(iteration),
                OPTION="APPUI",
                TYPE_MAILLE="{}D".format(mawbrest.getDimension()),
                GROUP_NO="ngrplas_{}".format(iteration),
                TYPE_APPUI="AU_MOINS_UN",
            )
            for iteration in liter
        )
        + tuple(
            [
                _F(
                    NOM="mgrplasfull",
                    TYPE_MAILLE="{}D".format(mawbrest.getDimension()),
                    UNION=tuple("mgrplas_{}".format(iteration) for iteration in liter),
                )
            ]
        ),
    )
