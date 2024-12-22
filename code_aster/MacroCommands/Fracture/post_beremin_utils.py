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
Utilitary functions for POST_BEREMIN
"""
from libaster import EntityType

from ...Objects import NonLinearResult


def get_beremin_properties(resusd, group_ma):
    """
    For each material field, returns eventual Beremin parameters

    Arguments:
        resusd (NonLinearResult): Resultat aster concept
        group_ma (list): List of names of mesh cells group

    Returns:
        dict: Keys are names of mesh cells groups, values are Beremin
        parameters
    """
    dwb = dict()
    for eltmater in resusd.getMaterialField().getVectorOfPartOfMaterialField():
        etttype = eltmater.getMeshEntity().getType()
        for elvema in eltmater.getVectorOfMaterial():
            if "WEIBULL" in elvema.getMaterialNames():
                nomres = ["M", "VOLU_REFE", "SEUIL_EPSP_CUMU", "SIGM_SEUIL"]
                lwbmat = [elvema.getValueReal("WEIBULL", param) for param in nomres]
                try:
                    lwbmat += [
                        elvema.getValueReal("WEIBULL", "SIGM_CNV"),
                        elvema.getFunction("WEIBULL", "SIGM_REFE"),
                    ]
                    nomres += ["SIGM_CNV", "SIGM_REFE"]
                except RuntimeError:
                    lwbmat += [elvema.getValueReal("WEIBULL", "SIGM_REFE")]
                    nomres += ["SIGM_REFE"]
                if etttype == EntityType.GroupOfCellsType:
                    for grpname in eltmater.getMeshEntity().getNames():
                        dwb.update({grpname: dict(zip(nomres, lwbmat))})
                elif etttype == EntityType.AllMeshEntitiesType:
                    dwb.update({group_ma: dict(zip(nomres, lwbmat))})

    return dwb


def get_resu_from_deftype(resusd, grmapb, defopb, numvs):
    """
    Result to consider to compute Beremin stress

    Arguments:
        resusd (NonLinearResult): Resultat aster input
        grmapb (str): Mesh cells group for Beremin (one group only).
        defopb (str): {"PETIT", "GDEF_LOG"} Deformation type
        numvs (dict): {"vari":LIST_NUME_VARI, "sief": LIST_NUME_SIEF}

    Returns:
        tuple:

            - First value is result to consider to compute Weibull stress
            - Second value is list of time steps (nume_ordre, time step), in plasticity only
            - Third value is list of values of maximum of plasticity, if
              strictly greater than 0

    """
    rvga = NonLinearResult()
    rvga.allocate(resusd.getNumberOfIndexes())
    rvga.userName = "rvga____"

    l_instplas = []
    l_epspmax = []
    for nume_inst in resusd.getIndexes():
        chvari = resusd.getField("VARI_ELGA", nume_inst)
        maxepspeq = max(
            chvari.getValuesWithDescription("V{}".format(numvs["vari"][0]), [grmapb])[0]
        )
        if maxepspeq > 0:
            l_instplas.append((nume_inst, resusd.getTime(nume_inst)))
            l_epspmax.append(maxepspeq)
        rvga.setField(chvari, "VARI_ELGA", nume_inst)
        rvga.setMaterialField(resusd.getMaterialField(nume_inst), nume_inst)
        rvga.setField(resusd.getField("COMPORTEMENT", nume_inst), "COMPORTEMENT", nume_inst)
        rvga.setTime(resusd.getTime(nume_inst), nume_inst)

    if len(l_instplas) > 0:

        if defopb == "GDEF_LOG":

            cmd_result = (get_resu_gdef_log(resusd, rvga, numvs["sief"]), l_instplas, l_epspmax)

        else:
            cmd_result = (resusd, l_instplas, l_epspmax)
    else:
        cmd_result = (None, l_instplas, l_epspmax)

    return cmd_result


def get_resu_gdef_log(resusd, rvga, numsief):
    """
    Result to consider to compute Beremin stress when DEFORMATION="GDEF_LOG"

    Arguments:
        resusd (NonLinearResult): Resultat aster concept
        rvga (NonLinearResult): Result of VARI_ELGA
        numsief (list): Indices of logarithmic stresses

    Returns:
        NonLinearResult: Result to consider to compute Weibull stress
    """
    dim_geom = resusd.getModel().getMesh().getDimension()

    reswbrest = NonLinearResult()
    reswbrest.allocate(resusd.getNumberOfIndexes())
    for rank in resusd.getIndexes():
        chvga = rvga.getField("VARI_ELGA", rank)

        reswbrest.setField(resusd.getField("COMPORTEMENT", rank), "COMPORTEMENT", rank)
        reswbrest.setField(
            chvga.asPhysicalQuantity(
                "SIEF_R",
                dict(
                    zip(
                        tuple(f"V{ncs}" for ncs in numsief),
                        {
                            2: ("SIXX", "SIYY", "SIZZ", "SIXY"),
                            3: ("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
                        }[dim_geom],
                    )
                ),
            ),
            "SIEF_ELGA",
            rank,
        )
        reswbrest.setField(chvga, "VARI_ELGA", rank)
        reswbrest.setModel(resusd.getModel(rank), rank)
        reswbrest.setTime(resusd.getTime(rank), rank)

    return reswbrest
