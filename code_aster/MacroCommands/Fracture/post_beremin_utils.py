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
import tempfile
import medcoupling as mc

from libaster import EntityType

from ...Cata.Syntax import _F
from ...CodeCommands import DEFI_FICHIER, IMPR_RESU, CREA_CHAMP

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
                lwbmat = [
                    elvema.getValueReal("WEIBULL", "M"),
                    elvema.getValueReal("WEIBULL", "VOLU_REFE"),
                    elvema.getValueReal("WEIBULL", "SEUIL_EPSP_CUMU"),
                ]
                nomres = ["M", "VOLU_REFE", "SEUIL_EPSP_CUMU"]
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
                        dwb.update({grpname: dict(list(zip(nomres, lwbmat)))})
                elif etttype == EntityType.AllMeshEntitiesType:
                    dwb.update({group_ma: dict(list(zip(nomres, lwbmat)))})

    return dwb


def get_resu_from_deftype(resusd, grmapb, defopb):
    """
    Result to consider to compute Beremin stress

    Arguments:
        resusd (NonLinearResult): Resultat aster concept
        grmapb (str): Name of mesh cells group for Beremin computation from
                      command file
        defopb (str): Deformation type ("PETIT" or "GDEF_LOG")

    Returns:
        tuple:

            - First value is result to consider to compute Weibull stress
            - Second value is indices of internal variables EPSPEQ and INDIPLAS
            - Third value is list of medcoupling time steps, followed by
              (iteration, order, time step)
            - Fourth value is list of values of maximum of plasticity, if
              strictly greater than 0

    """
    rvga = make_rvga(resusd)

    med_filename_in = tempfile.NamedTemporaryFile(dir=".", suffix=".med").name

    unite = DEFI_FICHIER(FICHIER=med_filename_in, ACTION="ASSOCIER", TYPE="LIBRE", ACCES="NEW")

    IMPR_RESU(UNITE=unite, PROC0="NON", RESU=_F(RESULTAT=rvga, NOM_CHAM="VARI_ELGA"))

    DEFI_FICHIER(ACTION="LIBERER", UNITE=unite)

    mclinst = []
    l_epspmax = []
    mcchamp = mc.MEDFileData(med_filename_in).getFields()["{}VARI_ELGA_NOMME".format(rvga.userName)]
    for mctime in mcchamp.getTimeSteps():
        epspmax = (
            mcchamp[(mctime[0], mctime[1])]
            .getFieldAtLevel(mc.ON_GAUSS_PT, 0)
            .getArray()[:, 0]
            .getMaxValue()[0]
        )
        if epspmax > 0:
            mclinst.append(mctime)
            l_epspmax.append(epspmax)

    if len(mclinst) > 0:
        for field in mc.MEDFileFields(med_filename_in, False).getFieldsNames():
            if field[8:] == "VARI_ELGA_NOMME":
                numv1v2 = [
                    1
                    + [
                        elt[0] for elt in mc.GetComponentsNamesOfField(med_filename_in, field)
                    ].index(comp)
                    for comp in ("EPSPEQ", "INDIPLAS")
                ]

        if defopb == "GDEF_LOG":
            cmd_result = (
                get_resu_gdef_log(resusd, med_filename_in, rvga, grmapb),
                numv1v2,
                mclinst,
                l_epspmax,
            )

        else:
            cmd_result = (resusd, numv1v2, mclinst, l_epspmax)
    else:
        cmd_result = (None, None, mclinst, l_epspmax)

    return cmd_result


def make_rvga(resusd):
    """
    Result to consider to compute Beremin stress

    Arguments:
        resusd (NonLinearResult): Resultat aster concept

    Returns:
        NonLinearResult: Restriction of resusd to VARI_ELGA
    """
    rvga = NonLinearResult()
    rvga.allocate(resusd.getNumberOfIndexes())
    rvga.userName = "rvga____"

    for nume_inst in resusd.getIndexes():
        rvga.setField(resusd.getField("VARI_ELGA", nume_inst), "VARI_ELGA", nume_inst)
        rvga.setMaterialField(resusd.getMaterialField(nume_inst), nume_inst)
        rvga.setField(resusd.getField("COMPORTEMENT", nume_inst), "COMPORTEMENT", nume_inst)
        rvga.setTime(resusd.getTime(nume_inst), nume_inst)

    return rvga


def get_resu_gdef_log(resusd, med_filename_in, rvga, grmapb):
    """
    Result to consider to compute Beremin stress when DEFORMATION="GDEF_LOG"

    Arguments:
        resusd (NonLinearResult): Resultat aster concept
        med_filename_in (str): Filename of temporary file
        rvga (NonLinearResult): Result of VARI_ELGA
        grmapb (str): Name of mesh cells group for Beremin computation from
                      command file

    Returns:
        NonLinearResult: Result to consider to compute Weibull stress
    """
    dim_geom = resusd.getModel().getMesh().getDimension()

    for field in mc.MEDFileFields(med_filename_in, False).getFieldsNames():
        if field[8:] == "VARI_ELGA_NOMME":
            numcmpsel = [
                1
                + [elt[0] for elt in mc.GetComponentsNamesOfField(med_filename_in, field)].index(
                    comp
                )
                for comp in {
                    2: ("TXX", "TYY", "TZZ", "TXY"),
                    3: ("TXX", "TYY", "TZZ", "TXY", "TXZ", "TYZ"),
                }[dim_geom]
            ]

    reswbrest = NonLinearResult()
    reswbrest.allocate(resusd.getNumberOfIndexes())

    for rank in resusd.getIndexes():
        chvga = rvga.getField("VARI_ELGA", rank)

        sglog = CREA_CHAMP(
            OPERATION="ASSE",
            TYPE_CHAM="ELGA_SIEF_R",
            MODELE=resusd.getModel(rank),
            PROL_ZERO="OUI",
            ASSE=(
                tuple(
                    _F(CHAM_GD=chvga, GROUP_MA=grmapb, NOM_CMP=nom_cmp, NOM_CMP_RESU=nom_cmp_resu)
                    for (nom_cmp, nom_cmp_resu) in zip(
                        tuple("V{}".format(ncs) for ncs in numcmpsel),
                        {
                            2: ("SIXX", "SIYY", "SIZZ", "SIXY"),
                            3: ("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
                        }[dim_geom],
                    )
                )
            ),
        )

        reswbrest.setField(resusd.getField("COMPORTEMENT", rank), "COMPORTEMENT", rank)
        reswbrest.setField(sglog, "SIEF_ELGA", rank)
        reswbrest.setField(chvga, "VARI_ELGA", rank)
        reswbrest.setModel(resusd.getModel(rank), rank)
        reswbrest.setTime(resusd.getTime(rank), rank)

    return reswbrest
