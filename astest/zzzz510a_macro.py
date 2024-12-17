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
POST_BEREMIN option SIGM_ELMOY codée par Sébastien Meunier
"""

from code_aster.Commands import *
from code_aster import CA

test = CA.TestCase()


def sigelmoy(
    reswbrest,
    chfmu,
    grmapb,
    # mclinst
):
    """
    Compute SIGM_ELMOY
    """
    dim_geom = reswbrest.getModel().getMesh().getDimension()

    rmelmoy = CA.NonLinearResult()
    rmelmoy.allocate(reswbrest.getNumberOfIndexes())

    poids = CALC_CHAM_ELEM(MODELE=reswbrest.getModel(), GROUP_MA=grmapb, OPTION="COOR_ELGA")

    pds, (cells, _, _, _) = poids.getValuesWithDescription("W", [grmapb])
    volume = dict()
    for maille in list(set(cells)):
        volume[maille] = 0

    for (num_maille, maille) in enumerate(cells):
        volume.update({maille: volume[maille] + pds[num_maille]})

    sief_elga = CA.NonLinearResult()
    # sief_elga.allocate(len(mclinst))
    sief_elga.allocate(reswbrest.getNumberOfIndexes())

    vari_elga = CA.NonLinearResult()
    # vari_elga.allocate(len(mclinst))
    vari_elga.allocate(reswbrest.getNumberOfIndexes())

    for iteration in reswbrest.getIndexes():

        sief_elga.setField(reswbrest.getField("SIEF_ELGA", iteration), "SIEF_ELGA", iteration)
        sief_elga.setTime(reswbrest.getTime(iteration), iteration)
        sief_elga.setModel(reswbrest.getModel(iteration), iteration)

        vari_elga.setField(reswbrest.getField("VARI_ELGA", iteration), "VARI_ELGA", iteration)
        vari_elga.setTime(reswbrest.getTime(iteration), iteration)
        vari_elga.setModel(reswbrest.getModel(iteration), iteration)

    rmelmoy2 = CALC_CHAMP(RESULTAT=reswbrest, GROUP_MA=grmapb, CONTRAINTE="SIMY_ELGA")
    import numpy as np

    for (nume_inst, inst) in enumerate(sief_elga.getAccessParameters()["INST"]):

        selga = CREA_CHAMP(
            TYPE_CHAM="ELGA_SIEF_R",
            OPERATION="EXTR",
            RESULTAT=sief_elga,
            NOM_CHAM="SIEF_ELGA",
            INST=inst,
        )

        velga = CREA_CHAMP(
            TYPE_CHAM="ELGA_VARI_R",
            OPERATION="EXTR",
            RESULTAT=vari_elga,
            NOM_CHAM="VARI_ELGA",
            INST=inst,
        )

        chbid = CREA_CHAMP(
            OPERATION="AFFE",
            TYPE_CHAM="ELGA_SIEF_R",
            MODELE=reswbrest.getModel(),
            PROL_ZERO="OUI",
            AFFE=_F(GROUP_MA=grmapb, NOM_CMP="SIXX", VALE=0),
        )

        sixx = CREA_CHAMP(
            OPERATION="ASSE",
            TYPE_CHAM="ELGA_NEUT_R",
            MODELE=sief_elga.getModel(),
            PROL_ZERO="OUI",
            ASSE=(
                _F(GROUP_MA=grmapb, CHAM_GD=poids, NOM_CMP="W", NOM_CMP_RESU="X1"),
                _F(GROUP_MA=grmapb, CHAM_GD=selga, NOM_CMP="SIXX", NOM_CMP_RESU="X2"),
            ),
        )

        siyy = CREA_CHAMP(
            OPERATION="ASSE",
            TYPE_CHAM="ELGA_NEUT_R",
            MODELE=sief_elga.getModel(),
            PROL_ZERO="OUI",
            ASSE=(
                _F(GROUP_MA=grmapb, CHAM_GD=poids, NOM_CMP="W", NOM_CMP_RESU="X1"),
                _F(GROUP_MA=grmapb, CHAM_GD=selga, NOM_CMP="SIYY", NOM_CMP_RESU="X2"),
            ),
        )

        sizz = CREA_CHAMP(
            OPERATION="ASSE",
            TYPE_CHAM="ELGA_NEUT_R",
            MODELE=sief_elga.getModel(),
            PROL_ZERO="OUI",
            ASSE=(
                _F(GROUP_MA=grmapb, CHAM_GD=poids, NOM_CMP="W", NOM_CMP_RESU="X1"),
                _F(GROUP_MA=grmapb, CHAM_GD=selga, NOM_CMP="SIZZ", NOM_CMP_RESU="X2"),
            ),
        )

        sixy = CREA_CHAMP(
            OPERATION="ASSE",
            TYPE_CHAM="ELGA_NEUT_R",
            MODELE=sief_elga.getModel(),
            PROL_ZERO="OUI",
            ASSE=(
                _F(GROUP_MA=grmapb, CHAM_GD=poids, NOM_CMP="W", NOM_CMP_RESU="X1"),
                _F(GROUP_MA=grmapb, CHAM_GD=selga, NOM_CMP="SIXY", NOM_CMP_RESU="X2"),
            ),
        )

        dsmoy = dict()
        for (_cmp, cham_para) in zip(("xx", "yy", "zz", "xy"), (sixx, siyy, sizz, sixy)):

            dsmoy.update(
                {
                    _cmp: CREA_CHAMP(
                        OPERATION="EVAL", TYPE_CHAM="ELGA_NEUT_R", CHAM_F=chfmu, CHAM_PARA=cham_para
                    )
                }
            )

        if dim_geom == 3:
            sixz = CREA_CHAMP(
                OPERATION="ASSE",
                TYPE_CHAM="ELGA_NEUT_R",
                MODELE=sief_elga.getModel(),
                PROL_ZERO="OUI",
                ASSE=(
                    _F(GROUP_MA=grmapb, CHAM_GD=poids, NOM_CMP="W", NOM_CMP_RESU="X1"),
                    _F(GROUP_MA=grmapb, CHAM_GD=selga, NOM_CMP="SIXZ", NOM_CMP_RESU="X2"),
                ),
            )

            siyz = CREA_CHAMP(
                OPERATION="ASSE",
                TYPE_CHAM="ELGA_NEUT_R",
                MODELE=sief_elga.getModel(),
                PROL_ZERO="OUI",
                ASSE=(
                    _F(GROUP_MA=grmapb, CHAM_GD=poids, NOM_CMP="W", NOM_CMP_RESU="X1"),
                    _F(GROUP_MA=grmapb, CHAM_GD=selga, NOM_CMP="SIYZ", NOM_CMP_RESU="X2"),
                ),
            )

            for (_cmp, cham_para) in zip(("xz", "yz"), (sixz, siyz)):

                dsmoy.update(
                    {
                        _cmp: CREA_CHAMP(
                            OPERATION="EVAL",
                            TYPE_CHAM="ELGA_NEUT_R",
                            CHAM_F=chfmu,
                            CHAM_PARA=cham_para,
                        )
                    }
                )

        sxxmoyvale, (cxx, _, _, _) = dsmoy["xx"].getValuesWithDescription("X1", [grmapb])
        syymoyvale, (cyy, _, _, _) = dsmoy["yy"].getValuesWithDescription("X1", [grmapb])
        szzmoyvale, (czz, _, _, _) = dsmoy["zz"].getValuesWithDescription("X1", [grmapb])
        sxymoyvale, (cxy, _, _, _) = dsmoy["xy"].getValuesWithDescription("X1", [grmapb])

        if dim_geom == 3:
            sxzmoyvale, (cxz, _, _, _) = dsmoy["xz"].getValuesWithDescription("X1", [grmapb])
            syzmoyvale, (cyz, _, _, _) = dsmoy["yz"].getValuesWithDescription("X1", [grmapb])

        sxxmoyfinal = dict()
        for maille in list(set(cxx)):
            sxxmoyfinal[maille] = 0

        for (num_maille, maille) in enumerate(cxx):
            sxxmoyfinal.update(
                {maille: sxxmoyfinal[maille] + 1 / volume[maille] * sxxmoyvale[num_maille]}
            )

        syymoyfinal = dict()
        for maille in list(set(cyy)):
            syymoyfinal[maille] = 0

        for (num_maille, maille) in enumerate(cyy):
            syymoyfinal.update(
                {maille: syymoyfinal[maille] + 1 / volume[maille] * syymoyvale[num_maille]}
            )

        szzmoyfinal = dict()
        for maille in list(set(czz)):
            szzmoyfinal[maille] = 0

        for (num_maille, maille) in enumerate(czz):
            szzmoyfinal.update(
                {maille: szzmoyfinal[maille] + 1 / volume[maille] * szzmoyvale[num_maille]}
            )

        sxymoyfinal = dict()
        for maille in list(set(cxy)):
            sxymoyfinal[maille] = 0

        for (num_maille, maille) in enumerate(cxy):
            sxymoyfinal.update(
                {maille: sxymoyfinal[maille] + 1 / volume[maille] * sxymoyvale[num_maille]}
            )

        if dim_geom == 3:
            sxzmoyfinal = dict()
            for maille in list(set(cxz)):
                sxzmoyfinal[maille] = 0

            for (num_maille, maille) in enumerate(cxz):
                sxzmoyfinal.update(
                    {maille: sxzmoyfinal[maille] + 1 / volume[maille] * sxzmoyvale[num_maille]}
                )

            syzmoyfinal = dict()
            for maille in list(set(cyz)):
                syzmoyfinal[maille] = 0

            for (num_maille, maille) in enumerate(cyz):
                syzmoyfinal.update(
                    {maille: syzmoyfinal[maille] + 1 / volume[maille] * syzmoyvale[num_maille]}
                )

        l_mailles_grma = chbid.getValuesWithDescription("SIXX", [grmapb])[1][0]
        l_mailles_total = chbid.getValuesWithDescription("SIXX", [])[1][0]
        nbrvale = len(chbid.getValues())
        champ = [0] * len(l_mailles_total) * 2 * dim_geom
        assert nbrvale == len(champ), f"{nbrvale}, {len(champ)}"

        for (indice, maille) in enumerate(l_mailles_total):
            if dim_geom == 2 and maille in l_mailles_grma:
                champ[4 * indice] = sxxmoyfinal[maille]
                champ[4 * indice + 1] = syymoyfinal[maille]
                champ[4 * indice + 2] = szzmoyfinal[maille]
                champ[4 * indice + 3] = sxymoyfinal[maille]
            elif dim_geom == 3 and maille in l_mailles_grma:
                champ[6 * indice] = sxxmoyfinal[maille]
                champ[6 * indice + 1] = syymoyfinal[maille]
                champ[6 * indice + 2] = szzmoyfinal[maille]
                champ[6 * indice + 3] = sxymoyfinal[maille]
                champ[6 * indice + 4] = sxzmoyfinal[maille]
                champ[6 * indice + 5] = syzmoyfinal[maille]
        chbid.setValues(champ)

        rmelmoy.setField(chbid, "SIEF_ELGA", nume_inst)
        rmelmoy.setField(velga, "VARI_ELGA", nume_inst)
        rmelmoy.setTime(inst, nume_inst)
        rmelmoy.setModel(reswbrest.getModel(nume_inst), nume_inst)
        rmelmoy.setField(reswbrest.getField("COMPORTEMENT", nume_inst), "COMPORTEMENT", nume_inst)
        if True:
            field = rmelmoy2.getField("SIMY_ELGA", nume_inst)
            tab1 = np.array(chbid.getValues())
            tab2 = np.array(field.getValues())
            if np.allclose(tab1, tab2, atol=1e-12) is not True:
                print(chbid.getValues())
                print(field.getValues())
                test.assertTrue(False)
            else:
                test.assertTrue(True)

    return rmelmoy
