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

"""
Calcul de proprietés homo
"""
from libaster import Physics, Modelings
from ...Cata.Syntax import _F
from ...CodeCommands import LIRE_MAILLAGE, LIRE_RESU, CREA_TABLE
from ...Objects import Model, ThermalResultDict, ElasticResultDict
from ...Messages import ASSERT, UTMESS
from ...Helpers.LogicalUnit import LogicalUnitFile
from ...Utilities import medcoupling as medc
from . import NameConverter


def parse_med(unit):
    """
    Check if MED file contains the homogeneisation corrector field.

    Arguments
    ---------
        unit (int) : Logical unit associated to the MED file.

    Returns
    -------
        meca (dict) : Name associaton between med elastic fields and CALC_MATE_HOMO
        ther (dict) : Name associaton between med thermal fields and CALC_MATE_HOMO

    """
    fname = LogicalUnitFile.filename_from_unit(unit)
    mdata = medc.MEDFileData(fname)
    mednames = [i.getName() for i in mdata.getFields()]
    ASSERT(all(n.startswith("Z") for n in mednames))

    meca = {n: NameConverter.fromMed(n) for n in mednames if n.endswith("DEPL")}
    ther = {n: NameConverter.fromMed(n) for n in mednames if n.endswith("TEMP")}

    astnames = [NameConverter.fromMed(n) for n in mednames]
    tabpara = CREA_TABLE(
        LISTE=(
            _F(PARA="NOM_MED", TYPE_K="K24", LISTE_K=mednames),
            _F(PARA="NOM_AST", TYPE_K="K24", LISTE_K=astnames),
        )
    )

    return meca, ther, tabpara


def load_elas_fields(fnames, mesh, unit):
    """
    Load elastic corrector fields from MED file.

    Arguments
    ---------
        fnames (dict) : Name associaton between med fields and CALC_MATE_HOMO
        mesh (Mesh) : Support mesh.
        unit (int) : Logical unit associated to the MED file.

    Returns
    -------
        meca (ElasticResultDict) : The elastic corrector fields.
    """

    ASSERT(len(fnames) > 0)
    mod_meca = Model(mesh)
    mod_meca.addModelingOnMesh(Physics.Mechanics, Modelings.Tridimensional)
    mod_meca.build()
    meca = ElasticResultDict()
    for medname, corrname in fnames.items():
        meca[corrname] = LIRE_RESU(
            TYPE_RESU="EVOL_ELAS",
            FORMAT="MED",
            MODELE=mod_meca,
            FORMAT_MED=(_F(NOM_CHAM_MED=medname, NOM_CHAM="DEPL"),),
            UNITE=unit,
            TOUT_ORDRE="OUI",
        )

    return meca


def load_ther_fields(fnames, mesh, unit):
    """
    Load thermal corrector fields from MED file.

    Arguments
    ---------
        fnames (dict) : Name associaton between med fields and CALC_MATE_HOMO
        mesh (Mesh) : Support mesh.
        unit (int) : Logical unit associated to the MED file.

    Returns
    -------
       ther (ThermalResultDict) : The thermal corrector fields.
    """
    ASSERT(len(fnames) > 0)
    mod_ther = Model(mesh)
    mod_ther.addModelingOnMesh(Physics.Thermal, Modelings.Tridimensional)
    mod_ther.build()
    ther = ThermalResultDict()
    for medname, corrname in fnames.items():
        ther[corrname] = LIRE_RESU(
            TYPE_RESU="EVOL_THER",
            FORMAT="MED",
            MODELE=mod_ther,
            FORMAT_MED=(_F(NOM_CHAM_MED=medname, NOM_CHAM="TEMP"),),
            UNITE=unit,
            TOUT_ORDRE="OUI",
        )

    return ther


def load_fields(unit):
    """Read corrector fields from med file

    Arguments
    ---------
        unit (int) : Logical unit associated to the MED file.

    Returns
    -------
       tabpara (Table) : List of MED fields found.
       meca (ElasticResultDict) : The elastic corrector fields.
       ther (ThermalResultDict) : The thermal corrector fields or None.
    """

    fnames_depl, fnames_ther, tabpara = parse_med(unit)
    mesh = LIRE_MAILLAGE(FORMAT="MED", UNITE=unit)

    corr_meca = None
    if fnames_depl:
        corr_meca = load_elas_fields(fnames_depl, mesh, unit)

    corr_ther = None
    if fnames_ther:
        corr_ther = load_ther_fields(fnames_ther, mesh, unit)

    return tabpara, corr_meca, corr_ther


def lire_corr_ops(self, **kwargs):
    """
    Main function for reading correctors from MED.

    """

    unit = kwargs.get("UNITE")
    tabpara, corr_meca, corr_ther = load_fields(unit)

    if kwargs.get("CORR_MECA") and corr_meca:
        self.register_result(corr_meca, kwargs.get("CORR_MECA"))

    if kwargs.get("CORR_THER") and corr_ther:
        self.register_result(corr_ther, kwargs.get("CORR_THER"))

    return tabpara
