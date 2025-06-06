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

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *
from .affe_cara_elem import AFFE_CARA_ELEM
from .affe_char_meca import AFFE_CHAR_MECA
from .affe_materiau import AFFE_MATERIAU
from .affe_modele import AFFE_MODELE


from ..Commons.c_comportement import compat_syntax


def pre_seisme_nonl_sdprod(self, RESULTAT, **args):
    if args.get("__all__"):
        return (
            [None],
            [None, evol_noli],
            [None, modele_sdaster],
            [None, maillage_sdaster],
            [None, cham_mater],
            [None, cara_elem],
            [None, mode_meca],
            [None, macr_elem_dyna],
            [None, char_meca],
        )

    if RESULTAT[0]["RESULTAT"]:
        self.type_sdprod(RESULTAT[0]["RESULTAT"], evol_noli)
    if RESULTAT[0]["MODELE"]:
        self.type_sdprod(RESULTAT[0]["MODELE"], modele_sdaster)
    if RESULTAT[0]["MAILLAGE"]:
        self.type_sdprod(RESULTAT[0]["MAILLAGE"], maillage_sdaster)
    if RESULTAT[0]["CHAM_MATER"]:
        self.type_sdprod(RESULTAT[0]["CHAM_MATER"], cham_mater)
    if RESULTAT[0]["CARA_ELEM"]:
        self.type_sdprod(RESULTAT[0]["CARA_ELEM"], cara_elem)
    if RESULTAT[0]["BASE_MODALE"]:
        if args["PRE_CALC_MISS"]:
            self.type_sdprod(RESULTAT[0]["BASE_MODALE"], mode_meca)
        else:
            raise CataError(
                "Le mot-clé PRE_CALC_MISS est obligatoire pour créer un concept de type BASE_MODALE"
            )
    if RESULTAT[0]["MACR_ELEM_DYNA"]:
        if args["PRE_CALC_MISS"]:
            self.type_sdprod(RESULTAT[0]["MACR_ELEM_DYNA"], macr_elem_dyna)
        else:
            raise CataError(
                "Le mot-clé PRE_CALC_MISS est obligatoire pour créer un concept de type MACR_ELEM_DYNA"
            )
    if RESULTAT[0]["CHARGE"]:
        for mcfact in RESULTAT[0]["CHARGE"]:
            if mcfact["OPTION"] == "LAPL_TEMPS" and args["PRE_CALC_MISS"]:
                raise CataError(
                    "Le mot-clé POST_CALC_MISS est obligatoire pour créer une charge de type LAPL_TEMPS"
                )
            self.type_sdprod(mcfact["NOM"], char_meca)


from ..Commons.c_comportement import compat_syntax


def affe_char_meca_keywords():
    """Keywords for AFFE_CHAR_MECA with small changes"""
    # DO NOT CHANGE ORIGINAL KEYWORDS!!!
    orig = AFFE_CHAR_MECA.entites.copy()
    orig["MODELE"] = SIMP(statut="f", typ=modele_sdaster)
    return orig


from ..Commons.c_comportement import compat_syntax


def affe_cara_elem_keywords():
    """Keywords for AFFE_CARA_ELEM with small changes"""
    # DO NOT CHANGE ORIGINAL KEYWORDS!!!
    orig = AFFE_CARA_ELEM.entites.copy()
    orig["MODELE"] = SIMP(statut="f", typ=modele_sdaster)
    return orig


PRE_SEISME_NONL = MACRO(
    nom="PRE_SEISME_NONL",
    op=OPS("code_aster.MacroCommands.pre_seisme_nonl_ops.pre_seisme_nonl_ops"),
    compat_syntax=compat_syntax,
    sd_prod=pre_seisme_nonl_sdprod,
    fr=tr("description"),
    reentrant="n",
    AFFE_MODELE=FACT(statut="f", regles=AFFE_MODELE.regles, **AFFE_MODELE.entites),
    AFFE_MATERIAU=FACT(statut="f", regles=AFFE_MATERIAU.regles, **AFFE_MATERIAU.entites),
    AFFE_CARA_ELEM=FACT(statut="f", regles=AFFE_CARA_ELEM.regles, **affe_cara_elem_keywords()),
    AFFE_CHAR_MECA=FACT(statut="f", regles=AFFE_CHAR_MECA.regles, **affe_char_meca_keywords()),
    PRE_CALC_MISS=FACT(
        statut="f",
        max=1,
        REDUC_DYNA_ISS=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
        REDUC_DYNA_IFS=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
        NMAX_MODE_ISS=SIMP(statut="o", typ="I"),
        b_ISFS=BLOC(
            condition=""" equal_to("CALC_MISS_OPTION", 'ISFS') """,
            NMAX_MODE_IFS=SIMP(statut="o", typ="I"),
        ),
        GROUP_NO_CENT=SIMP(statut="f", typ=grno, max="**"),
        CALC_MISS_OPTION=SIMP(statut="o", typ="TXM", into=("ISS", "ISFS")),
        GROUP_MA_INTERF=SIMP(statut="o", typ=grma, max="**"),
    ),
    POST_CALC_MISS=FACT(
        statut="f",
        max=1,
        MACR_ELEM_DYNA=SIMP(
            statut="o", typ=macr_elem_dyna, fr=tr("Macro élément produit en amont")
        ),
        GROUP_NO_CENT=SIMP(statut="o", typ=grno, max="**"),
        GROUP_MA_INTERF=SIMP(statut="o", typ=grma, max="**"),
        #  b_unite_resu = BLOC(regles= AU_MOINS_UN('UNITE_RESU_RIGI','UNITE_RESU_AMOR','UNITE_RESU_MASS'),
        UNITE_RESU_RIGI=SIMP(statut="f", typ=UnitType(), inout="in"),
        UNITE_RESU_MASS=SIMP(statut="f", typ=UnitType(), inout="in"),
        UNITE_RESU_AMOR=SIMP(statut="f", typ=UnitType(), inout="in"),
        #                       ),
    ),
    STAT_DYNA=FACT(
        statut="f",
        min=1,
        max=1,
        COMPORTEMENT=C_COMPORTEMENT("MECA_NON_LINE"),
        CONVERGENCE=C_CONVERGENCE("MECA_NON_LINE"),
        RESULTAT=SIMP(statut="o", typ=evol_noli),
        EXCIT=FACT(
            statut="o",
            max="**",
            CHARGE=SIMP(statut="o", typ=(char_meca, char_cine_meca)),
            FONC_MULT=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
            TYPE_CHARGE=SIMP(
                statut="f",
                typ="TXM",
                defaut="FIXE_CSTE",
                into=("FIXE_CSTE", "FIXE_PILO", "SUIV", "DIDI"),
            ),
        ),
        BASE_MODALE=SIMP(statut="o", typ=mode_meca),
        UNITE_IMPE_TEMPS=FACT(
            statut="o",
            max="**",
            #         b_unite_resu = BLOC(regles= AU_MOINS_UN('UNITE_RESU_RIGI','UNITE_RESU_AMOR','UNITE_RESU_MASS'),
            UNITE_RESU_RIGI=SIMP(statut="f", typ=UnitType(), inout="in"),
            UNITE_RESU_MASS=SIMP(statut="f", typ=UnitType(), inout="in"),
            UNITE_RESU_AMOR=SIMP(statut="f", typ=UnitType(), inout="in"),
            #                    ),
        ),
        UNITE_IMPE_FREQ=SIMP(statut="o", typ=UnitType(), inout="in"),
        FORCE_SOL=SIMP(statut="o", typ=char_meca),
        COEF_AMOR=SIMP(statut="f", typ="R", defaut=0.0),
        NB_INST=SIMP(statut="f", typ="R", defaut=100.0),
    ),
    RESULTAT=FACT(
        statut="o",
        min=1,
        max=1,
        RESULTAT=SIMP(statut="f", typ=CO),
        MODELE=SIMP(statut="f", typ=CO),
        MAILLAGE=SIMP(statut="f", typ=CO),
        CHAM_MATER=SIMP(statut="f", typ=CO),
        CARA_ELEM=SIMP(statut="f", typ=CO),
        # ONLY when PRE_CALC_MISS is not None
        BASE_MODALE=SIMP(statut="f", typ=CO),
        # ONLY when PRE_CALC_MISS is not None
        MACR_ELEM_DYNA=SIMP(statut="f", typ=CO),
        CHARGE=FACT(
            statut="f",
            max="**",
            OPTION=SIMP(statut="o", typ="TXM", into=("COND_LIM", "LAPL_TEMPS")),
            NOM=SIMP(statut="o", typ=CO),
        ),
    ),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
