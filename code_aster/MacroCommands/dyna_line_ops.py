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
# person_in_charge: yannick.tampango at edf.fr
#
from functools import reduce

import numpy as np

import aster
from libaster import AsterError
from ..Cata.Syntax import _F
from ..CodeCommands import (
    AFFE_CHAR_MECA,
    ASSE_MATRICE,
    ASSE_VECTEUR,
    CALC_CHAR_SEISME,
    CALC_FONCTION,
    CALC_MATR_ELEM,
    CALC_MISS,
    CALC_MODES,
    CALC_VECT_ELEM,
    CO,
    COMB_MATR_ASSE,
    CREA_CHAMP,
    CREA_RESU,
    DEFI_BASE_MODALE,
    DEFI_INTERF_DYNA,
    DYNA_VIBRA,
    FACTORISER,
    MACR_ELEM_DYNA,
    MACRO_ELAS_MULT,
    MACRO_MATR_AJOU,
    MODE_STATIQUE,
    NUME_DDL,
    NUME_DDL_GENE,
    PROJ_MATR_BASE,
    PROJ_VECT_BASE,
    RESOUDRE,
    REST_GENE_PHYS,
)
from ..Utilities import force_list
from ..Messages import UTMESS


class DynaLineFEM:
    """hold the FEM model used for dyna_line"""

    def __init__(
        self,
        parent,
        MODELE,
        TYPE_CALCUL,
        BASE_CALCUL,
        CHARGE=None,
        AMORTISSEMENT=None,
        CHAM_MATER=None,
        CARA_ELEM=None,
        ISS="NON",
        **args
    ):
        self.parent = parent
        self.keywords = {"MODELE": MODELE}
        try:
            charge = tuple([x for x in CHARGE if x.getType()[:9] == "CHAR_MECA"])
        except TypeError:
            charge = None
        if charge:
            self.keywords["CHARGE"] = charge
        self.char_cine = {}
        try:
            char_cine = tuple([x for x in CHARGE if x.getType()[:14] == "CHAR_CINE_MECA"])
        except TypeError:
            char_cine = None
        if char_cine:
            self.char_cine["CHAR_CINE"] = char_cine
        if CHAM_MATER:
            self.keywords["CHAM_MATER"] = CHAM_MATER
        if CARA_ELEM:
            self.keywords["CARA_ELEM"] = CARA_ELEM
        self.__isTypeTran = TYPE_CALCUL == "TRAN"
        self.__isTypeHarm = TYPE_CALCUL == "HARM"
        assert self.__isTypeTran ^ self.__isTypeHarm
        self.__isBaseGene = BASE_CALCUL == "GENE"
        self.__isBasePhys = BASE_CALCUL == "PHYS"
        assert self.__isBaseGene ^ self.__isBasePhys
        if AMORTISSEMENT:
            self.amortissement = force_list(AMORTISSEMENT)[0]
        else:
            self.amortissement = None
        self.iss = ISS == "OUI"
        if self.iss:
            self.__init_iss(**args)

    def __init_iss(self, GROUP_NO_INTERF=None, **args):
        """initialize complementary parameters for iss"""
        self.group_no_interf = GROUP_NO_INTERF
        self.__add_constraint_iss()

    def __add_constraint_iss(self):
        """Get global bounderis conditions"""
        __cliss = AFFE_CHAR_MECA(
            MODELE=self.getModele(),
            DDL_IMPO=_F(GROUP_NO=self.group_no_interf, BLOCAGE=("DEPLACEMENT", "ROTATION")),
        )
        self.keywords["CHARGE"] = (
            tuple(list(self.keywords["CHARGE"]) + [__cliss])
            if "CHARGE" in self.keywords
            else (__cliss,)
        )

    def getModele(self):
        """return the model"""
        return self.keywords["MODELE"]

    def getMaillage(self):
        """return the mesh"""
        return self.getModele().getMesh()

    def getNumeddl(self):
        """return dof numbering associated to elementary stiffnesses"""
        if hasattr(self, "_DynaLineFEM__numeddl"):
            return self.__numeddl
        __numeddl = NUME_DDL(MATR_RIGI=self.__getRigielem())
        self.__numeddl = __numeddl
        return self.__numeddl

    def __getRigielem(self):
        """return elementary stiffnesses"""
        if hasattr(self, "_DynaLineFEM__rigielem"):
            return self.__rigielem
        __rigielem = CALC_MATR_ELEM(OPTION="RIGI_MECA", **self.keywords)
        if self.amortissement and self.amortissement["TYPE_AMOR"] == "HYST":
            __rigielem = CALC_MATR_ELEM(
                OPTION="RIGI_MECA_HYST", RIGI_MECA=__rigielem, **self.keywords
            )
        self.__rigielem = __rigielem
        return self.__rigielem

    def __getMasselem(self):
        """return elementary masses"""
        if hasattr(self, "_DynaLineFEM__masselem"):
            return self.__masselem
        __masselem = CALC_MATR_ELEM(OPTION="MASS_MECA", **self.keywords)
        self.__masselem = __masselem
        return self.__masselem

    def __getAmorelem(self):
        """return elementary damping"""
        if hasattr(self, "_DynaLineFEM__amorelem"):
            return self.__amorelem
        if "CHAM_MATER" in self.keywords:
            # RIGI_MECA and MASS_MECA needed if CHAM_MATER for CALC_MATR_ELEM option AMOR_MECA
            keywords = self.keywords.copy()
            keywords["RIGI_MECA"] = self.__getRigielem()
            keywords["MASS_MECA"] = self.__getMasselem()
        else:
            keywords = self.keywords
        if self.amortissement["MATR_IMPE_PHI"] == "OUI":
            amor_flui = "NON"
        else:
            amor_flui = "OUI"
        __amorelem = CALC_MATR_ELEM(
            OPTION="AMOR_MECA", **keywords, AMOR_FLUI=amor_flui, VNOR=self.amortissement["VNOR"]
        )
        self.__amorelem = __amorelem
        return self.__amorelem

    def __getImpeelem(self):
        """return elementary impedances"""
        if hasattr(self, "_DynaLineFEM__impeelem"):
            return self.__impeelem
        if self.amortissement != None:
            if self.amortissement["MATR_IMPE_PHI"]:
                if self.amortissement["MATR_IMPE_PHI"] == "OUI":
                    __impeelem = CALC_MATR_ELEM(
                        OPTION="IMPE_MECA", **self.keywords, VNOR=self.amortissement["VNOR"]
                    )
                    self.__impeelem = __impeelem
        return self.__impeelem

    def getRigiPhy(self):
        """return assembled stiffness matrix over physical basis"""
        if hasattr(self, "_DynaLineFEM__rigiPhy"):
            return self.__rigiPhy
        __rigiPhy = ASSE_MATRICE(
            MATR_ELEM=self.__getRigielem(), NUME_DDL=self.getNumeddl(), **self.char_cine
        )
        self.__rigiPhy = __rigiPhy
        return self.__rigiPhy

    def getMassPhy(self):
        """return assembled mass matrix over physical basis"""
        if hasattr(self, "_DynaLineFEM__massPhy"):
            return self.__massPhy
        __massPhy = ASSE_MATRICE(
            MATR_ELEM=self.__getMasselem(), NUME_DDL=self.getNumeddl(), **self.char_cine
        )
        self.__massPhy = __massPhy
        return self.__massPhy

    def __getAmorPhy(self):
        """return assembled damping matrix over physical basis"""
        if hasattr(self, "_DynaLineFEM__amorPhy"):
            return self.__amorPhy
        if not self.amortissement or self.amortissement["TYPE_AMOR"] != "RAYLEIGH":
            # damping matrix only if RAYLEIGH is defined
            self.__amorPhy = None
            return self.__amorPhy
        __amorPhy = ASSE_MATRICE(
            MATR_ELEM=self.__getAmorelem(), NUME_DDL=self.getNumeddl(), **self.char_cine
        )
        self.__amorPhy = __amorPhy
        return self.__amorPhy

    def __getImpePhy(self):
        """return assembled impedance matrix over physical basis"""
        if hasattr(self, "_DynaLineFEM__impePhy"):
            return self.__impePhy
        self.__impePhy = None
        # if "CHARGE" in self.keywords:
        # for charge in self.keywords["CHARGE"]:
        # if charge.hasLoadField("IMPED"):
        if self.amortissement != None:
            if self.amortissement["MATR_IMPE_PHI"]:
                if self.amortissement["MATR_IMPE_PHI"] == "OUI":
                    __impePhy = ASSE_MATRICE(
                        MATR_ELEM=self.__getImpeelem(), NUME_DDL=self.getNumeddl(), **self.char_cine
                    )
                    self.__impePhy = __impePhy
        return self.__impePhy

    def __getRigiGen(self):
        """return assembled stiffness matrix over gene basis"""
        if hasattr(self, "_DynaLineFEM__rigiGen"):
            return self.__rigiGen
        __rigiGen = self.dynaLineBasis.matrPhyToGen(self.getRigiPhy())
        self.__rigiGen = __rigiGen
        return self.__rigiGen

    def __getMassGen(self):
        """return assembled mass matrix over gene basis"""
        if hasattr(self, "_DynaLineFEM__massGen"):
            return self.__massGen
        __massGen = self.dynaLineBasis.matrPhyToGen(self.getMassPhy())
        if self.dynaLineBasis.getAddedMass():
            __massGen = COMB_MATR_ASSE(
                COMB_R=(
                    _F(MATR_ASSE=__massGen, COEF_R=1.0),
                    _F(MATR_ASSE=self.dynaLineBasis.getAddedMass(), COEF_R=1.0),
                )
            )
        self.__massGen = __massGen
        return self.__massGen

    def __getAmorGen(self):
        """return assembled damping matrix over gene basis"""
        if hasattr(self, "_DynaLineFEM__amorGen"):
            return self.__amorGen
        __amorGen = self.dynaLineBasis.matrPhyToGen(self.__getAmorPhy())
        self.__amorGen = __amorGen
        return self.__amorGen

    def __getImpeGen(self):
        """return assembled impedance matrix over gene basis"""
        if hasattr(self, "_DynaLineFEM__impeGen"):
            return self.__impeGen
        __impeGen = self.dynaLineBasis.matrPhyToGen(self.__getImpePhy())
        self.__impeGen = __impeGen
        return self.__impeGen

    def getMass(self):
        """return assembled mass matrix"""
        if hasattr(self, "_DynaLineFEM__mass"):
            return self.__mass
        if self.__isBasePhys:
            self.__mass = self.getMassPhy()
        if self.__isBaseGene:
            self.__mass = self.__getMassGen()
        return self.__mass

    def getRigi(self):
        """return assembled stiffness matrix"""
        if hasattr(self, "_DynaLineFEM__rigi"):
            return self.__rigi
        if self.__isBasePhys:
            self.__rigi = self.getRigiPhy()
        if self.__isBaseGene:
            self.__rigi = self.__getRigiGen()
        return self.__rigi

    def getAmor(self):
        """return assembled damping matrix"""
        if hasattr(self, "_DynaLineFEM__amor"):
            return self.__amor
        if self.__isBasePhys:
            self.__amor = self.__getAmorPhy()
        if self.__isBaseGene:
            self.__amor = self.__getAmorGen()
        return self.__amor

    def getImpe(self):
        """return assembled impedance matrix"""
        if hasattr(self, "_DynaLineFEM__impe"):
            return self.__impe
        if self.__isBasePhys:
            self.__impe = self.__getImpePhy()
        if self.__isBaseGene:
            self.__impe = self.__getImpeGen()
        return self.__impe

    def getAmorModal(self):
        """return modal damping"""
        if hasattr(self, "_DynaLineFEM__amorModal"):
            return self.__amorModal
        if not self.amortissement or self.amortissement["TYPE_AMOR"] != "MODAL":
            self.__amorModal = None
            return self.__amorModal
        self.__amorModal = self.amortissement["AMOR_REDUIT"]
        return self.__amorModal

    def convertBlocToVect(self):
        """compute and return the assembled vector corresponding to blockage"""
        keywords = self.keywords.copy()
        if "CHARGE" in keywords:
            del keywords["MODELE"]
            keywords["CHARGE"] = [x for x in keywords["CHARGE"] if not x.hasLoadField("IMPED")]
            if len(keywords["CHARGE"]) == 0:
                return None
            __vectelem = CALC_VECT_ELEM(OPTION="CHAR_MECA", **keywords)
            __vect = ASSE_VECTEUR(VECT_ELEM=__vectelem, NUME_DDL=self.getNumeddl())
            return __vect
        else:
            return None

    def convertChargeToVect(self, charge):
        """compute and return the assembled vector corresponding to 'charge'"""
        keywords = self.keywords.copy()
        keywords["CHARGE"] = [charge]
        del keywords["MODELE"]
        __vectelem = CALC_VECT_ELEM(OPTION="CHAR_MECA", **keywords)
        __vect = ASSE_VECTEUR(VECT_ELEM=__vectelem, NUME_DDL=self.getNumeddl())
        return __vect

    def asseToNumeddl(self, vect):
        """return the vect result of CREA_CHAMP(ASSE) using self.__numeddl"""
        __vect = CREA_CHAMP(
            TYPE_CHAM="NOEU_DEPL_R",
            MODELE=self.getModele(),
            NUME_DDL=self.getNumeddl(),
            OPERATION="ASSE",
            ASSE=_F(TOUT="OUI", CHAM_GD=vect),
        )
        return __vect

    def setBasis(self, dynaLineBasis):
        """set a basis to the model"""
        self.dynaLineBasis = dynaLineBasis

    def __getRigiPhyInv(self):
        """call FACTORISER to return the inverse stiffness matrix"""
        if hasattr(self, "_DynaLineFEM__rigiPhyInv"):
            return self.__rigiPhyInv
        __rigiPhyInv = FACTORISER(
            reuse=self.getRigiPhy(),
            MATR_ASSE=self.getRigiPhy(),
            STOP_SINGULIER="OUI",
            METHODE="MULT_FRONT",
            NPREC=8,
        )
        self.__rigiPhyInv = __rigiPhyInv
        return self.__rigiPhyInv

    def getActiveDofsForNodesOrNodeGroups(
        self, nodes_or_node_groups, filter_dofs=("DX", "DY", "DZ")
    ):
        """retrieve active dofs for nodes or node groups containded in nodes_or_node_groups
        A dof is considered as active for the current model, if a force_nodale imply a non
        null displacement in the dof direction for all nodes or node groups contained in nodes_or_node_groups
        """
        active_dofs = []
        for dof in filter_dofs:
            force_nodale = nodes_or_node_groups.copy()
            if "R" in dof:
                force_nodale["M" + dof[-1]] = 1.0
            else:
                force_nodale["F" + dof[-1]] = 1.0
            try:
                __ch = AFFE_CHAR_MECA(MODELE=self.getModele(), FORCE_NODALE=force_nodale)
            except AsterError as err:
                continue

            __vectelem = CALC_VECT_ELEM(OPTION="CHAR_MECA", CHARGE=__ch)
            __vect = ASSE_VECTEUR(VECT_ELEM=__vectelem, NUME_DDL=self.getNumeddl())

            __champ = RESOUDRE(MATR=self.__getRigiPhyInv(), CHAM_NO=__vect)
            vWithD = __champ.getValuesWithDescription()
            values = vWithD[0]
            cmps = vWithD[1][1]
            values_f = []
            for icount, cmp in enumerate(cmps):
                if cmp[:5] != "LAGR:":
                    values_f.append(values[icount])
            py = np.array(values_f)
            if max(abs(py)) > np.finfo(float).eps:
                active_dofs.append(dof)
        assert active_dofs
        return active_dofs


class DynaLineExcit:
    def __init__(self, dynaLineFEM, TYPE_CALCUL, BASE_CALCUL, EXCIT=None, ISS="NON", **args):
        self.dynaLineFEM = dynaLineFEM
        if EXCIT:
            self.charges = force_list(EXCIT)
        else:
            self.charges = []
        self.__isTypeTran = TYPE_CALCUL == "TRAN"
        self.__isTypeHarm = TYPE_CALCUL == "HARM"
        assert self.__isTypeTran ^ self.__isTypeHarm
        self.__isBaseGene = BASE_CALCUL == "GENE"
        self.__isBasePhys = BASE_CALCUL == "PHYS"
        # Excit always computed over physical basis if ISS
        if ISS == "OUI":
            self.__isBaseGene = False
            self.__isBasePhys = True
        assert self.__isBaseGene ^ self.__isBasePhys

    def __getCharMecaLoadings(self):
        """return the list of all loadings with CHARGE defined"""
        if hasattr(self, "_DynaLineExcit__charMecaLoadings"):
            return self.__charMecaLoadings
        self.__charMecaLoadings = [x for x in self.charges if "CHARGE" in x]
        return self.__charMecaLoadings

    def __getMultiAppuiLoadings(self):
        """return the list of all loadings with TYPE_APPUI='MULTI'"""
        if hasattr(self, "_DynaLineExcit__multiAppuiLoadings"):
            return self.__multiAppuiLoadings
        self.__multiAppuiLoadings = [
            x for x in self.charges if "TYPE_APPUI" in x and x["TYPE_APPUI"] == "MULTI"
        ]
        return self.__multiAppuiLoadings

    def getMonoAppuiLoadings(self):
        """return the list of all loadings with TYPE_APPUI='MONO'"""
        if hasattr(self, "_DynaLineExcit__monoAppuiLoadings"):
            return self.__monoAppuiLoadings
        self.__monoAppuiLoadings = [
            x for x in self.charges if "TYPE_APPUI" in x and x["TYPE_APPUI"] == "MONO"
        ]
        return self.__monoAppuiLoadings

    def __check(self):
        """check consistency of EXCIT, raise error if not consistent"""
        if len(self.getMonoAppuiLoadings()) > 1:
            monoAppuiLoadingsWithoutGrNoOrNodeKey = [
                x for x in self.getMonoAppuiLoadings() if not ("GROUP_NO" in x or "NOEUD" in x)
            ]
            monoAppuiLoadingsWithGrNoKey = [
                x for x in self.getMonoAppuiLoadings() if "GROUP_NO" in x
            ]
            grNoValues = [x["GROUP_NO"] for x in monoAppuiLoadingsWithGrNoKey]
            monoAppuiLoadingsWithNodeKey = [x for x in self.getMonoAppuiLoadings() if "NOEUD" in x]
            nodeValues = [x["NOEUD"] for x in monoAppuiLoadingsWithNodeKey]
            # verify that there is only one loading if NODE and GROUP_NO not defined
            if monoAppuiLoadingsWithoutGrNoOrNodeKey:
                if len(monoAppuiLoadingsWithoutGrNoOrNodeKey) > 1 or grNoValues or nodeValues:
                    raise Exception(
                        "More than one loading with TYPE_APPUI='MONO', "
                        "without GROUP_NO or NOEUD defined, not expected"
                    )
            # verify that all mono appui loadings share the same nodes
            if (
                grNoValues
                and nodeValues
                or grNoValues
                and reduce(lambda x, y: x if x == y else False, grNoValues) is False
                or nodeValues
                and reduce(lambda x, y: x if x == y else False, nodeValues) is False
            ):
                raise Exception(
                    "The nodes or group_no must not change between "
                    "the different keywords factors of the exitation "
                    "in the case of TYPE_APPUI = 'MONO'"
                )
        if len(self.__getMultiAppuiLoadings()) > 0 and len(self.getMonoAppuiLoadings()) > 0:
            raise NotImplementedError("found MULTI_APPUI and MONO_APPUI in EXCIT, not expected")
        if len(self.__getCharMecaLoadings()) + len(self.__getMultiAppuiLoadings()) + len(
            self.getMonoAppuiLoadings()
        ) != len(self.charges):
            raise NotImplementedError(
                "at least one charge in EXCIT has not been taken into account, not expected"
            )

    def __setVectOrVectGeneToCharge(self, charge, __vect):
        if self.__isBasePhys:
            charge["VECT_ASSE"] = __vect
        if self.__isBaseGene:
            __vectGen = self.dynaLineFEM.dynaLineBasis.vectPhyToGen(__vect)
            charge["VECT_ASSE_GENE"] = __vectGen

    def get(self):
        """return all loading after converted to 'vect_asse' or 'vect_asse_gene'"""
        if hasattr(self, "_DynaLineExcit__charges"):
            return self.__charges
        self.__check()
        self.__charges = []
        # add blockages
        __vect = self.dynaLineFEM.convertBlocToVect()
        if __vect is not None:
            charge = {"COEF_MULT": 1.0}
            self.__setVectOrVectGeneToCharge(charge, __vect)
            self.__charges.append(charge)
        # kinematic loads
        if self.dynaLineFEM.char_cine != {}:
            charge = {"COEF_MULT": 1.0, "CHARGE": self.dynaLineFEM.char_cine["CHAR_CINE"]}
            self.__charges.append(charge)

        for charMecaLoading in self.__getCharMecaLoadings():
            charge = charMecaLoading.copy()
            __vect = self.dynaLineFEM.convertChargeToVect(charge["CHARGE"])
            self.__setVectOrVectGeneToCharge(charge, __vect)
            del charge["CHARGE"]
            self.__charges.append(charge)
        for multiAppuiLoading in self.__getMultiAppuiLoadings():
            charge = multiAppuiLoading.copy()
            keywords = {}
            for key in ["GROUP_NO", "NOEUD", "DIRECTION"]:
                if key in charge:
                    keywords[key] = charge[key]
            __vect = CALC_CHAR_SEISME(
                MATR_MASS=self.dynaLineFEM.getMassPhy(),
                MODE_STAT=self.dynaLineFEM.dynaLineBasis.getStaticMultiModes(),
                **keywords
            )
            self.__setVectOrVectGeneToCharge(charge, __vect)
            if self.__isTypeTran:
                charge["MULT_APPUI"] = "OUI"
            if self.__isTypeHarm:
                for key in list(keywords.keys()):
                    del charge[key]
            del charge["TYPE_APPUI"]
            self.__charges.append(charge)

        for monoAppuiLoading in self.getMonoAppuiLoadings():
            charge = monoAppuiLoading.copy()
            if "GROUP_NO" in charge or "NOEUD" in charge:
                keywords = {}
                for key in ["GROUP_NO", "NOEUD", "DIRECTION"]:
                    if key in charge:
                        keywords[key] = charge[key]
                __vect = CALC_CHAR_SEISME(
                    MATR_MASS=self.dynaLineFEM.getMassPhy(),
                    MODE_STAT=self.dynaLineFEM.dynaLineBasis.getStaticMonoModes(),
                    **keywords
                )
            else:
                __vect = CALC_CHAR_SEISME(
                    MATR_MASS=self.dynaLineFEM.getMassPhy(),
                    MONO_APPUI="OUI",
                    DIRECTION=charge["DIRECTION"],
                )

            self.__setVectOrVectGeneToCharge(charge, __vect)

            for key in ["GROUP_NO", "NOEUD", "DIRECTION"]:
                if key in charge:
                    del charge[key]
            del charge["TYPE_APPUI"]
            self.__charges.append(charge)
        self.__charges.extend(self.dynaLineFEM.dynaLineBasis.getAddedForces())
        return self.__charges


class DynaLineFrequencyBand:
    """hold frequency band parameters"""

    def __init__(self, PCENT_COUP, TYPE_CALCUL, EXCIT=None, BANDE_ANALYSE=None, **args):
        if EXCIT:
            self.charges = force_list(EXCIT)
        else:
            self.charges = []
        self.pcent_coup = PCENT_COUP
        self.bande_analyse = BANDE_ANALYSE
        self.type_calcul = TYPE_CALCUL

    def __isTypeCalculHarm(self):
        """return True if type_calcul is HARM"""
        return self.type_calcul == "HARM"

    def getFc(self):
        """return cutoff frequency.
        Cutoff frequency will be set to F_max given by self.bande_analyse if self.bande_analyse is defined.
        Else, cutoff frequency will automaticaly be computed using self.pcent_coup"""
        if hasattr(self, "_DynaLineIncrement__fc"):
            return self.__fc
        if self.bande_analyse:
            self.__fc = self.bande_analyse[-1]
            return self.__fc
        if not self.charges:
            raise Exception(
                "unable to compute automatic cutoff frequency without EXCIT factor keyword"
            )
        l_freq_fc = []
        for charge in self.charges:
            if self.__isTypeCalculHarm():
                if "FONC_MULT" in charge:
                    fonction = charge["FONC_MULT"]
                elif "FONC_MULT_C" in charge:
                    fonction = charge["FONC_MULT_C"]
                else:
                    raise NotImplementedError(
                        "Expected to have FONC_MULT or FONC_MULT_C defined in charge for automatic computation of cutoff frequency"
                    )
                if fonction.getType()[:7] == "FORMULE":
                    raise Exception(
                        "automatic computation of cutoff frequency is not possible if a CHARGE is defined with a FORMULE. "
                        "A FONCTION should be pass as FONC_MULT instead, you can use CALC_FONC_INTERP to convert your FORMULE to a discret FONCTION"
                    )
                l_freq_fc.append(max(fonction.Absc()))
            else:
                if "FONC_MULT" in charge:
                    fonction = charge["FONC_MULT"]
                elif "ACCE" in charge:
                    fonction = charge["ACCE"]
                else:
                    raise NotImplementedError(
                        "Expected to have FONC_MULT or ACCE defined in charge"
                    )
                if fonction.getType()[:7] == "FORMULE":
                    raise Exception(
                        "automatic computation of cutoff frequency is not possible if a CHARGE is defined with a FORMULE. "
                        "A FONCTION should be pass as FONC_MULT instead, you can use CALC_FONC_INTERP to convert your FORMULE to a discret FONCTION"
                    )
                # compute fft
                __fonc = CALC_FONCTION(FFT=_F(FONCTION=fonction))
                len_fonc = len(__fonc.Absc()) // 2  # remove negative frequencies of the spectrum
                l_freq = __fonc.Absc()[:len_fonc]
                l_real2 = [x**2 for x in __fonc.Ordo()[:len_fonc]]
                l_imag2 = [x**2 for x in __fonc.OrdoImg()[:len_fonc]]
                l_norm2 = [x + y for x, y in zip(l_real2, l_imag2)]
                # retrieve cutoff value
                l_serie = []
                last_value = 0
                for norm2 in l_norm2:
                    last_value += norm2
                    l_serie.append(last_value)
                cutoff_value = last_value * self.pcent_coup / 100.0
                # retrieve corresponding cutoff frequency
                i = 0
                while l_serie[i] < cutoff_value and i < len(l_serie) - 1:
                    i += 1
                l_freq_fc.append(l_freq[i])
        # retrieve the max of cutoff frequency
        self.__fc = max(l_freq_fc)
        return self.__fc

    def __getFmax(self):
        """return Fmax as last argument of self.bande_analyse if given, else return cutoff frequency times 2"""
        if self.bande_analyse:
            return self.bande_analyse[-1]
        return 2 * self.getFc()

    def __getFmin(self):
        """return Fmin as first value of self.bande_analyse if given, else return 0"""
        if self.bande_analyse and len(self.bande_analyse) == 2:
            return self.bande_analyse[0]
        return 0.0

    def get(self):
        """return the frequency band"""
        return (self.__getFmin(), self.__getFmax())


class DynaLineBasis:
    """hold the calculation basis"""

    dofs = ["DX", "DY", "DZ", "DRX", "DRY", "DRZ"]
    group_no_keys = ["GROUP_NO", "GROUP_NO_1", "GROUP_NO_2", "GROUP_NO_D", "GROUP_NO_G"]
    noeud_keys = ["NOEUD", "NOEUD_1", "NOEUD_2", "NOEUD_D", "NOEUD_G"]

    def __init__(
        self,
        parent,
        dynaLineFEM,
        dynaLineFrequencyBand,
        EXCIT=None,
        BASE_RESU=None,
        ENRI_STAT="NON",
        ORTHO="NON",
        FREQ_COUP_STAT=None,
        ISS="NON",
        IFS="NON",
        COMPORTEMENT=None,
        **args
    ):
        self.parent = parent
        self.base_resu = BASE_RESU
        self.dynaLineFEM = dynaLineFEM
        # cross reference dynaLineBasis and dynaLineFEM
        self.dynaLineFEM.setBasis(self)
        self.dynaLineFrequencyBand = dynaLineFrequencyBand
        if EXCIT:
            self.charges = force_list(EXCIT)
        else:
            self.charges = []
        self.enri_stat = ENRI_STAT == "OUI"
        self.ortho = ORTHO == "OUI"
        self.fc_stat = FREQ_COUP_STAT
        self.iss = ISS == "OUI"
        self.ifs = IFS == "OUI"
        if COMPORTEMENT:
            self.comportement = force_list(COMPORTEMENT)
        else:
            self.comportement = None
        if self.iss:
            self.__init_iss(**args)
        if self.ifs:
            self.__init_ifs(**args)

    def __init_iss(self, GROUP_NO_INTERF=None, TYPE_MODE=None, NB_MODE_INTERF=None, **args):
        """initialize complementary parameters for iss"""
        self.group_no_interf = GROUP_NO_INTERF
        self.type_mode_interf = TYPE_MODE
        self.nb_mode_interf = NB_MODE_INTERF

    def __init_ifs(
        self,
        MODELISATION_FLU=None,
        GROUP_MA_FLUIDE=None,
        GROUP_MA_INTERF=None,
        RHO_FLUIDE=None,
        PRESSION_FLU_IMPO=None,
        FORC_AJOU=None,
        **args
    ):
        """initialize complementary parameters for ifs"""
        self.modelisation_flu = MODELISATION_FLU
        self.group_ma_fluide = GROUP_MA_FLUIDE
        self.group_ma_interf = GROUP_MA_INTERF
        self.rho_fluide = force_list(RHO_FLUIDE)[0]
        self.pression_flu_impo = force_list(PRESSION_FLU_IMPO)[0]
        self.forc_ajou = FORC_AJOU == "OUI"

    def matrPhyToGen(self, matrix):
        """project a given matrix to the gen basis"""
        if matrix is None:
            return None
        __matrixGen = PROJ_MATR_BASE(
            BASE=self.get(), NUME_DDL_GENE=self.__getNumeDdlGene(), MATR_ASSE=matrix
        )
        return __matrixGen

    def vectPhyToGen(self, vector, typVect="FORC"):
        """project a given vector to the gen basis"""
        __vectorGen = PROJ_VECT_BASE(
            BASE=self.get(),
            NUME_DDL_GENE=self.__getNumeDdlGene(),
            VECT_ASSE=vector,
            TYPE_VECT=typVect,
        )
        return __vectorGen

    def __getNumeDdlGene(self):
        """retrieve NUME_DDL_GENE corresponding to this basis"""
        if hasattr(self, "_DynaLineBasis__numeDdlGene"):
            return self.__numeDdlGene
        __numeDdlGene = NUME_DDL_GENE(BASE=self.get(), STOCKAGE=self.getStockageType())
        self.__numeDdlGene = __numeDdlGene
        return self.__numeDdlGene

    def __getPseudoModesForLoadingType(self, loading_type):
        """perform a MODE_STATIQUE with PSEUDO_MODE to retrieve the static modes\
        associated to loading_type TYPE_APPUI. loading_type should be in ['MONO', 'MULTI']"""
        if not self.enri_stat:
            return None
        charges = [x for x in self.charges if "TYPE_APPUI" in x and x["TYPE_APPUI"] == loading_type]
        if len(charges) == 0:
            return None
        pseudo_mode = []
        keywords = {}
        keywords["NOM_DIR"] = "AxePseudo_Mode"
        for charge in charges:
            keywords["DIRECTION"] = charge["DIRECTION"]
            pseudo_mode.append(keywords)
        __pseudoModes = MODE_STATIQUE(
            MATR_MASS=self.dynaLineFEM.getMassPhy(),
            MATR_RIGI=self.dynaLineFEM.getRigiPhy(),
            PSEUDO_MODE=pseudo_mode,
        )
        return __pseudoModes

    def __getPseudoMonoModes(self):
        """return the static pseudo modes associated to loading TYPE_APPUI='MONO'"""
        if hasattr(self, "_DynaLineBasis__pseudoMonoModes"):
            return self.__pseudoMonoModes
        self.__pseudoMonoModes = self.__getPseudoModesForLoadingType("MONO")
        return self.__pseudoMonoModes

    def __getPseudoMultiModes(self):
        """return the static pseudo modes associated to loading TYPE_APPUI='MULTI'"""
        if hasattr(self, "_DynaLineBasis__pseudoMultiModes"):
            return self.__pseudoMultiModes
        self.__pseudoMultiModes = self.__getPseudoModesForLoadingType("MULTI")
        return self.__pseudoMultiModes

    def __getStaticModesForLoadingType(self, loading_type):
        """perform a MODE_STATIQUE with MODE_STAT to retrieve the static modes\
        associated to loading_type TYPE_APPUI. loading_type should be in ['MONO', 'MULTI']"""
        charges = [x for x in self.charges if "TYPE_APPUI" in x and x["TYPE_APPUI"] == loading_type]
        if len(charges) == 0:
            return None
        mode_stat = []
        for charge in charges:
            keywords = {}
            for key in ["GROUP_NO", "NOEUD"]:
                if key in charge:
                    keywords[key] = charge[key]
            avec_cmp = []
            for i, dof_value in enumerate(charge["DIRECTION"]):
                if dof_value != 0:
                    avec_cmp.append(self.dofs[i])
            if len(avec_cmp) == 0:
                raise Exception(
                    "No component found in DIRECTION for a TYPE_APPUI='%s' load" % loading_type
                )
            keywords["AVEC_CMP"] = avec_cmp
            mode_stat.append(keywords)
        __staticModes = MODE_STATIQUE(
            MATR_MASS=self.dynaLineFEM.getMassPhy(),
            MATR_RIGI=self.dynaLineFEM.getRigiPhy(),
            MODE_STAT=mode_stat,
        )
        return __staticModes

    def getStaticMonoModes(self):
        """return the static modes associated to loading TYPE_APPUI='MONO'"""
        if hasattr(self, "_DynaLineBasis__staticMonoModes"):
            return self.__staticMonoModes
        self.__staticMonoModes = self.__getStaticModesForLoadingType("MONO")
        return self.__staticMonoModes

    def getStaticMultiModes(self):
        """return the static modes associated to loading TYPE_APPUI='MULTI'"""
        if hasattr(self, "_DynaLineBasis__staticMultiModes"):
            return self.__staticMultiModes
        self.__staticMultiModes = self.__getStaticModesForLoadingType("MULTI")
        return self.__staticMultiModes

    def __orthogonalize(self, __modalBasis):
        """orthogonalize __modalBasis using the following procedure and return the new modal basis:
        1) PROJ_MATR_BASE : condensate mass and stiffness over given __modalBasis
        2) CALC_MODE : compute new modes using condensate matrices
        3) REST_GENE_PHY : getting back new modes over physical basis"""
        __numeDdlGene = NUME_DDL_GENE(BASE=__modalBasis)
        __massGen = PROJ_MATR_BASE(
            BASE=__modalBasis, NUME_DDL_GENE=__numeDdlGene, MATR_ASSE=self.dynaLineFEM.getMassPhy()
        )
        __rigiGen = PROJ_MATR_BASE(
            BASE=__modalBasis, NUME_DDL_GENE=__numeDdlGene, MATR_ASSE=self.dynaLineFEM.getRigiPhy()
        )

        if self.fc_stat:
            freq = self.dynaLineFrequencyBand.get()
            args = _F(OPTION="BANDE", CALC_FREQ=_F(FREQ=(freq[0], self.fc_stat)))
        else:
            args = _F(OPTION="TOUT", SOLVEUR_MODAL=_F(METHODE="QZ"))
        __calc_modes = CALC_MODES(
            MATR_MASS=__massGen,
            MATR_RIGI=__rigiGen,
            VERI_MODE=_F(STOP_ERREUR="NON"),
            STOP_BANDE_VIDE="NON",
            **args
        )
        if self.base_resu:
            nombase = REST_GENE_PHYS(RESU_GENE=__calc_modes, MODE_MECA=__modalBasis)
            self.parent.register_result(nombase, self.base_resu)
            __modalBasis = nombase
        else:
            __modalBasis = REST_GENE_PHYS(RESU_GENE=__calc_modes, MODE_MECA=__modalBasis)
        return __modalBasis

    def __getDynaModes(self):
        """perform a modal analysis to retrieve the dynamic modes"""
        if hasattr(self, "_DynaLineBasis__dynaModes"):
            return self.__dynaModes
        __dynaModes = CALC_MODES(
            OPTION="BANDE",
            VERI_MODE=_F(STOP_ERREUR="NON"),
            STOP_BANDE_VIDE="NON",
            CALC_FREQ=_F(FREQ=self.dynaLineFrequencyBand.get()),
            MATR_MASS=self.dynaLineFEM.getMassPhy(),
            MATR_RIGI=self.dynaLineFEM.getRigiPhy(),
        )
        if self.ifs:
            __numeDdlGene = NUME_DDL_GENE(BASE=__dynaModes, STOCKAGE=self.getStockageType())
            res = MACRO_MATR_AJOU(
                __use_namedtuple__=True,
                MAILLAGE=self.dynaLineFEM.getMaillage(),
                MODELISATION=self.modelisation_flu,
                GROUP_MA_FLUIDE=self.group_ma_fluide,
                GROUP_MA_INTERF=self.group_ma_interf,
                FLUIDE=self.rho_fluide,
                DDL_IMPO=self.pression_flu_impo,
                MODE_MECA=__dynaModes,
                NUME_DDL_GENE=__numeDdlGene,
                MATR_MASS_AJOU=CO("addedMass"),
            )
            __rigiGen = PROJ_MATR_BASE(
                BASE=__dynaModes,
                NUME_DDL_GENE=__numeDdlGene,
                MATR_ASSE=self.dynaLineFEM.getRigiPhy(),
            )
            __massGen = PROJ_MATR_BASE(
                BASE=__dynaModes,
                NUME_DDL_GENE=__numeDdlGene,
                MATR_ASSE=self.dynaLineFEM.getMassPhy(),
            )
            __massGen = COMB_MATR_ASSE(
                COMB_R=(
                    _F(MATR_ASSE=__massGen, COEF_R=1.0),
                    _F(MATR_ASSE=res.addedMass, COEF_R=1.0),
                )
            )
            __dynaModes = CALC_MODES(
                OPTION="BANDE",
                VERI_MODE=_F(STOP_ERREUR="NON"),
                STOP_BANDE_VIDE="NON",
                CALC_FREQ=_F(FREQ=self.dynaLineFrequencyBand.get()),
                MATR_MASS=__massGen,
                MATR_RIGI=__rigiGen,
            )
            if self.base_resu:
                nombase = REST_GENE_PHYS(RESU_GENE=__dynaModes, TOUT_CHAM="OUI")
                self.parent.register_result(nombase, self.base_resu)
                __dynaModes = nombase
            else:
                __dynaModes = REST_GENE_PHYS(RESU_GENE=__dynaModes, TOUT_CHAM="OUI")
        self.__dynaModes = __dynaModes
        return self.__dynaModes

    def getAddedMass(self):
        """return the added mass computed on self.__modalBasis"""
        return self.__getAddedMassAndForces()[0]

    def getAddedForces(self):
        """return the added forces for IFS with sismic load
        computed on self.__modalBasis"""
        return self.__getAddedMassAndForces()[1]

    def __getAddedMassAndForces(self):
        """return the added mass and forces for IFS with sismic load
        computed on self.__modalBasis"""
        if hasattr(self, "_DynaLineBasisaddedMass"):
            return self.addedMass, self.__excitForcAjou
        if not self.ifs:
            self.addedMass = None
            self.__excitForcAjou = []
            return self.addedMass, self.__excitForcAjou
        keywords = {}
        tmp = []
        if self.forc_ajou:
            monoAppuiLoadings = [
                x for x in self.charges if "TYPE_APPUI" in x and x["TYPE_APPUI"] == "MONO"
            ]
            multiAppuiLoadings = [
                x for x in self.charges if "TYPE_APPUI" in x and x["TYPE_APPUI"] == "MULTI"
            ]
            if monoAppuiLoadings:
                AppuiLoadings = monoAppuiLoadings
                keywords["MONO_APPUI"] = "OUI"
            elif multiAppuiLoadings:
                AppuiLoadings = multiAppuiLoadings
                keywords["MODE_STAT"] = self.getStaticMultiModes()
            if keywords:
                keywords["FORC_AJOU"] = []
                count = 0
                for loading in AppuiLoadings:
                    charge = loading.copy()
                    nameToSave = "adfx" + str(count)
                    __addeforce = CO(nameToSave)
                    d_forc_ajou = {"DIRECTION": charge["DIRECTION"], "VECTEUR": __addeforce}
                    for key in ["NOEUD", "GROUP_NO"]:
                        if key in charge:
                            if "MODE_STAT" in keywords:
                                d_forc_ajou[key] = charge[key]
                            else:
                                del charge[key]
                    if not "MODE_STAT" in keywords:
                        del charge["TYPE_APPUI"]
                        del charge["DIRECTION"]
                    d_excit = {"attrname": nameToSave}
                    for key in charge:
                        d_excit[key] = charge[key]
                    keywords["FORC_AJOU"].append(d_forc_ajou)
                    tmp.append(d_excit)
                    count = count + 1
        res = MACRO_MATR_AJOU(
            __use_namedtuple__=True,
            MAILLAGE=self.dynaLineFEM.getMaillage(),
            MODELISATION=self.modelisation_flu,
            GROUP_MA_FLUIDE=self.group_ma_fluide,
            GROUP_MA_INTERF=self.group_ma_interf,
            FLUIDE=self.rho_fluide,
            DDL_IMPO=self.pression_flu_impo,
            MODE_MECA=self.get(),
            NUME_DDL_GENE=self.__getNumeDdlGene(),
            MATR_MASS_AJOU=CO("addedMass"),
            **keywords
        )
        self.__excitForcAjou = []
        for d_excit in tmp:
            d_excit["VECT_ASSE_GENE"] = getattr(res, d_excit["attrname"])
            del d_excit["attrname"]
            self.__excitForcAjou.append(d_excit)
        self.addedMass = res.addedMass
        return self.addedMass, self.__excitForcAjou

    def __getElasModes(self):
        """perform static linear analysis to retrieve the static elas modes"""
        if hasattr(self, "_DynaLineBasis__elasModes"):
            return self.__elasModes
        if not self.enri_stat:
            self.__elasModes = None
            return self.__elasModes
        if self.iss:
            if self.type_mode_interf == "CRAIG_BAMPTON":
                __elasModes = MODE_STATIQUE(
                    MATR_RIGI=self.dynaLineFEM.getRigiPhy(),
                    MODE_STAT=_F(GROUP_NO=self.group_no_interf, TOUT_CMP="OUI"),
                )
            else:
                if self.type_mode_interf == "PSEUDO":
                    nb_mode_interf = 6
                elif self.type_mode_interf == "INTERF":
                    nb_mode_interf = self.nb_mode_interf
                __elasModes = MODE_STATIQUE(
                    MATR_RIGI=self.dynaLineFEM.getRigiPhy(),
                    MODE_STAT=_F(
                        GROUP_NO=self.group_no_interf, NB_MODE=nb_mode_interf, TOUT_CMP="OUI"
                    ),
                )
        else:
            keywords = self.dynaLineFEM.keywords.copy()
            try:
                char_meca_global = keywords.pop("CHARGE")
            except KeyError:
                char_meca_global = None
            if char_meca_global:
                keywords["CHAR_MECA_GLOBAL"] = char_meca_global
            else:
                keywords["LIAISON_DISCRET"] = "OUI"
            elasCharges = [x.copy() for x in self.charges if "CHARGE" in x]
            if len(elasCharges) == 0:
                self.__elasModes = None
                return self.__elasModes
            for i, charge in enumerate(elasCharges):
                charge["CHAR_MECA"] = charge.pop("CHARGE")
                charge["OPTION"] = "SANS"
                charge["NOM_CAS"] = "CAS_%i" % i
                for key in [
                    "FONC_MULT_C",
                    "COEF_MULT_C",
                    "FONC_MULT",
                    "COEF_MULT",
                    "PHAS_DEG",
                    "PUIS_PULS",
                ]:
                    if key in charge:
                        del charge[key]
            __elasModes = MACRO_ELAS_MULT(
                __use_namedtuple__=True,
                CAS_CHARGE=elasCharges,
                SOLVEUR=_F(METHODE="MULT_FRONT"),
                NUME_DDL=self.dynaLineFEM.getNumeddl(),
                **keywords
            )
        self.__elasModes = __elasModes
        return self.__elasModes

    def __getStaticNonLinearModes(self):
        """return a mode_statique corresponding to non linear behaviours. A MODE_STAT is performed
        with FORCE_NODALE for all components of NOEUD and GROUP_NO involved in self.comportement"""
        if hasattr(self, "_DynaLineBasis__staticNonLinearModes"):
            return self.__staticNonLinearModes
        if not self.enri_stat or not self.comportement:
            self.__staticNonLinearModes = None
            return self.__staticNonLinearModes
        force_nodale = []
        for comportement in self.comportement:
            aster.affiche("MESSAGE", "%s" % comportement)
            nodes = []
            node_groups = []
            for key in self.noeud_keys:
                if key in comportement:
                    nodes.append(comportement[key])
            for key in self.group_no_keys:
                if key in comportement:
                    node_groups.append(comportement[key])

            if nodes:
                dic_node = {}
                dic_node["NOEUD"] = nodes
                if "NOM_CMP" in comportement:
                    dic_node["AVEC_CMP"] = comportement["NOM_CMP"]
                else:
                    dic_node["AVEC_CMP"] = self.dynaLineFEM.getActiveDofsForNodesOrNodeGroups(
                        dic_node
                    )
                force_nodale.append(dic_node)
                aster.affiche("MESSAGE", "%s" % dic_node)

            if node_groups:
                dic_gno = {}
                dic_gno["GROUP_NO"] = node_groups
                if "NOM_CMP" in comportement:
                    dic_gno["AVEC_CMP"] = comportement["NOM_CMP"]
                else:
                    dic_gno["AVEC_CMP"] = self.dynaLineFEM.getActiveDofsForNodesOrNodeGroups(
                        dic_gno
                    )
                force_nodale.append(dic_gno)
                aster.affiche("MESSAGE", "%s" % dic_gno)

        __staticNonLinearModes = MODE_STATIQUE(
            MATR_MASS=self.dynaLineFEM.getMassPhy(),
            MATR_RIGI=self.dynaLineFEM.getRigiPhy(),
            FORCE_NODALE=force_nodale,
        )
        self.__staticNonLinearModes = __staticNonLinearModes
        return self.__staticNonLinearModes

    def getStockageType(self):
        """get stockage type of basis"""
        if self.enri_stat and not self.ortho or self.ifs:
            return "PLEIN"
        return "DIAG"

    def __assembleModeIntf(self, l_mode_intf):
        """Assemble the interface modes in one concept"""
        ichamp = 0
        affe_crea_resu = []
        for mode_intf in l_mode_intf:
            dico = mode_intf.LIST_PARA()
            for nume_ordre in dico["NUME_ORDRE"]:
                ichamp += 1
                __champ = CREA_CHAMP(
                    NUME_ORDRE=nume_ordre,
                    TYPE_CHAM="NOEU_DEPL_R",
                    RESULTAT=mode_intf,
                    NOM_CHAM="DEPL",
                    OPERATION="EXTR",
                )
                affe_crea_resu.append(
                    _F(
                        CHAM_GD=__champ,
                        NOM_CHAM="DEPL",
                        MODELE=self.dynaLineFEM.getModele(),
                        NOM_CAS="CHAMP_%i" % ichamp,
                        # NUME_ORDRE = ichamp,
                        # NUME_MODE = ichamp,
                    )
                )
        __modeintf = CREA_RESU(OPERATION="AFFE", TYPE_RESU="MULT_ELAS", AFFE=affe_crea_resu)
        return __modeintf

    def get(self):
        """return the modal basis"""
        if hasattr(self, "_DynaLineBasis__modalBasis"):
            return self.__modalBasis
        mode_meca = self.__getDynaModes()
        if self.ifs:
            # for ifs, DEFI_BASE_MODALE not working with basis comming from REST_GENE_PHYS
            # workarround: retrieve the basis comming from REST_GENE_PHYS and skip DEFI_BASE_MODALE
            self.__modalBasis = mode_meca
            return self.__modalBasis

        mode_intf = []
        if self.__getPseudoMultiModes():
            mode_intf.append(self.__getPseudoMultiModes())
        if self.__getPseudoMonoModes():
            mode_intf.append(self.__getPseudoMonoModes())
        if self.__getElasModes():
            mode_intf.append(self.__getElasModes())
        if self.__getStaticNonLinearModes():
            mode_intf.append(self.__getStaticNonLinearModes())

        ritz = [{"MODE_MECA": mode_meca}]
        args = {}
        if len(mode_intf) > 1:
            mode_intf = self.__assembleModeIntf(mode_intf)
        if mode_intf:
            ritz.append({"MODE_INTF": mode_intf})
            if not self.iss:
                args["ORTHO"] = "OUI"
                args["MATRICE"] = self.dynaLineFEM.getRigiPhy()
        if self.base_resu and not (self.ortho and mode_intf):
            nombase = DEFI_BASE_MODALE(RITZ=ritz, NUME_REF=self.dynaLineFEM.getNumeddl(), **args)
            self.parent.register_result(nombase, self.base_resu)
            __modalBasis = nombase
        else:
            __modalBasis = DEFI_BASE_MODALE(
                RITZ=ritz, NUME_REF=self.dynaLineFEM.getNumeddl(), **args
            )
            if self.ortho and mode_intf:
                __modalBasis = self.__orthogonalize(__modalBasis)

        self.__modalBasis = __modalBasis
        return self.__modalBasis


class DynaLineIncrement:
    """hold increment parameters"""

    def __init__(
        self, dynaLineFrequencyBand, TYPE_CALCUL, FREQ=None, LIST_FREQ=None, INCREMENT=None, **args
    ):
        self.dynaLineFrequencyBand = dynaLineFrequencyBand
        self.__isTypeTran = TYPE_CALCUL == "TRAN"
        self.__isTypeHarm = TYPE_CALCUL == "HARM"
        assert self.__isTypeTran ^ self.__isTypeHarm
        if INCREMENT:
            assert self.__isTypeTran
            self.increment = force_list(INCREMENT)
        else:
            assert self.__isTypeHarm
            self.increment = {}
            assert FREQ or LIST_FREQ
            if LIST_FREQ:
                FREQ = LIST_FREQ.getValues()
            self.increment["FREQ"] = FREQ

    def __getStep(self):
        """return the time step, according to cutoff frequency, step = 1/(5*fc)"""
        if hasattr(self, "_DynaLineIncrement__step"):
            return self.__step
        self.__step = 1.0 / (5 * self.dynaLineFrequencyBand.getFc())
        return self.__step

    def get(self):
        if self.__isTypeTran and "PAS" not in self.increment[0]:
            # compute and add step to self.increment if not given before returning it
            self.increment[0]["PAS"] = self.__getStep()
        return self.increment


class DynaLineInitialState:
    """hold initial state parameters"""

    def __init__(self, dynaLineFem, BASE_CALCUL, ETAT_INIT=None, **args):
        self.dynaLineFem = dynaLineFem
        if ETAT_INIT:
            self.etat_init = force_list(ETAT_INIT)[0]
        else:
            self.etat_init = None
        self.__isBaseGene = BASE_CALCUL == "GENE"

    def __project(self):
        """project all self.etat_init to gene basis"""
        for key in ["DEPL", "VITE"]:
            if key in self.etat_init and self.etat_init[key].getType()[:7] == "CHAM_NO":
                if self.__isBaseGene:
                    __vectGen = self.dynaLineFem.dynaLineBasis.vectPhyToGen(
                        self.etat_init[key], key
                    )
                    self.etat_init[key] = __vectGen
                else:
                    __vect = self.dynaLineFem.asseToNumeddl(self.etat_init[key])
                    self.etat_init[key] = __vect
        if "RESULTAT" in self.etat_init:
            raise Exception("keyword RESULTAT should be implemented here")
        self.isProjected = True

    def get(self):
        if self.etat_init:
            self.__project()
        return self.etat_init


class DynaLineResu:
    def __init__(
        self,
        parent,
        dynaLineFEM,
        dynaLineExcit,
        dynaLineIncrement,
        dynaLineInitialState,
        TYPE_CALCUL,
        BASE_CALCUL,
        SCHEMA_TEMPS=None,
        ARCHIVAGE=None,
        BASE_RESU=None,
        RESU_GENE=None,
        ISS="NON",
        COMPORTEMENT=None,
        SOLVEUR=None,
        **args
    ):
        self.parent = parent
        self.dynaLineFEM = dynaLineFEM
        self.dynaLineExcit = dynaLineExcit
        self.dynaLineIncrement = dynaLineIncrement
        self.dynaLineInitialState = dynaLineInitialState
        self.type_calcul = TYPE_CALCUL
        self.base_calcul = BASE_CALCUL
        self.base_resu = BASE_RESU
        self.schema_temps = SCHEMA_TEMPS
        self.archivage = ARCHIVAGE
        self.resu_gene = RESU_GENE
        self.iss = ISS == "OUI"
        self.solveur = SOLVEUR

        if self.iss:
            self.__init_iss(**args)
        self.comportement = COMPORTEMENT

    def __init_iss(
        self,
        TABLE_SOL=None,
        MATER_SOL=None,
        UNITE_RESU_IMPE=None,
        UNITE_RESU_FORC=None,
        PARAMETRE=None,
        GROUP_MA_INTERF=None,
        CALC_IMPE_FORC=None,
        VERSION_MISS=None,
        **args
    ):
        """initialize complementary parameters for iss"""
        self.table_sol = TABLE_SOL
        self.mater_sol = MATER_SOL
        assert not (self.table_sol and self.mater_sol)
        self.unite_resu_impe = UNITE_RESU_IMPE
        self.unite_resu_forc = UNITE_RESU_FORC
        if PARAMETRE:
            self.parametre = force_list(PARAMETRE)[0]
        else:
            self.parametre = {}
        self.parametre["TYPE"] = "BINAIRE"
        self.group_ma_interf = GROUP_MA_INTERF
        self.calc_impe_forc = CALC_IMPE_FORC == "OUI"
        self.version_miss = VERSION_MISS

    def __getDynaVibraKeywords(self):
        """return common keywords used for calling DYNA_VIBRA"""
        keywords = {}
        if self.type_calcul == "TRAN":
            keywords["INCREMENT"] = self.dynaLineIncrement.get()
            keywords["INCREMENT"][0]["VERI_PAS"] = "NON"
            keywords["SCHEMA_TEMPS"] = self.schema_temps
            if self.archivage:
                keywords["ARCHIVAGE"] = self.archivage
            if self.dynaLineFEM.dynaLineBasis.getStaticMultiModes():
                # DYNA_VIBRA complains for MODE_STAT keyword in case of MULT_APPUI in EXCIT
                keywords["MODE_STAT"] = self.dynaLineFEM.dynaLineBasis.getStaticMultiModes()
        else:
            keywords.update(self.dynaLineIncrement.get())
        if self.dynaLineFEM.getAmor():
            keywords["MATR_AMOR"] = self.dynaLineFEM.getAmor()
        if self.dynaLineFEM.getImpe():
            keywords["MATR_IMPE_PHI"] = self.dynaLineFEM.getImpe()
        if self.dynaLineFEM.getAmor() and self.dynaLineFEM.getImpe():
            keywords["AMOR_FLUI"] = "NON"
        if self.dynaLineFEM.getAmorModal():
            keywords["AMOR_MODAL"] = {"AMOR_REDUIT": self.dynaLineFEM.getAmorModal()}
        if self.dynaLineInitialState.get():
            keywords["ETAT_INIT"] = self.dynaLineInitialState.get()
        if self.comportement:
            keywords["COMPORTEMENT"] = self.comportement
        if self.solveur:
            keywords["SOLVEUR"] = self.solveur
        if self.dynaLineExcit.get():
            keywords["EXCIT"] = self.dynaLineExcit.get()
        return keywords

    def __getDynaVibra(self):
        """return result from DYNA_VIBRA (no iss)"""
        keywords = self.__getDynaVibraKeywords()
        if self.base_calcul == "PHYS":
            nomresu = DYNA_VIBRA(
                MODELE=self.dynaLineFEM.getModele(),
                TYPE_CALCUL=self.type_calcul,
                BASE_CALCUL=self.base_calcul,
                MATR_MASS=self.dynaLineFEM.getMass(),
                MATR_RIGI=self.dynaLineFEM.getRigi(),
                **keywords
            )
        else:
            if self.resu_gene:
                resgene = DYNA_VIBRA(
                    TYPE_CALCUL=self.type_calcul,
                    BASE_CALCUL=self.base_calcul,
                    MATR_MASS=self.dynaLineFEM.getMass(),
                    MATR_RIGI=self.dynaLineFEM.getRigi(),
                    # TOUT_CHAM='OUI',
                    **keywords
                )
                # retrieve also result in gene basis
                self.parent.register_result(resgene, self.resu_gene)
                __dyna_vibra = resgene
            else:
                __dyna_vibra = DYNA_VIBRA(
                    TYPE_CALCUL=self.type_calcul,
                    BASE_CALCUL=self.base_calcul,
                    MATR_MASS=self.dynaLineFEM.getMass(),
                    MATR_RIGI=self.dynaLineFEM.getRigi(),
                    # TOUT_CHAM='OUI',
                    **keywords
                )
            keywords = {}
            if self.dynaLineFEM.dynaLineBasis.getStaticMultiModes():
                keywords["MULT_APPUI"] = "OUI"
            if len(self.dynaLineExcit.getMonoAppuiLoadings()) == 1:
                monoAppuiLoading = self.dynaLineExcit.getMonoAppuiLoadings()[0]
                keywords["DIRECTION"] = monoAppuiLoading["DIRECTION"]
                if "FONC_MULT" in monoAppuiLoading:
                    keywords["ACCE_MONO_APPUI"] = monoAppuiLoading["FONC_MULT"]
                elif "ACCE" in monoAppuiLoading:
                    keywords["ACCE_MONO_APPUI"] = monoAppuiLoading["ACCE"]
                else:
                    assert False
            nomresu = REST_GENE_PHYS(
                RESU_GENE=__dyna_vibra,
                MODE_MECA=self.dynaLineFEM.dynaLineBasis.get(),
                TOUT_CHAM="OUI",
                **keywords
            )
        return nomresu

    def __getCalcMissKeywords(self):
        """return common keywords used for calling CALC_MISS"""
        keywords = {}
        if self.table_sol:
            keywords["TABLE_SOL"] = self.table_sol
        if self.mater_sol:
            keywords["MATER_SOL"] = self.mater_sol
        if self.unite_resu_impe:
            keywords["UNITE_RESU_IMPE"] = self.unite_resu_impe
        if self.unite_resu_forc:
            keywords["UNITE_RESU_FORC"] = self.unite_resu_forc
        if self.parametre:
            if self.type_calcul == "HARM":
                l_is_freq_defined = [
                    x in list(self.parametre.keys()) for x in ["FREQ_MIN", "LIST_FREQ", "FREQ_IMAG"]
                ]
                is_freq_defined = reduce(lambda x, y: x or y, l_is_freq_defined)
                if not is_freq_defined:
                    # recuperate frequencies from increment definition if not defined
                    self.parametre["LIST_FREQ"] = self.dynaLineIncrement.get()["FREQ"]
            keywords["PARAMETRE"] = self.parametre
        if self.group_ma_interf:
            keywords["GROUP_MA_INTERF"] = self.group_ma_interf
        keywords["BASE_MODALE"] = self.dynaLineFEM.dynaLineBasis.get()
        if self.version_miss:
            keywords["VERSION"] = self.version_miss
        return keywords

    def __getCalcMissAdditionalKeywords(self):
        """return common keywords used only for the second CALC_MISS call"""
        keywords = {}
        if self.dynaLineFEM.getAmorModal():
            keywords["AMOR_REDUIT"] = self.dynaLineFEM.getAmorModal()
        if self.type_calcul == "HARM":
            # recuperate frequencies from increment definition
            self.parametre["LIST_FREQ"] = self.dynaLineIncrement.get()["FREQ"]
            for key in ["FREQ_MIN", "FREQ_MAX", "FREQ_PAS", "FREQ_IMAG"]:
                if key in self.parametre:
                    del self.parametre[key]
            keywords["PARAMETRE"] = self.parametre
        if self.type_calcul == "HARM" and len(self.dynaLineExcit.get()) > 0:
            keywords["EXCIT_HARMO"] = self.dynaLineExcit.get()
        return keywords

    def __getCalcMiss(self):
        """return a result from CALC_MISS (if iss)"""
        keywords = self.__getCalcMissKeywords()
        # first calc_miss with FICHIER
        if self.calc_impe_forc and self.unite_resu_impe and self.unite_resu_forc:
            CALC_MISS(TYPE_RESU="FICHIER", **keywords)
        keywords.update(self.__getCalcMissAdditionalKeywords())
        if self.resu_gene:
            resu_gene = CALC_MISS(
                TYPE_RESU=self.type_calcul + "_GENE",
                MODELE=self.dynaLineFEM.getModele(),
                MATR_RIGI=self.dynaLineFEM.getRigiPhy(),
                MATR_MASS=self.dynaLineFEM.getMassPhy(),
                **keywords
            )
            # retrieve also result in gene basis
            self.parent.register_result(resgene, self.resu_gene)
            __calc_miss = resu_gene
        else:
            __calc_miss = CALC_MISS(
                TYPE_RESU=self.type_calcul + "_GENE",
                MODELE=self.dynaLineFEM.getModele(),
                MATR_RIGI=self.dynaLineFEM.getRigiPhy(),
                MATR_MASS=self.dynaLineFEM.getMassPhy(),
                **keywords
            )
        keywords = {}
        nomresu = REST_GENE_PHYS(
            RESU_GENE=__calc_miss,
            MODE_MECA=self.dynaLineFEM.dynaLineBasis.get(),
            TOUT_CHAM="OUI",
            **keywords
        )
        return nomresu

    def get(self):
        """return the result for DYNA_LINE"""
        if not self.iss:
            return self.__getDynaVibra()
        else:
            return self.__getCalcMiss()


def dyna_line_ops(self, **args):
    """Ecriture de la macro DYNA_LINE"""

    dynaLineFEM = DynaLineFEM(self, **args)
    dynaLineFrequencyBand = DynaLineFrequencyBand(**args)
    dynaLineBasis = DynaLineBasis(self, dynaLineFEM, dynaLineFrequencyBand, **args)
    dynaLineExcit = DynaLineExcit(dynaLineFEM, **args)
    dynaLineIncrement = DynaLineIncrement(dynaLineFrequencyBand, **args)
    dynaLineInitialState = DynaLineInitialState(dynaLineFEM, **args)

    dynaLineResu = DynaLineResu(
        self, dynaLineFEM, dynaLineExcit, dynaLineIncrement, dynaLineInitialState, **args
    )
    return dynaLineResu.get()
