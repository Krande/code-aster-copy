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

# person_in_charge: mathieu.courtois at edf.fr

"""Commande CALC_FONCTION"""

import math
import os
import traceback

import aster_fonctions
import numpy as NP
from libaster import AsterError

from ..Cata.Syntax import _F
from ..CodeCommands import DEFI_FONCTION, DEFI_NAPPE, IMPR_FONCTION
from ..Messages import ASSERT, UTMESS
from ..Objects import Function, Function2D, FunctionComplex
from ..Objects.function_py import (
    enveloppe,
    fractile,
    homo_support_nappe,
    moyenne,
    t_fonction,
    t_fonction_c,
    t_nappe,
)
from ..Utilities import force_list
from .defi_inte_spec_ops import tocomplex
from .Utils import liss_enveloppe as LISS
from .Utils.calc_coherency import calc_cohefromdata
from .Utils.random_signal_utils import ACCE2SRO, DSP2SRO, SRO2DSP, acce_filtre_CP, f_phase_forte


def calc_fonction_prod(
    DERIVE=None,
    EXTRACTION=None,
    INTEGRE=None,
    INVERSE=None,
    COMB=None,
    COMB_C=None,
    MULT=None,
    ENVELOPPE=None,
    FRACTILE=None,
    PROL_SPEC_OSCI=None,
    SPEC_OSCI=None,
    ASSE=None,
    FFT=None,
    COMPOSE=None,
    CORR_ACCE=None,
    COHERENCE=None,
    PUISSANCE=None,
    LISS_ENVELOP=None,
    ABS=None,
    REGR_POLYNOMIALE=None,
    DSP=None,
    MOYENNE=None,
    INTEGRE_FREQ=None,
    DERIVE_FREQ=None,
    INTERPOL_FFT=None,
    **args
):
    if INTEGRE is not None:
        return Function
    if DERIVE is not None:
        return Function
    if INVERSE is not None:
        return Function
    if COMB is not None:
        comb = COMB[0]["FONCTION"]
        if type(comb) not in (list, tuple):
            type_vale = type(COMB[0]["FONCTION"])
        else:
            type_vale = type(COMB[0]["FONCTION"][0])
        for mcfact in COMB:
            if type(mcfact["FONCTION"]) is not type_vale:
                raise TypeError("CALC_FONCTION/COMB : pas de types hétérogènes nappe/fonction")
        return type_vale
    if COMB_C is not None:
        vale = COMB_C[0]["FONCTION"]
        if type(vale) is list:
            vale = vale[0]
        if type(vale) is Function2D:
            for mcfact in COMB_C[1:]:
                if type(mcfact["FONCTION"]) is not Function2D:
                    raise TypeError(
                        "CALC_FONCTION/COMB_C : pas de types hétérogènes nappe/fonction"
                    )
            return Function2D
        else:
            for mcfact in COMB_C:
                if type(mcfact["FONCTION"]) is Function2D:
                    raise TypeError(
                        "CALC_FONCTION/COMB_C : pas de types hétérogènes nappe/fonction"
                    )
            return FunctionComplex
    if ENVELOPPE is not None:
        if type(ENVELOPPE[0]["FONCTION"]) not in (list, tuple):
            return type(ENVELOPPE[0]["FONCTION"])
        else:
            return type(ENVELOPPE[0]["FONCTION"][0])
    if FRACTILE is not None:
        if type(FRACTILE[0]["FONCTION"]) not in (list, tuple):
            return type(FRACTILE[0]["FONCTION"])
        else:
            return type(FRACTILE[0]["FONCTION"][0])
    if MOYENNE is not None:
        if type(MOYENNE[0]["FONCTION"]) not in (list, tuple):
            return type(MOYENNE[0]["FONCTION"])
        else:
            return type(MOYENNE[0]["FONCTION"][0])
    if EXTRACTION is not None:
        return Function
    if PROL_SPEC_OSCI is not None:
        return Function
    if SPEC_OSCI is not None:
        if SPEC_OSCI[0]["TYPE_RESU"] == "NAPPE":
            return Function2D
        else:
            if SPEC_OSCI[0]["AMOR_REDUIT"] is not None:
                if len(SPEC_OSCI[0]["AMOR_REDUIT"]) == 1:
                    return Function
                else:
                    return Function2D
            else:
                return Function2D
    if DSP is not None:
        return Function
    if COMPOSE is not None:
        return Function
    if ASSE is not None:
        return Function
    if MULT is not None:
        comb = MULT[0]["FONCTION"]
        if type(comb) not in (list, tuple):
            type_vale = type(MULT[0]["FONCTION"])
        else:
            type_vale = type(MULT[0]["FONCTION"][0])
        for mcfact in MULT:
            if type(MULT[0]["FONCTION"]) not in (list, tuple):
                type_test = type(MULT[0]["FONCTION"])
            else:
                type_test = type(MULT[0]["FONCTION"][0])
            if type_test != type_vale:
                raise TypeError("CALC_FONCTION/MULT : pas de types hétérogènes nappe/fonction")
        return type_vale
    if FFT is not None:
        vale = FFT[0]["FONCTION"]
        if type(vale) is Function:
            return FunctionComplex
        if type(vale) is FunctionComplex:
            return Function
    if INTERPOL_FFT is not None:
        return Function
    if CORR_ACCE is not None:
        return Function
    if COHERENCE is not None:
        return Function
    if LISS_ENVELOP is not None:
        return Function2D
    if REGR_POLYNOMIALE is not None:
        return Function
    if PUISSANCE is not None:
        if type(PUISSANCE[0]["FONCTION"]) not in (list, tuple):
            return type(PUISSANCE[0]["FONCTION"])
        else:
            return type(PUISSANCE[0]["FONCTION"][0])
    if ABS is not None:
        return Function
    if INTEGRE_FREQ is not None:
        return Function
    if DERIVE_FREQ is not None:
        return Function
    raise TypeError("type de concept resultat non prevu")


def calc_fonction_ops(self, **args):
    """Corps de la macro CALC_FONCTION"""
    args = _F(args)
    # éléments de contexte
    ctxt = Context()

    operation = CalcFonctionOper.factory(self, ctxt, args)
    try:
        result = operation.run()
    except AsterError as msg:
        UTMESS("F", "FONCT0_30", valk=(ctxt.f, str(msg), traceback.format_exc()))
    return result


class CalcFonctionOper:
    """Base of all CALC_FONCTION operations.

    Subclasses must implement the '_run' method and, if necessary, may
    overload the '_build_data' method.
    """

    @staticmethod
    def factory(macro, ctxt, kwargs):
        """Factory that returns the operation object
        Subclasses must be named 'CalcFonction_KEYWORD'
        """
        types = CalcFonctionOper.__subclasses__()
        try:
            for class_ in types:
                keyw = class_.__name__.replace("CalcFonction_", "")
                if kwargs[keyw]:
                    return class_(macro, keyw, ctxt, kwargs)
        except KeyError:
            UTMESS("F", "DVP_1")

    def __init__(self, macro, oper, ctxt, kwargs):
        """Initialization"""
        self.macro = macro
        self.oper = oper
        self.ctxt = ctxt
        self.args = kwargs
        if type(self.args[self.oper]) in (list, tuple):
            if len(self.args[self.oper]) == 1:
                self.kw = self.args[self.oper][0]
            else:
                self.kw = self.args[self.oper]
        else:
            self.kw = self.args[self.oper]
        self.resu = None
        self._lf = []
        self._dat = None
        self.typres = calc_fonction_prod(**kwargs)

    def _build_data(self):
        """Read keywords to build the data"""
        self._build_list_fonc()

    def _run(self):
        """Use the input functions from 'self._lf' and the data from
        'self._dat' to compute the result function 'self.resu'."""
        raise NotImplementedError("must be defined in a subclass")

    def run(self):
        """Run the operator"""
        self._build_data()
        self.ctxt.f = [func.nom for func in self._lf]
        self._run()
        return self.build_result()

    def build_result(self):
        """Create the result function"""
        # common keywords to DEFI_FONCTION & DEFI_NAPPE
        para = self.resu.para
        for p in ("NOM_PARA", "NOM_RESU", "PROL_DROITE", "PROL_GAUCHE", "INTERPOL"):
            if self.args[p] is not None:
                para[p] = self.args[p]

        if self.typres is not Function2D:
            if self.typres is FunctionComplex:
                mcval = "VALE_C"
            else:
                mcval = "VALE"
            para[mcval] = self.resu.tabul()
            result = DEFI_FONCTION(**para)
        else:
            intf = self.args["INTERPOL_FONC"]
            prdf = self.args["PROL_DROITE_FONC"]
            prgf = self.args["PROL_GAUCHE_FONC"]
            def_fonc = []
            for f_i in self.resu.l_fonc:
                def_fonc.append(
                    _F(
                        VALE=f_i.tabul(),
                        INTERPOL=intf or f_i.para["INTERPOL"],
                        PROL_DROITE=prdf or f_i.para["PROL_DROITE"],
                        PROL_GAUCHE=prgf or f_i.para["PROL_GAUCHE"],
                    )
                )
            npf = "NOM_PARA_FONC"
            if self.args[npf] is not None:
                para[npf] = self.args[npf]
            result = DEFI_NAPPE(PARA=self.resu.vale_para.tolist(), DEFI_FONCTION=def_fonc, **para)
        if self.args["INFO"] > 1:
            IMPR_FONCTION(FORMAT="TABLEAU", UNITE=6, COURBE=_F(FONCTION=result))
        return result

    # utilities
    def _use_list_para(self):
        """Interpolate using LIST_PARA."""
        if self.args["LIST_PARA"] is not None:
            self.resu = self.resu.evalfonc(self.args["LIST_PARA"].getValues())

    def _get_mcsimp(self, mcsimp):
        """Return the list of mcsimp values for all occurrences of mcfact."""
        # only one occurrence of MCFACT or only one value in MCSIMP
        value = []
        if type(self.kw) in (list, tuple):
            kw = self.kw
        else:
            kw = (self.kw,)
        try:
            nbmf = len(kw)
        except TypeError:
            nbmf = 1
        for mcf in kw:
            val = force_list(mcf[mcsimp])
            assert nbmf == 1 or len(val) == 1, (nbmf, val)
            value.extend(val)
        return value

    def _build_list_fonc(self, arg="real", mcsimp="FONCTION"):
        """Return the list of functions under mcfact/mcsimp converted
        as t_fonction objects.
        nappe_sdaster objects are interpolated on the same abscissa."""
        lf_in = self._get_mcsimp(mcsimp)
        all_nap = min([int(i.getType() == "NAPPE_SDASTER") for i in lf_in]) == 1
        if all_nap:
            list_fonc = [tf.convert() for tf in lf_in]
            list_fonc = homo_support_nappe(list_fonc)
        else:
            list_fonc = [tf.convert(arg) for tf in lf_in]
        self._lf = list_fonc


class CalcFonction_ABS(CalcFonctionOper):
    """Return absolute value"""

    def _run(self):
        """ABS"""
        self.resu = self._lf[0].abs()


class CalcFonction_ASSE(CalcFonctionOper):
    """Concatenate two functions"""

    def _run(self):
        """ASSE"""
        assert len(self._lf) == 2, "exactly 2 functions required"
        fo0, fo1 = self._lf
        self.resu = fo0.cat(fo1, self.kw["SURCHARGE"])


class CalcFonction_COMB(CalcFonctionOper):
    """Combinate real functions."""

    def _run(self):
        """COMB_C"""
        coef = self._get_mcsimp("COEF")
        self.resu = 0.0
        for item, cfr in zip(self._lf, coef):
            self.ctxt.f = item.nom
            self.resu = item * cfr + self.resu
        # take the parameters of the first function
        self.resu.para = self._lf[0].para.copy()
        self._use_list_para()


class CalcFonction_COMB_C(CalcFonctionOper):
    """Combinate complex functions."""

    def _build_data(self):
        """Read keywords to build the data"""
        self._build_list_fonc(arg="complex")
        coefr = self._get_mcsimp("COEF_R")
        coefc = self._get_mcsimp("COEF_C")
        self._dat = {"R": coefr, "C": coefc}

    def _run(self):
        """COMB_C"""
        self.resu = 0.0
        for item, cfr, cfc in zip(self._lf, self._dat["R"], self._dat["C"]):
            coef = 1.0
            if cfr is not None:
                coef = complex(cfr)
            elif cfc is not None:
                coef = cfc
                if type(cfc) in (list, tuple):
                    coef = tocomplex(cfc)
            self.ctxt.f = item.nom
            self.resu = item * coef + self.resu
        # take the parameters of the first function
        self.resu.para = self._lf[0].para.copy()
        self._use_list_para()


class CalcFonction_COMPOSE(CalcFonctionOper):
    """Compose two functions"""

    def _build_data(self):
        """Read keywords to build the data"""
        self._dat = (self.kw["FONC_RESU"].convert(), self.kw["FONC_PARA"].convert())

    def _run(self):
        """COMPOSE"""
        fo1, fo2 = self._dat
        self.resu = fo1[fo2]
        self.resu.para["NOM_PARA"] = fo2.para["NOM_PARA"]


class CalcFonction_PROL_SPEC_OSCI(CalcFonctionOper):
    """PROL_SPEC_OSCI"""

    def _run(self):
        """run PROL_SPEC_OSCI"""
        kw = self.kw
        f_in = self._lf[0]
        dmax = kw["DEPL_MAX"]
        vale_freq = self._lf[0].vale_x
        vale_sro_acce = self._lf[0].vale_y / kw["NORME"]
        vale_sro_depl = vale_sro_acce * vale_freq ** (-2)
        f_min = vale_freq[0]
        assert f_min > 0.0
        d_fmin = vale_sro_depl[0]
        pente = (d_fmin - vale_sro_depl[1]) / (f_min - vale_freq[1])
        freq_dmax = f_min + (dmax - d_fmin) / pente

        if d_fmin < dmax:
            if freq_dmax < 0.00:
                dc = d_fmin - pente * f_min
                UTMESS("F", "FONCT0_78", valr=[dmax, 0.0, dc])
            elif pente > 0.00:
                UTMESS("F", "FONCT0_78", valr=[dmax, f_min, d_fmin])
            else:
                pass
        elif d_fmin > dmax:
            if pente < 0.00:
                dc = d_fmin - pente * f_min
                UTMESS("F", "FONCT0_78", valr=[dmax, f_min, d_fmin])
            elif freq_dmax < 0.00:
                dc = d_fmin - pente * f_min
                UTMESS("F", "FONCT0_78", valr=[dmax, 0.0, dc])
            else:
                pass

        vale_sro_depl = list(vale_sro_depl)
        vale_freq = list(vale_freq)
        vale_sro_depl.insert(0, dmax)
        vale_freq.insert(0, freq_dmax)
        vale_sro_depl.insert(0, dmax)
        vale_freq.insert(0, 0.0)
        sro_prol_acce = NP.array(vale_sro_depl) * NP.array(vale_freq) ** (2) * kw["NORME"]
        para = f_in.para.copy()

        self.resu = t_fonction(vale_freq, sro_prol_acce, para)


class CalcFonction_CORR_ACCE(CalcFonctionOper):
    """CORR_ACCE"""

    def _run(self):
        """run CORR_ACCE"""
        f_in = self._lf[0]
        kw = self.kw
        para = f_in.para.copy()
        assert kw["METHODE"] in ("FILTRAGE", "POLYNOME")
        if kw["METHODE"] == "POLYNOME":
            # suppression de la tendance de l accelero
            fres = f_in.suppr_tend()
            # calcul de la vitesse
            fres = fres.trapeze(0.0)
            # calcul de la tendance de la vitesse : y = a1*x +a0
            fres = fres.suppr_tend()
            if self.kw["CORR_DEPL"] == "OUI":
                # suppression de la tendance deplacement
                # calcul du deplacement : integration
                fres = fres.trapeze(0.0)
                # calcul de la tendance du déplacement : y = a1*x +a0
                fres = fres.suppr_tend()
                # regeneration de la vitesse : derivation
                fres = fres.derive()
            # regeneration de l accelero : derivation
            self.resu = fres.derive()
            self.resu.para = para
        elif kw["METHODE"] == "FILTRAGE":
            dt = f_in.vale_x[1] - f_in.vale_x[0]
            acce_filtre = acce_filtre_CP(f_in.vale_y, dt, kw["FREQ_FILTRE"])
            self.resu = t_fonction(f_in.vale_x, acce_filtre, para)


class CalcFonction_DERIVE(CalcFonctionOper):
    """Derivation"""

    def _run(self):
        """DERIVE"""
        self.resu = self._lf[0].derive()


class CalcFonction_ENVELOPPE(CalcFonctionOper):
    """Return the envelop function"""

    def _run(self):
        """ENVELOPPE"""
        crit = self.kw["CRITERE"]
        if self.typres is Function2D:
            nap0 = self._lf[0]
            vale_para = nap0.vale_para
            para = nap0.para
            l_fonc_f = []
            for i in range(len(vale_para)):
                env = nap0.l_fonc[i]
                for nap in self._lf[1:]:
                    self.ctxt.f = nap.l_fonc[i].nom
                    env = enveloppe([env, nap.l_fonc[i]], crit)
                l_fonc_f.append(env)
            self.resu = t_nappe(vale_para, l_fonc_f, para)
        else:
            self.resu = enveloppe(self._lf, crit)


class CalcFonction_EXTRACTION(CalcFonctionOper):
    """Extract real/imaginary part"""

    def _build_data(self):
        dconv = {"REEL": "real", "IMAG": "imag", "MODULE": "modul", "PHASE": "phase"}
        arg = dconv[self.kw["PARTIE"]]
        self._build_list_fonc(arg=arg)

    def _run(self):
        """EXTRACTION"""
        self.resu = self._lf[0]


class CalcFonction_FFT(CalcFonctionOper):
    """Fast Fourier Transform"""

    def _build_data(self):
        """Read keywords to build the data"""
        opts = {}
        if self.typres is Function:
            opts["arg"] = "complex"
        self._build_list_fonc(**opts)

    def _run(self):
        """FFT"""
        kw = self.kw
        if self.typres is FunctionComplex:
            self.resu = self._lf[0].fft(kw["METHODE"])
        else:
            self.resu = self._lf[0].fft(kw["METHODE"], kw["SYME"])


class CalcFonction_INTERPOL_FFT(CalcFonctionOper):
    """Zero padding method"""

    def _build_data(self):
        """Read keywords to build the data"""
        opts = {}
        # if self.typres is fonction_sdaster:
        # opts['arg'] = 'complex'
        self._build_list_fonc(**opts)

    def _run(self):
        """INTERPOL_FFT"""
        kw = self.kw
        t0 = self._lf[0].vale_x[0]

        dt_init = self._lf[0].vale_x[1] - t0
        N_init = len(self._lf[0].vale_x)

        dt_cible = kw["PAS_INST"]
        if dt_init < dt_cible:
            UTMESS("F", "FONCT0_35")
        # nombre d'intervalles
        N_init -= 1
        N_sortie = int((N_init) * dt_init / dt_cible)

        if N_init * dt_init / dt_cible - N_sortie >= 0.5:
            N_sortie += 1
        # retour au nombre de valeurs
        N_sortie += 1

        # FFT
        ft = self._lf[0].fft("COMPLET")

        # suppression de la partie symetrique du signal
        N = len(ft.vale_x)
        valex = list(ft.vale_x[: N // 2 + 1])
        valey = list(ft.vale_y[: N // 2 + 1])

        # zero padding
        dfreq = (valex[1] - valex[0]).real
        last_freq = valex[-1]
        N_pad = N_sortie // 2 + 1 - N // 2 - 1
        for i in range(N_pad):
            freq = last_freq + (i + 1) * dfreq
            valex.append(freq)
            valey.append(0.0)
        ft.vale_x = NP.array(valex)
        ft.vale_y = NP.array(valey)

        # IFFT
        self.resu = ft.fft("COMPLET", "NON")
        self.resu.vale_x = self.resu.vale_x + t0

        # dt fin reel
        dt_fin = self.resu.vale_x[1] - self.resu.vale_x[0]

        # normalisation
        coef_norm = dt_init / dt_fin
        self.resu.vale_y = self.resu.vale_y * coef_norm

        ecart = abs(dt_fin - kw["PAS_INST"]) / kw["PAS_INST"]
        if ecart > kw["PRECISION"]:
            UTMESS("A", "FONCT0_51", valr=[dt_fin, kw["PAS_INST"], 100 * ecart])


class CalcFonction_FRACTILE(CalcFonctionOper):
    """Compute the fractile of functions"""

    def _run(self):
        """FRACTILE"""
        fract = self.kw["FRACT"]
        if self.typres is Function2D:
            nap0 = self._lf[0]
            vale_para = nap0.vale_para
            para = nap0.para
            l_fonc_f = []
            for i in range(len(vale_para)):
                self.ctxt.f = [nap.l_fonc[i].nom for nap in self._lf]
                lfr = fractile([nap.l_fonc[i] for nap in self._lf], fract)
                l_fonc_f.append(lfr)
            self.resu = t_nappe(vale_para, l_fonc_f, para)

        else:
            self.resu = fractile(self._lf, fract)


class CalcFonction_MOYENNE(CalcFonctionOper):
    """Compute the mean of functions"""

    def _run(self):
        """MOYENNE"""
        if self.typres is Function2D:
            nap0 = self._lf[0]
            vale_para = nap0.vale_para
            para = nap0.para
            l_fonc_f = []
            for i in range(len(vale_para)):
                self.ctxt.f = [nap.l_fonc[i].nom for nap in self._lf]
                lfr = moyenne([nap.l_fonc[i] for nap in self._lf])
                l_fonc_f.append(lfr)
            self.resu = t_nappe(vale_para, l_fonc_f, para)
        else:
            self.resu = moyenne(self._lf)


class CalcFonction_COHERENCE(CalcFonctionOper):
    """Compute the coherency function of two sets of signals"""

    def _build_data(self):
        """Read keywords to build the data"""
        self._build_list_fonc(mcsimp="NAPPE_1")

    def _run(self):
        """COHERENCE"""
        Mm = self.kw["NB_FREQ_LISS"]
        FREQ_COUP = self.kw["FREQ_COUP"]
        para = {
            "INTERPOL": ["LIN", "LIN"],
            "NOM_PARA": "FREQ",
            "PROL_DROITE": "CONSTANT",
            "PROL_GAUCHE": "EXCLU",
            "NOM_RESU": "ACCE",
        }
        nap1 = self._lf[0]
        assert nap1.para["NOM_PARA"] == "NUME_ORDRE"
        vale_para1 = nap1.vale_para
        nap2 = self.kw["NAPPE_2"]
        vale_para2, lfonc2 = nap2.Valeurs()
        assert len(vale_para1) == len(vale_para2), "NAPPE_1 and NAPPE_2 must have same length."
        assert set(vale_para2) == set(vale_para1), "Data lists are not ordered as pairs."

        acce1 = []
        acce2 = []
        for ii, fonc2 in enumerate(lfonc2):
            lt = nap1.l_fonc[ii].vale_x
            fonc1 = nap1.l_fonc[ii].vale_y
            assert len(lt) == len(
                fonc2[0]
            ), "Signals with same length required for NUME_ORDRE " + str(vale_para1[ii])
            assert (fonc2[0][1] - fonc2[0][0]) == (lt[1] - lt[0]), "same time steps required"
            if self.kw["OPTION"] == "DUREE_PHASE_FORTE":
                if ii == 0:
                    p1 = self.kw["BORNE_INF"]
                    p2 = self.kw["BORNE_SUP"]
                    N1, N2 = f_phase_forte(lt, fonc1, p1, p2)
                    UTMESS("I", "SEISME_79", valr=(lt[N1], lt[N2]))
                acce2.append(fonc2[1][N1:N2])
                acce1.append(fonc1[N1:N2])
            else:
                acce2.append(fonc2[1])
                acce1.append(fonc1)
        acce1 = NP.array(acce1)
        acce2 = NP.array(acce2)
        dt = lt[1] - lt[0]
        lfreq, fcohe = calc_cohefromdata(acce1, acce2, dt, Mm)
        N1 = NP.searchsorted(lfreq, 0.0)
        N2 = len(lfreq)
        if FREQ_COUP is not None:
            if lfreq[-1] > FREQ_COUP:
                N2 = NP.searchsorted(lfreq, FREQ_COUP)
        f_cohe = fcohe[N1:N2]
        l_freq = lfreq[N1:N2]
        self.resu = t_fonction(l_freq, f_cohe.real, para)


class CalcFonction_INTEGRE(CalcFonctionOper):
    """Integration"""

    def _run(self):
        """INTEGRE"""
        f_in = self._lf[0]
        kw = self.kw
        assert kw["METHODE"] in ("TRAPEZE", "SIMPSON")
        if kw["METHODE"] == "TRAPEZE":
            self.resu = f_in.trapeze(kw["COEF"])
        elif kw["METHODE"] == "SIMPSON":
            self.resu = f_in.simpson(kw["COEF"])


class CalcFonction_INVERSE(CalcFonctionOper):
    """Reverse"""

    def _run(self):
        """INVERSE"""
        self.resu = self._lf[0].inverse()


class CalcFonction_MULT(CalcFonctionOper):
    """Multiply the given functions."""

    def _build_data(self):
        """Read keywords to build the data"""
        opts = {}
        if self.typres is FunctionComplex:
            opts["arg"] = "complex"
        self._build_list_fonc(**opts)

    def _run(self):
        """MULT"""
        self.resu = 1.0
        for item in self._lf:
            self.ctxt.f = item.nom
            self.resu = item * self.resu
        # take the parameters of the first function
        self.resu.para = self._lf[0].para.copy()
        self._use_list_para()


class CalcFonction_PUISSANCE(CalcFonctionOper):
    """Compute f^n"""

    def _run(self):
        """PUISSANCE"""
        self.resu = self._lf[0]
        for i in range(self.kw["EXPOSANT"] - 1):
            self.resu = self.resu * self._lf[0]


class CalcFonction_SPEC_OSCI(CalcFonctionOper):
    """SPEC_OSCI"""

    def _build_data(self):
        """Read keywords to build the data"""
        CalcFonctionOper._build_list_fonc(self)
        kw = self.kw
        self._dat = {}
        # amor
        if kw["AMOR_REDUIT"] is None:
            l_amor = [0.02, 0.05, 0.1]
            UTMESS("I", "FONCT0_31", valr=l_amor)
        else:
            l_amor = force_list(kw["AMOR_REDUIT"])
        eps = 1.0e-6
        for amor in l_amor:
            if amor > (1 - eps):
                UTMESS("S", "FONCT0_36")
        self._dat["AMOR"] = l_amor
        # freq
        if kw["LIST_FREQ"] is not None:
            l_freq = kw["LIST_FREQ"].getValues()
        elif kw["FREQ"] is not None:
            l_freq = force_list(kw["FREQ"])
        else:
            l_freq = []
            for i in range(56):
                l_freq.append(0.2 + 0.050 * i)
            for i in range(8):
                l_freq.append(3.0 + 0.075 * i)
            for i in range(14):
                l_freq.append(3.6 + 0.100 * i)
            for i in range(24):
                l_freq.append(5.0 + 0.125 * i)
            for i in range(28):
                l_freq.append(8.0 + 0.250 * i)
            for i in range(6):
                l_freq.append(15.0 + 0.500 * i)
            for i in range(4):
                l_freq.append(18.0 + 1.0 * i)
            for i in range(10):
                l_freq.append(22.0 + 1.500 * i)
            texte = []
            for i in range(len(l_freq) // 5):
                texte.append(" %f %f %f %f %f" % tuple(l_freq[i * 5 : i * 5 + 5]))
            UTMESS("I", "FONCT0_32", vali=len(l_freq), valk=os.linesep.join(texte))
        if min(l_freq) < 1.0e-10:
            UTMESS("F", "FONCT0_43")
        self._dat["FREQ"] = l_freq
        # check
        if abs(kw["NORME"]) < 1.0e-10:
            UTMESS("S", "FONCT0_33")
        if kw["NATURE_FONC"] == "DSP":
            ASSERT(kw["METHODE"] == "RICE")

    def _run(self):
        """SPEC_OSCI"""
        f_in = self._lf[0]
        l_freq, l_amor = self._dat["FREQ"], self._dat["AMOR"]
        kw = self.kw
        l_fonc_f = []
        # construction de la nappe ou de la fonction
        vale_para = l_amor
        para = {
            "INTERPOL": ["LIN", "LOG"],
            "NOM_PARA_FONC": "FREQ",
            "NOM_PARA": "AMOR",
            "PROL_DROITE": "EXCLU",
            "PROL_GAUCHE": "EXCLU",
            "NOM_RESU": kw["NATURE"],
        }
        para_fonc = {
            "INTERPOL": ["LOG", "LOG"],
            "NOM_PARA": "FREQ",
            "PROL_DROITE": "CONSTANT",
            "PROL_GAUCHE": "EXCLU",
            "NOM_RESU": kw["NATURE"],
        }
        if kw["NATURE"] == "DEPL":
            ideb = 0
        elif kw["NATURE"] == "VITE":
            ideb = 1
        else:
            ASSERT(kw["NATURE"] == "ACCE")
            ideb = 2
        if kw["METHODE"] == "RICE":
            # appel à DSP2SRO
            ASSERT(kw["NATURE_FONC"] == "DSP")
            deuxpi = 2.0 * math.pi
            f_dsp = t_fonction(f_in.vale_x * deuxpi, f_in.vale_y / deuxpi, f_in.para)
            for iamor in l_amor:
                spectr = DSP2SRO(f_dsp, iamor, kw["DUREE"], l_freq, ideb)
                vale_y = spectr.vale_y / kw["NORME"]
                l_fonc_f.append(t_fonction(l_freq, vale_y, para_fonc))
        elif kw["METHODE"] == "NIGAM":
            # appel à SPEC_OSCI
            ASSERT(kw["NATURE_FONC"] == "ACCE")
            spectr = aster_fonctions.SPEC_OSCI(f_in.vale_x, f_in.vale_y, l_freq, l_amor)
            for iamor in range(len(l_amor)):
                vale_y = spectr[iamor, ideb, :] / kw["NORME"]
                l_fonc_f.append(t_fonction(l_freq, vale_y, para_fonc))
        elif kw["METHODE"] == "HARMO":
            # appel à ACCE2DSP
            ASSERT(kw["NATURE_FONC"] == "ACCE")
            for iamor in l_amor:
                spectr = ACCE2SRO(f_in, iamor, l_freq, ideb)
                vale_y = spectr.vale_y / kw["NORME"]
                l_fonc_f.append(t_fonction(l_freq, vale_y, para_fonc))
        if self.typres == Function2D:
            self.resu = t_nappe(vale_para, l_fonc_f, para)
        else:
            self.resu = l_fonc_f[0]


class CalcFonction_DSP(CalcFonctionOper):
    """DSP"""

    def _run(self):
        """DSP"""
        kw = self.kw
        f_in = self._lf[0]
        vale_freq = f_in.vale_x
        vale_sro = f_in.vale_y
        f_min = f_in.vale_x[0]
        f_in = t_fonction(NP.insert(vale_freq, 0, 0.0), NP.insert(vale_sro, 0, 0.0), para=f_in.para)
        deuxpi = 2.0 * math.pi
        freq_coup = kw["FREQ_COUP"]
        SRO_args = {
            "DUREE_PHASE_FORTE": kw["DUREE"],
            "FREQ_COUP": freq_coup,
            "NORME": kw["NORME"],
            "AMORT": kw["AMOR_REDUIT"],
            "FMIN": f_min,
            "FONC_SPEC": f_in,
            "NITER": kw["NB_ITER"],
            "FREQ_FILTRE_ZPA": kw["FREQ_FILTRE_ZPA"],
            "NB_FREQ_LISS": kw["NB_FREQ_LISS"],
        }
        if kw["FREQ_PAS"] is not None:
            SRO_args["PAS"] = kw["FREQ_PAS"]
        elif kw["LIST_FREQ"] is not None:
            l_freq = kw["LIST_FREQ"].getValues()
            if l_freq[0] <= 0.0:
                UTMESS("F", "FONCT0_43")
            SRO_args["LIST_FREQ"] = l_freq
            SRO_args["PAS"] = None
        f_dsp, f_sro_ref = SRO2DSP(**SRO_args)
        self.resu = t_fonction(f_dsp.vale_x / deuxpi, f_dsp.vale_y * deuxpi, para=f_in.para)


class CalcFonction_LISS_ENVELOP(CalcFonctionOper):
    """LISS_ENVELOP"""

    def _build_data(self):
        """Read keywords to build the data"""
        # CalcFonctionOper._build_list_fonc(self)
        kw = self.kw
        if self.kw["NAPPE"] is not None:
            self._build_list_fonc(mcsimp="NAPPE")
        elif self.kw["FONCTION"] is not None:
            self._build_list_fonc(mcsimp="FONCTION")
        elif self.kw["TABLE"] is not None:
            lf_in = self._get_mcsimp("TABLE")
            para_fonc = {
                "PROL_DROITE": "EXCLU",
                "INTERPOL": ["LIN", "LIN"],
                "PROL_GAUCHE": "EXCLU",
                "NOM_RESU": "TOUTRESU",
                "NOM_PARA": "FREQ",
            }
            para_napp = {
                "PROL_DROITE": "EXCLU",
                "INTERPOL": ["LIN", "LIN"],
                "PROL_GAUCHE": "EXCLU",
                "NOM_RESU": "TOUTRESU",
                "NOM_PARA": "AMOR",
                "NOM_PARA_FONC": "FREQ",
            }

            # conversion des tables en nappes
            l_nappe = []
            for tab in lf_in:
                nom_para = tab.get_nom_para()
                if kw["LIST_AMOR"] is not None:
                    amor = kw["LIST_AMOR"]
                else:
                    amor = list(range(1, len(nom_para)))
                # error
                assert "FREQ" in nom_para, nom_para
                nom_para.remove("FREQ")
                dico = tab.EXTR_TABLE().values()
                l_fonc_f = []
                for para in nom_para:
                    freq = dico["FREQ"]
                    vale = dico[para]
                    l_fonc_f.append(t_fonction(freq, vale, para_fonc))
                l_nappe.append(t_nappe(amor, l_fonc_f, para_napp))
            self._lf = l_nappe

    def _run(self):
        """LISS_ENVELOP"""
        kw = self.kw
        l_sp_nappe = []
        # formatage selon les donnes d'entrees
        if kw["FONCTION"] is not None:
            f_in = self._lf[0]
            if kw["LIST_AMOR"] is not None:
                amor = kw["LIST_AMOR"][0]
            else:
                amor = 0.0
            sp_nappe = LISS.nappe(
                listFreq=f_in.vale_x, listeTable=[f_in.vale_y], listAmor=[amor], entete=""
            )
            para_fonc = f_in.para
            para = f_in.para.copy()
            para["NOM_PARA"] = "AMOR"
            para["NOM_PARA_FONC"] = para_fonc["NOM_PARA"]
            l_sp_nappe = [sp_nappe]
        elif kw["NAPPE"] is not None or kw["TABLE"] is not None:
            for i_nappe in range(len(self._lf)):
                f_in = self._lf[i_nappe]
                sp_nappe = LISS.nappe(
                    listFreq=f_in.l_fonc[0].vale_x,
                    listeTable=[f.vale_y for f in f_in.l_fonc],
                    listAmor=f_in.vale_para,
                    entete="",
                )
                # verification que les nappes ont les memes amortissements
                if i_nappe == 0:
                    l_amor = f_in.vale_para
                    l_amor.sort()
                    erreur_amor = 0
                else:
                    if len(l_amor) == len(f_in.vale_para):
                        l_amor2 = f_in.vale_para.copy()
                        l_amor2.sort()
                        d_amor = l_amor - l_amor2
                        if max(abs(d_amor)) > 1e-6:
                            erreur_amor = 1
                    else:
                        erreur_amor = 1
                if erreur_amor:
                    UTMESS("F", "FONCT0_74")
                l_sp_nappe.append(sp_nappe)
            para_fonc = f_in.l_fonc[0].para
            para = f_in.para

        if kw["OPTION"] == "CONCEPTION":
            sp_lisse = LISS.liss_enveloppe(
                l_sp_nappe,
                option=kw["OPTION"],
                fmin=kw["FREQ_MIN"],
                fmax=kw["FREQ_MAX"],
                l_freq=list(kw["LIST_FREQ"] or []),
                nb_pts=kw["NB_FREQ_LISS"],
                zpa=kw["ZPA"],
                precision=1e-5,
                critere="RELATIF",
            )
        else:
            sp_lisse = LISS.liss_enveloppe(
                l_sp_nappe,
                option=kw["OPTION"],
                coef_elarg=kw["ELARG"],
                fmin=kw["FREQ_MIN"],
                fmax=kw["FREQ_MAX"],
                l_freq=list(kw["LIST_FREQ"] or []),
                nb_pts=kw["NB_FREQ_LISS"],
                zpa=kw["ZPA"],
                precision=1e-5,
                critere="RELATIF",
            )

        l_fonc_f = []
        for spec in sp_lisse.listSpec:
            l_fonc_f.append(t_fonction(sp_lisse.listFreq, spec.dataVal, para_fonc))
        self.resu = t_nappe(sp_lisse.listAmor, l_fonc_f, para)


class CalcFonction_REGR_POLYNOMIALE(CalcFonctionOper):
    """Polynomial regression"""

    def _run(self):
        """REGR_POLYNOMIALE"""
        f_in = self._lf[0]
        deg = self.kw["DEGRE"]
        coef = NP.polyfit(f_in.vale_x, f_in.vale_y, deg)
        if coef is None:
            raise AsterError("La régression polynomiale n'a pas convergé.")
        # interpolation sur une liste d'abscisses
        absc = f_in.vale_x
        if self.args["LIST_PARA"] is not None:
            absc = self.args["LIST_PARA"].getValues()
        vale = NP.polyval(coef, absc)
        # paramètres
        para = f_in.para.copy()
        para["INTERPOL"] = ["LIN", "LIN"]
        self.resu = t_fonction(absc, vale, para)
        coef_as_str = os.linesep.join(["   a[%d] = %f" % (i, ci) for i, ci in enumerate(coef)])
        UTMESS("I", "FONCT0_57", coef_as_str)


class CalcFonction_INTEGRE_FREQ(CalcFonctionOper):
    """INTEGRE_FREQ"""

    def _run(self):
        """INTEGRE_FREQ"""
        kw = self.kw
        f_in = self._lf[0]
        para = f_in.para.copy()
        para.update(_F(PROL_GAUCHE="CONSTANT", PROL_DROITE="CONSTANT"))
        vale_t = f_in.vale_x
        vale_s = f_in.vale_y
        dt = vale_t[1] - vale_t[0]
        nbdt = len(vale_t)
        tfin = vale_t[nbdt - 1]
        fmax = 0.5 / dt  # facteur 1/2 car FFT calculee avec SYME='NON'
        df = 2.0 * fmax / nbdt
        freq_coup = kw["FREQ_COUP"]
        freq_filt = kw["FREQ_FILTRE"]

        if freq_filt > 0.0:
            # CORR_ACCE / METHODE="FILTRAGE"
            vale_s = acce_filtre_CP(vale_s, dt, freq_filt)

        acc0 = t_fonction(vale_t, vale_s, para)
        xff0 = acc0.fft(methode="COMPLET")

        lfreq = NP.arange(0.0, fmax, df)

        para_filt = para.copy()
        para_filt.update(_F(INTERPOL=["LIN", "LIN"]))
        filtre = t_fonction_c([0.0, df, freq_coup, freq_coup + df], [0.0, 1.0, 1.0, 0.0], para_filt)
        xf1 = (xff0 * filtre).evalfonc(lfreq)
        if kw["NIVEAU"] == 2:
            vale = -1.0 / (2.0 * math.pi * lfreq[1:]) ** 2 * xf1.vale_y[1:]
        else:
            vale = 1.0 / (2.0 * math.pi * lfreq[1:]) * xf1.vale_y[1:]
            vale = -1j * vale

        xf0 = t_fonction_c(lfreq, NP.concatenate(([0.0], vale)), para)
        xf0.para["NOM_PARA"] = "FREQ"
        depl = xf0.fft(methode="COMPLET", syme="NON")
        dep0 = -depl(0.0)

        linst = NP.arange(0.0, tfin + dt, dt)
        result = (depl + dep0).evalfonc(linst)
        result.para = f_in.para.copy()
        self.resu = result


class CalcFonction_DERIVE_FREQ(CalcFonctionOper):
    """DERIVE_FREQ"""

    def _run(self):
        """DERIVE_FREQ"""
        kw = self.kw
        f_in = self._lf[0]
        para = f_in.para.copy()
        para.update(_F(PROL_GAUCHE="CONSTANT", PROL_DROITE="CONSTANT"))
        vale_t = f_in.vale_x
        dt = vale_t[1] - vale_t[0]
        nbdt = len(vale_t)
        tfin = vale_t[nbdt - 1]
        fmax = 0.5 / dt  # facteur 1/2 car FFT calculee avec SYME='NON'
        df = 2.0 * fmax / nbdt
        freq_coup = kw["FREQ_COUP"]

        xff0 = f_in.fft(methode="COMPLET")

        lfreq = NP.arange(0.0, fmax, df)

        para_filt = para.copy()
        para_filt.update(_F(INTERPOL=["LIN", "LIN"]))
        filtre = t_fonction_c([0.0, df, freq_coup, freq_coup + df], [0.0, 1.0, 1.0, 0.0], para_filt)
        xf1 = (xff0 * filtre).evalfonc(lfreq)
        if kw["NIVEAU"] == 2:
            vale = -1.0 * (2.0 * math.pi * lfreq[1:]) ** 2 * xf1.vale_y[1:]
        else:
            vale = 2.0 * math.pi * lfreq[1:] * xf1.vale_y[1:]
            vale = 1j * vale

        xf0 = t_fonction_c(lfreq, NP.concatenate(([0.0], vale)), para)
        xf0.para["NOM_PARA"] = "FREQ"
        depl = xf0.fft(methode="COMPLET", syme="NON")

        linst = NP.arange(0.0, tfin + dt, dt)
        result = depl.evalfonc(linst)
        result.para = f_in.para.copy()
        self.resu = result


class Context:
    """Permet de stocker des éléments de contexte pour aider au
    diagnostic lors de l'émission de message.
    usage :
       context = Context()
       context.f = 'nomfon'
       print context.f
    """

    def __init__(self):
        self.__nomf = None

    def get_val(self):
        """Retourne le texte formatté."""
        nomf = self.__nomf
        if type(nomf) not in (list, tuple):
            nomf = [nomf]
        pluriel = ""
        if len(nomf) > 1:
            pluriel = "s"
        res = """Fonction%(s)s concernée%(s)s : %(nomf)s""" % {
            "s": pluriel,
            "nomf": ", ".join(nomf),
        }
        return res

    def set_val(self, value):
        """Set function"""
        self.__nomf = value

    def del_val(self):
        """Remove value"""
        del self.__nomf

    f = property(get_val, set_val, del_val, "")
