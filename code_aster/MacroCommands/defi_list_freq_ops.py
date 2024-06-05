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

from math import fmod, sqrt

from ..Cata.Syntax import _F
from ..Cata.SyntaxUtils import remove_none
from ..CodeCommands import DEFI_LIST_REEL
from ..Messages import UTMESS
from ..Utilities import no_new_attributes


class FrequencyDefinition:
    """Definition of a list of frequencies.

    Arguments:
        keywords (dict): Keywords passed to the command.
    """

    list_args = refine_args = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, keywords) -> None:
        args = _F(keywords)
        remove_none(args)
        self.refine_args = args.pop("RAFFINEMENT", [None])[0]
        self.list_args = args

    def create_initial_list(self):
        """Create the list of frequencies

        Returns:
            list[float]: List of frequencies as float values.
        """
        listr8 = DEFI_LIST_REEL(**self.list_args)
        return listr8.getValues()

    def refine(self, list_init):
        """Refine the given list."""
        # 2. Récuperation des données liées au raffinement
        if self.refine_args.get("CRITERE") in ("RELATIF", "ABSOLU"):
            dispersion = self.refine_args.get("DISPERSION")

        haveAmor = False
        if self.refine_args.get("CRITERE") == "LARGEUR_3DB":
            haveAmor = True
            if self.refine_args.get("AMOR_REDUIT") is not None:
                l_amor = list(self.refine_args.get("AMOR_REDUIT"))
            else:
                l_amor = self.refine_args.get("LIST_AMOR").getValues()

        dfMin = self.refine_args.get("PAS_MINI")
        nbPtsRaf = self.refine_args.get("NB_POINTS")
        l_freq = self.refine_args.get("LIST_RAFFINE")

        # Si le nombre d'amortissements donnés est inférieur au nombre de fréquences données
        # dans LIST_RAFFINE, les amortissements des fréquences supplémentaires sont
        # pris égaux au dernier amortissement de la liste.
        if haveAmor:
            if len(l_amor) < len(l_freq):
                for i in range(len(l_freq) - len(l_amor)):
                    l_amor.append(l_amor[-1])

        # 3. Elimination des modes multiples dans la liste des fréquences l_freq
        l_freq_clean = [l_freq[0]]
        if haveAmor:
            lamor_clean = [l_amor[0]]

        for i in range(1, len(l_freq)):
            if (l_freq[i] - l_freq[i - 1]) > dfMin:
                l_freq_clean.append(l_freq[i])
                if haveAmor:
                    lamor_clean.append(l_amor[i])

        l_freq = l_freq_clean[:]
        if haveAmor:
            l_amor = lamor_clean[:]

        # 4. Raffinement de la liste des fréquences

        # On stocke toutes les fréquences des intervalles raffinés, dans la liste l_raf.
        # Celle-ci contient automatiquement toutes les valeurs de la liste de base.
        l_raf = []

        for i in range(0, len(l_freq)):
            if haveAmor:
                if l_amor[i] > 1e-12:
                    # on décale la fréquence l_freq[i] pour se centrer sur la
                    # résonance d'amplitude
                    l_freq[i] = l_freq[i] * sqrt(1 - 2 * l_amor[i] ** 2)
                    # largeur de l'intervalle i à raffiner
                    df = 2 * l_amor[i] * l_freq[i]
                else:
                    df = 0.01 * l_freq[i]
            elif self.refine_args.get("CRITERE") == "RELATIF":
                df = dispersion * l_freq[i]
            elif self.refine_args.get("CRITERE") == "ABSOLU":
                df = dispersion

            ltemp = [l_freq[i] - df / 2.0 + j * df / (nbPtsRaf - 1) for j in range(0, nbPtsRaf)]

            if i > 1:
                # on vérifie s'il y a un recouvrement d'intervalle
                if ltemp[0] < max(l_raf):
                    freq_i = l_freq[i]
                    freq_imoins1 = l_freq[i - 1]
                    if haveAmor:
                        freq_i = l_freq[i] / sqrt(1 - 2 * l_amor[i] ** 2)
                        freq_imoins1 = l_freq[i - 1] / sqrt(1 - 2 * l_amor[i - 1] ** 2)
                    UTMESS("I", "DYNAMIQUE_26", valr=(freq_imoins1, freq_i))

            if haveAmor:
                # si le nombre de points à ajouter est pair :
                if not fmod(nbPtsRaf, 2):
                    # on crée un point au milieu (résonance d'amplitude)
                    ltemp.append(l_freq[i])
                    # résonance de phase
                    ltemp.append(l_freq[i] / sqrt(1 - 2 * l_amor[i] ** 2))

            l_raf += ltemp

        l_raf += list_init
        l_raf.sort()

        # 5. Elimination des fréquences trop proches (écart < dfMin)

        l_freq_ok = l_freq + list_init

        i = 1
        while i < len(l_raf):
            if (l_raf[i] - l_raf[i - 1]) < dfMin:
                if l_raf[i - 1] in l_freq_ok and l_raf[i] not in l_freq_ok:
                    l_raf.remove(l_raf[i])
                    i -= 1
                elif l_raf[i] in l_freq_ok and l_raf[i - 1] not in l_freq_ok:
                    l_raf.remove(l_raf[i - 1])
                    i -= 1
                elif l_raf[i] in l_freq_ok and l_raf[i - 1] in l_freq_ok:
                    UTMESS("I", "DYNAMIQUE_27", valr=(l_raf[i - 1], l_raf[i]))
                else:
                    l_raf.remove(l_raf[i])
                    i -= 1
            i += 1
        return l_raf


def defi_list_freq_ops(self, **args):
    """Definition of a list of frequencies."""
    builder = FrequencyDefinition(args)

    list_freq = builder.create_initial_list()
    refined_list = builder.refine(list_freq)

    return DEFI_LIST_REEL(VALE=refined_list)
