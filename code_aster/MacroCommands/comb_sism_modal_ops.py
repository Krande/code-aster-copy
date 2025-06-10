# coding: utf-8

# Copyright (C) 1991 - 2025 EDF R&D                www.code-aster.org
#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.

from itertools import count, product
from math import sqrt

import numpy as np
from libaster import setFortranLoggingLevel

from ..Messages import UTMESS
from ..Objects import MultipleElasticResult
from ..Supervis.ExecuteCommand import _get_object_repr


def get_nodes(mesh, GROUP_NO=None):
    """Get the number of nodes.

    Arguments:
        mesh: the mesh
        GROUP_NO (list[str]): list name of groups of nodes

    Returns:
        tuple: Tuple containing the list of number of nodes and the list of
        the nodes name for each group.
    """

    nodes_num = []
    nodes_name = []
    if GROUP_NO is not None:
        nodes_num = mesh.getNodes(GROUP_NO)
        nodes_name = [mesh.getNodeName(i) for i in nodes_num]
    return nodes_num, nodes_name


def get_spectres_mono_appui(spectre_in):
    """Get the input for SPECTRE for mono-appui

    Arguments:
        spectre_in: input for key_factor SPECTRE

    Returns:
        Spectres: list of [directions, nappes, coefficients, cor_freqs, natures]
        directions: list of OX, OY, OZ
        nappes: list of aster objects as nappe
        coefficients: list of coefficient for each direction
        cor_freqs: list of option to correct freq as "OUI" or "NON"
        natures: list of nature of spectra
    """
    coefficients = []
    directions = []
    nappes = []
    cor_freqs = []
    natures = []
    # for each key_factor
    for i in range(len(spectre_in)):
        # extraction of directions
        list_axe = spectre_in[i].get("LIST_AXE")
        for axe in list_axe:
            if axe == "X":
                dir = "OX"
            elif axe == "Y":
                dir = "OY"
            elif axe == "Z":
                dir = "OZ"
            # save
            directions.append(dir)
            # extraction of parameters for associated direction
            # scale factor
            coefficients.append(spectre_in[i].get("ECHELLE"))
            # aster object: nappe
            nappes.append(spectre_in[i].get("SPEC_OSCI"))
            # correction of frequency
            cor_freqs.append(spectre_in[i].get("CORR_FREQ"))
            # nature of spectre
            natures.append(spectre_in[i].get("NATURE"))
    if i > 3:
        raise Exception(
            "En cas MONO_APPUI, Il faut renseigner au maximun " + "trois axes globaux X, Y, Z"
        )
    # save
    spectres = [directions, nappes, coefficients, cor_freqs, natures]
    # return
    return spectres


def get_spectres_mult_appui(spectre_in, mesh):
    """Get the input for SPECTRE in case of mult-appui

    Arguments:
        spectre_in: input for key_factor SPECTRE

    Returns:
        Spectres: list of [directions, nappes, coefficients, cor_freqs, nom_appui]
        directions: list of X, Y, Z
        nappes: list of aster objects as nappe
        coefficients: list of coefficient for each direction
        cor_freqs: list of option to correct freq as "OUI" or "NON"
        natures: list of nature of spectra
        nom_appui: name of support associated to the spectra
    """
    # extraction all nom_appui in key_factor SPECTRE
    nom_appuis_all = []
    for i in range(len(spectre_in)):
        nom_appuis_all.append(spectre_in[i].get("NOM_APPUI"))
    # For each support: get all dirs, coef, nappes, cor_freqs et natures
    l_nom_appui = list(set(nom_appuis_all))
    # get info spectre for each support
    spectres = []
    for nom_appui in l_nom_appui:
        coefficients = []
        directions = []
        nappes = []
        cor_freqs = []
        natures = []
        # serach for all indexes of key_factor corresponding the nom_appui
        indexes = [i for i, x in enumerate(nom_appuis_all) if x == nom_appui]
        # extraction data from coef, dir, nappe, cor_freq, nature
        for i in range(len(indexes)):
            # -- extraction
            list_axe = spectre_in[indexes[i]].get("LIST_AXE")
            for axe in list_axe:
                if axe == "X":
                    dir = "OX"
                elif axe == "Y":
                    dir = "OY"
                elif axe == "Z":
                    dir = "OZ"
                # save
                directions.append(dir)
                # get scale factor
                coefficients.append(spectre_in[indexes[i]].get("ECHELLE"))
                # get nappe object
                nappes.append(spectre_in[indexes[i]].get("SPEC_OSCI"))
                # get corr_freq option
                cor_freqs.append(spectre_in[indexes[i]].get("CORR_FREQ"))
                # get nature of spectrum
                natures.append(spectre_in[indexes[i]].get("NATURE"))
        # check if more than 3 key_factor is used for nom_appui
        if i > 3:
            raise Exception(
                "En cas MULT_APPUI, Il faut renseigner au maximun "
                + "trois axes globaux X, Y, Z par APPUI"
            )
        # save
        spectres.append([directions, nappes, coefficients, cor_freqs, natures, nom_appui])
    # return
    return spectres


def get_appuis(appuis_in, mesh):
    """Get the input for APPUIS in case of Multi Appui

    Arguments:
        appuis_in: input for key_factor APPUIS

    Returns:
        nom_appuis_all: list of all nom_appui
        l_nodes_num_all : list of all node numbers for each nom_appui
        l_nodes_name_all : list of all node names for each nom_appui
    """
    # preparing of saving listes
    nom_appuis_all = []
    l_nodes_num_all = []
    l_nodes_name_all = []
    l_group_no_all = []
    # get information on definition of each APPUI
    for i in range(len(appuis_in)):
        # get nom appui
        nom_appuis_all.append(appuis_in[i].get("NOM"))
        # get group_no to build APPUI
        group_no = appuis_in[i].get("GROUP_NO")
        l_group_no_all.append(group_no)
        # get information on nodes to build APPUI
        l_nodes_num, l_nodes_name = get_nodes(mesh, GROUP_NO=group_no)
        # save list of all node numbers
        l_nodes_num_all.append(l_nodes_num)
        # save list of all node names
        l_nodes_name_all.append(l_nodes_name)
    # traitement
    # check if one node belonging to two APPUIS
    all_node_num = []
    for i_appui in range(len(l_nodes_num_all)):
        all_node_num += l_nodes_num_all[i_appui]
    if len(all_node_num) != len(set(all_node_num)):
        raise Exception("Il y a au moins un noeud appartenant a deux appuis")
    # return
    return nom_appuis_all, l_nodes_num_all, l_nodes_name_all, l_group_no_all


def get_group_appuis(spectres, group_appui_correle=None):
    """Get information about group_appuis

    Arguments:
        spectres: input of operande SPECTRE
        group_appui_correle: input for operande GROUP_APPUI_CORRELE

    Returns:
        group_appui_all
        nom_group_appui_all
    """
    # preparing for group_appui_correle
    group_appui_all = []
    nom_group_appui_all = []
    # add group_appui by users
    if group_appui_correle:
        # get information of input for GROUP_APPUI_CORRELE
        for i_group_appui in range(len(group_appui_correle)):
            # add APPUIS to group_appui
            if group_appui_correle[i_group_appui].get("LIST_APPUI") is not None:
                group_appui_all.append(group_appui_correle[i_group_appui].get("LIST_APPUI"))
                nom_group_appui_all.append(group_appui_correle[i_group_appui].get("NOM"))
            # add all APPUIS to group_appui
            elif group_appui_correle[i_group_appui].get("TOUT") == "OUI":
                group_appui_all.append([spectres[i_appui][5] for i_appui in range(len(spectres))])
                nom_group_appui_all.append(
                    group_appui_correle[i_group_appui].get("NOM_GROUP_APPUI")
                )
        # all APPUIS not mentionned will be regrouped in separed groupe named as nom_appui
        for i_appui in range(len(spectres)):
            if all(spectres[i_appui][5] not in i for i in group_appui_all):
                group_appui_all.append(spectres[i_appui][5])
                nom_group_appui_all.append(str(spectres[i_appui][5]))
    else:  # comb_mult_appui_corr is not mentionned --> all appuis are in the same group_appui_correle
        for i_appui in range(len(spectres)):
            group_appui_all.append(spectres[i_appui][5])
            nom_group_appui_all.append(str(spectres[i_appui][5]))
    # return
    return group_appui_all, nom_group_appui_all


def get_amor_reduit(mode_meca, l_nume_ordre, amor_reduit, list_amor, amor_gene):
    """Compute array of modal damping for each mode

    Arguments:
        mode_meca: modal basis (resu as mode_meca aster SD)
        l_nume_ordre: list of all numeros of orders (nume_ordre)
        amor_reduit: a damping coefficient mentionned by users
        list_amor : list of damping coefficients by users
        amor_gene : generalized damping matrix by users

    Returns:
        amor_reduit: a scalar or a liste of damping coefficients
    """
    # List of parameters in mode_meca
    list_para = mode_meca.LIST_PARA()
    # get values of damping coefficient
    if list_amor is not None:
        # get values of damping in list_amor
        amor_reduit = list_amor.getValues()
    elif amor_gene is not None:
        # computing damping coefficient from generalized damping matrix
        l_amors = amor_gene.getUpperValues()
        amor_reduit = []
        for amor, frq, ordre, mass in zip(
            l_amors, list_para["FREQ"], list_para["NUME_ORDRE"], list_para["MASS_GENE"]
        ):
            if ordre in l_nume_ordre:
                amor_reduit.append(amor / (4 * np.pi * frq * mass))
    # return
    return amor_reduit


def filter_ordre_freq(list_para, increment):
    """Filter the order of mode from mode_meca

    Arguments:
        list_para    : list_para from modal basis
        increment    : increment keyword

    Returns:
        freqs: ordered list of filtered frequencies 
        nume_ordres: ordered list of filtered number order
        nume_modes: ordered list of filtered number mode
        gene_masses: ordered list of masse gene
    """
    l_nume_ordre = list_para["NUME_ORDRE"]
    l_freq = list_para["FREQ"]
    l_nume_mode = list_para["NUME_MODE"]
    l_gene_masses = list_para["MASS_GENE"]

    tout_ordre = increment.get("TOUT_ORDRE")
    list_ordre = increment.get("LIST_ORDRE")
    nume_mode = increment.get("NUME_MODE")
    nume_ordre = increment.get("NUME_ORDRE")
    list_freq = increment.get("LIST_FREQ")
    freq = increment.get("FREQ")
    precision = increment.get("PRECISION")
    critere = increment.get("CRITERE")

    # Gerer le cas où TOUT_ORDRE est la valeur par default si tout est à None
    if all(key is None for key in (tout_ordre, nume_ordre, freq, nume_mode, list_freq, list_ordre)):
        tout_ordre = "OUI"
    # all modes in modal basis are selected
    if tout_ordre == "OUI":
        l_nume_ordre_filtered = l_nume_ordre
        l_freq_filtered = l_freq
        l_nume_mode_filtered = l_nume_mode
    else:
        # selecting mode by nume_mode
        nume_mode = nume_mode if nume_mode is not None else []
        # selecting mode by nume_ordre
        nume_ordre = nume_ordre if nume_ordre is not None else []
        # selecting mode by list of nume_ordre
        if list_ordre is not None:
            nume_ordre = list_ordre.getValues()
        # selecting mode by list of frequencies
        elif list_freq is not None:
            freq = list_freq.getValues()
        # filtering frequencies
        if freq is not None:
            freqs = np.array(freq)
            if critere == "RELATIF":
                borne_inf, borne_sup = freqs * (1 - precision), freqs * (1 + precision)
            else:
                borne_inf, borne_sup = freqs - precision, freqs + precision
        else:
            borne_inf, borne_sup = -1, -1
        # filtered order modes
        l_nume_ordre_filtered, l_freq_filtered, l_nume_mode_filtered = [], [], []
        # Filtration des ordres
        for i_ordre, f, i_mode in zip(l_nume_ordre, l_freq, l_nume_mode):
            is_in_interval = np.any((borne_inf <= f) & (f <= borne_sup))
            if i_ordre in nume_ordre or i_mode in nume_mode or is_in_interval:
                l_nume_ordre_filtered.append(i_ordre)
                l_freq_filtered.append(f)
                l_nume_mode_filtered.append(i_mode)
    
    # Sort frequencies in increasing order
    _, _, argsort_freq = zip(*sorted(zip(l_freq_filtered, l_nume_ordre_filtered, count())))
    argsort_freq = list(argsort_freq)
    freqs = np.array(l_freq_filtered)[argsort_freq]
    # Sorted number of order
    nume_ordres = np.array(l_nume_ordre_filtered)[argsort_freq]
    # Sorted number of modes
    nume_modes = np.array(l_nume_mode_filtered)[argsort_freq]
    # Sorted generalized masses
    gene_masses = np.array(l_gene_masses)[argsort_freq]

    return freqs, nume_ordres, nume_modes, gene_masses



class CombModalResponse:

    def __init__(self, comb_mode, type_analyse, amors, freqs):
        """Combines modals responses

        Args:
            comb_mode: input for key_factor COMB_MODE
            type_analyse: input for key TYPE_ANALYSE (MONO_APPUI or MULT_APPUI)
            amors: list of all damping coefficients
            freqs: list of all frequencies
    
        Returns:
            R_qs : quasi-statique combined response (for Gupta rule only)
            R_m2 : dynamic combined response of all oscillators in square
    
        """
        self._type_comb = comb_mode["TYPE"]
        self._amors = amors
        self._freqs = freqs
        self._H = None

        if self._type_comb == "GUPTA":
            if type_analyse != "MONO_APPUI":
                raise Exception("Pour la combinaison type Gupta n'existe qu'en cas de mono-appui")
            self._freq1, self._freq2 = comb_mode["FREQ_1"], comb_mode["FREQ_2"]
            self._alpha_r = None

        elif self._type_comb in ("DSC", "NRC_DSA"):
            self._duree = comb_mode["DUREE"]
            if self._duree is None:
                raise Exception("Il faut renseigner DUREE")
        
        elif self._type_comb == "NRC_GROUPING":
            self._groups = None

    @staticmethod
    def cqc_array(amors, freqs):
        """Matrix des coefficients CQC
    
        Arguments:
            amors           : list of damping coefficients
            freqs           : list of frequencies
    
        Returns:
            H: matrix of correlation between modes for CQC rule
        """
        nbmode = len(freqs)
        H = np.zeros((nbmode, nbmode))
        if np.any(amors == 0):
            H[np.diag_indices_from(H)] = 1
            return H
        omega = 2 * np.pi * freqs
        for i, (amor_i, w_i) in enumerate(zip(amors, omega)):
            for j, (amor_j, w_j) in enumerate(zip(amors, omega)):
                wij = w_i * w_j
                amor_ij = amor_i * amor_j
                H[i, j] = (8 * sqrt(amor_ij * wij) * (amor_i * w_i + amor_j * w_j) * wij) / (
                    (w_i**2 - w_j**2) ** 2
                    + 4 * amor_ij * wij * (w_i**2 + w_j**2)
                    + 4 * (amor_i**2 + amor_j**2) * wij**2
                )
        # return
        return H
    
    @staticmethod
    def dsc_array(amors, freqs, s):
        """Matrix des coefficients DSC
    
        Arguments:
            amors           : list of damping coefficients
            freqs           : list of frequencies
            s               : duration in second
    
        Returns:
            H: matrix of correlation between modes for DSC rule
        """
        if np.any(amors == 0) or np.any(amors > 1):
            raise Exception("Impossible de faire TYPE_COMB DSC avec des amortissements nuls ou > 1")
        nbmode = len(freqs)
        H = np.zeros((nbmode, nbmode))
        omega = 2 * np.pi * freqs
        # correlation coefficients for DSC rule
        for i, (amor_i, w_i) in enumerate(zip(amors, omega)):
            for j, (amor_j, w_j) in enumerate(zip(amors, omega)):
                H[i, j] = 1 / (
                    1
                    + (w_i * sqrt(1 - amor_i**2) - w_j * sqrt(1 - amor_j**2)) ** 2
                    / (w_i * (amor_i + 2 / (s * w_i)) + w_j * (amor_j + 2 / (s * w_j))) ** 2
                )
        # return
        return H
    
    @staticmethod
    def nrc_ten_percent_array(freqs):
        """Matrix des coefficients DPC (selon NRC_TEN_PERCENT)
    
        Arguments:
            freqs           : list of frequencies
    
        Returns:
            H: matrix of correlation between modes for NRC_TEN_PERCENT rule
        """
    
        nbmode = len(freqs)
        H = np.zeros((nbmode, nbmode))
        # correlation coefficients for DSC rule
        # using round() to invoid approximated values of frequencies
        for i in range(nbmode):
            for j in range(nbmode):
                if freqs[i] <= freqs[j] :
                    if freqs[j] <= round(1.1*freqs[i], 5):
                        H[i, j] = 1
                else:
                    if freqs[i] <= round(1.1*freqs[j], 5):
                        H[i, j] = 1
        # return
        return H

    @property
    def alpha_r(self):
        if self._alpha_r is None:
            f1 = min(self._freq1, self._freq2)
            f2 = max(self._freq1, self._freq2)
            # check if f1 < f2
            if f1 > f2:
                raise Exception(
                    "Pour la combinaison type Gupta '{f1}' doit être inférieure à '{f2}'".format(
                        f1=f1, f2=f2
                    )
                )
            self._alpha_r = np.expand_dims(np.log(self._freqs / f1) / np.log(f2 / f1), axis=1)
            self._alpha_r[self._freqs < f1] = 0
            self._alpha_r[self._freqs > f2] = 1
        return self._alpha_r

    @property
    def H(self):
        if self._H is None:
            if self._type_comb in ("CQC", "GUPTA"):
                self._H = self.cqc_array(self._amors, self._freqs)
                if self._type_comb == "GUPTA":
                    coeff_p = np.sqrt(1 - self.alpha_r**2)
                    self._H  = self._H * coeff_p.T * coeff_p
            elif self._type_comb in ("DSC", "NRC_DSA"):
                self._H = self.dsc_array(self._amors, self._freqs, self._duree)
            elif self._type_comb == "NRC_TEN_PERCENT":
                self._H = self.nrc_ten_percent_array(self._freqs)
            self._H[np.diag_indices_from(self._H)] /= 2
        return self._H

    @property
    def groups(self):
        """
        Create the modal group list by analyzing the frequencies
        **bottom-up** as required by RG 1.92 Rev.1
        """
        if self._groups is None:
            self._groups = []
            for fidx, ff in enumerate(self._freqs):
                if fidx == 0:
                    # Initialize first group with the first mode
                    self._groups.append([fidx])
                    continue
                #
                # "Current" group reference frequency
                ff_ref = self._freqs[self._groups[-1][0]]
                #
                if (ff-ff_ref)/ff_ref > 0.1:
                    # Start a new group if the 10% is exceeded
                    self._groups.append([fidx])
                else:
                    # Else, just append the mode index to the preceding group
                    self._groups[-1].append(fidx)
        return self._groups

    def get(self, R_mi):
        """
        Args:
            R_mi: list of all vectors corresponding to all modal responses (mode by mode)

        Returns:
            R_qs : quasi-statique combined response (for Gupta rule only)
            R_m2 : dynamic combined response of all oscillators in square
        """

        if self._type_comb == "SRSS":
            R_m2 = np.sum(R_mi**2, axis=0)
        elif self._type_comb == "ABS":
            R_m2 = np.sum(np.abs(R_mi), axis=0) ** 2
        elif self._type_comb in ("CQC", "DSC", "GUPTA"):
            H = self.H
            R_m2 = 0
            for i, r_i in enumerate(R_mi):
                for j, r_j in enumerate(R_mi[i:], i):
                    R_m2 += 2 * H[i, j] * r_i * r_j
            R_m2 = np.maximum(R_m2, 0)
        elif self._type_comb in ("NRC_DSA", "NRC_TEN_PERCENT"):
            H = self.H
            R_m2 = 0
            for i, r_i in enumerate(R_mi):
                for j, r_j in enumerate(R_mi[i:], i):
                    R_m2 += 2 * H[i, j] * np.abs(r_i * r_j)
            R_m2 = np.maximum(R_m2, 0)
        elif self._type_comb == "DPC":
            # neighbor modes (frequence close until 10%) will be combined by abs
            f0, r_r = self._freqs[0], R_mi[0]
            l_abs_R_r = [np.abs(r_r)]
            for i_freq in range(1, len(self._freqs)):
                freq = self._freqs[i_freq]
                r_r = R_mi[i_freq]
                fm_ = 0.5 * (f0 + freq)
                if abs(freq - f0) / fm_ >= 0.10:  # discrepance in freq bigger than 10%
                    f0 = freq
                    l_abs_R_r.append(np.abs(r_r))
                else:  # # discrepance in freq bigger than 10% --> sum by abs
                    f0 = freq
                    # combinied by abs value
                    l_abs_R_r[-1] += np.abs(r_r)
            return sum(r_x**2 for r_x in l_abs_R_r), 0
    
        elif self._type_comb == "NRC_GROUPING":
            # The regular Squared-Sum contribution (no correlation part)
            sum_rk_sq = np.sum(R_mi**2, axis=0)
            # The correlated part
            sum_rl_rm = np.zeros_like(sum_rk_sq)
            for group in self.groups:
                nbmodes = len(group)
                if nbmodes == 1:
                    continue
                for ii in range(0, nbmodes-1):
                    rl = R_mi[group[ii]]
                    for jj in range(ii+1, nbmodes):
                        rm = R_mi[group[jj]]
                        sum_rl_rm += 2.*np.abs(rl*rm)

            R_m2 = sum_rk_sq + sum_rl_rm
        else:
            raise Exception("Le TYPE '{type_comb}' n'est pas reconnu".format(type_comb=self._type_comb))

        if self._type_comb == "GUPTA":
            R_qs = np.sum(R_mi * self.alpha_r, axis=0)
        else:
            R_qs = np.zeros((np.shape(R_m2)))

        return R_m2, R_qs


def comb_appui_corr(COMB_MULT_APPUI_CORR, R_mi):
    """Combines modals responses

    Args:
        COMB_MULT_APPUI_CORR: rule for combinaison of correlated supports
        R_mi: list of all vectors corresponding to all correlated supports

    Returns:
        R_m : combined response

    """
    # combinaison rule for correlated supports
    type_comb = COMB_MULT_APPUI_CORR
    # preparing for combinaison
    R_m = 0
    if type_comb == "QUAD":
        R_m = np.sqrt(np.sum(R_mi**2, axis=0))
    elif type_comb == "LINE":
        R_m = np.sum(R_mi, axis=0)
    elif type_comb == "ABS":
        R_m = np.sum(np.abs(R_mi), axis=0)
    # return
    return R_m


def comb_directions(type_comb_dir, l_R_x):
    """Combine directional responses according to the given type of combination

    Args:
        type_comb_dir (str): type of combination
        l_R_x (list): list of the directional responses

    Returns:
        R_xyz (ndarray): combined array
    """
    # ("lancer combi directionnelle")
    nb_direction = len(l_R_x)
    if nb_direction > 1:
        if type_comb_dir == "QUAD":
            R_xyz = np.sqrt(sum(r_x**2 for r_x in l_R_x))
            R_newmark_all = []

        elif type_comb_dir == "ABS":
            R_xyz = np.sum(np.abs(l_R_x), axis=0)
            R_newmark_all = []

        else:  # type_comb_dir == "NEWMARK":
            # 24/8/2 combinations if 3/2/1 directions :
            # R = [± EX  ± 0, 4 * EY  ± 0, 4 * EZ]
            # R = [± 0, 4 * EX  ± EY  ± 0, 4 * EZ]
            # R = [± 0, 4 * EX  ± 0, 4 * EY  ± EZ]
            newmark = np.zeros(len(l_R_x[0]))
            R_newmark_all = []

            def circular_shifts(p):
                return [np.roll(p, -i) for i in range(len(p))]  # idem as in more_itertools

            for ijk in circular_shifts(range(nb_direction)):
                r_xyz = [l_R_x[idx] for idx in ijk]
                for exponents in product([0, 1], repeat=nb_direction):
                    comb_nk = sum(
                        (-1) ** expo * coef * r_x
                        for expo, coef, r_x in zip(exponents, [1, 0.4, 0.4], r_xyz)
                    )
                    newmark = np.maximum(newmark, comb_nk)
                    R_newmark_all.append(comb_nk)
            R_xyz = newmark
    elif nb_direction == 1:
        R_xyz = l_R_x[0]
        R_newmark_all = [l_R_x[0], -1.0 * l_R_x[0]]
    # return
    return R_xyz, R_newmark_all


def corr_pseudo_mode_mono(
    option,
    pseudo_mode,
    amors,
    freq_coup,
    pr_wr2_phi,
    w_r,
    spectre_dir,
    spectre_nappe,
    spectre_coeff,
    corr_freq,
    spectre_nature,
):
    """correction by pseudo-mode/mode statique

    Args:
        option      : option of field to combine
        pseudo_mode : pseudo_mode (mode_meca)
        amors       : list of damping coefficients
        freq_coup   : scalar value of cutting frequency at which ZPA is read
        pr_wr2_phi  : list of produit rho*phi/omega^2
        w_r         : list of omega = 2*pi*freq
        spectre_dir : direction
        spectre_nappe: "nappe" objects for spectrum
        spectre_coeff: scale factor for spectrum
        corr_freq   : corr_freq option (OUI or NON)
        spectre_nature: nature of spectrum

    Returns:
        R_c (ndarray): response by correction of pseudo-mode
    """
    # Run corr_pseudo_mode_mono"
    # ZPA at cut frequency
    S_r_freq_coup = spectre_nappe(amors[-1], freq_coup) * spectre_coeff
    # correction of spectrum if corr_freq is "OUI"
    if corr_freq == "OUI":
        S_r_freq_coup = S_r_freq_coup * np.sqrt(1 - amors[-1] ** 2)
    # nature of spectrum to ACCE spectrum
    if spectre_nature == "DEPL":
        S_r_freq_coup = ((2 * np.pi * freq_coup) ** 2) * S_r_freq_coup
    elif spectre_nature == "VITE":
        S_r_freq_coup = 2 * np.pi * freq_coup * S_r_freq_coup
    elif spectre_nature == "ACCE":
        pass
    # search for index (NUME_CMP) in MODE_STATIQUE corresponding to direction
    ps_noeud_cmp = pseudo_mode.getAccessParameters()["NOEUD_CMP"]
    ps_nume_mode = pseudo_mode.getAccessParameters()["NUME_MODE"]
    if "OX" in spectre_dir:
        index_pseudo_mode = ps_nume_mode[ps_noeud_cmp.index("ACCE    X")]
    elif "OY" in spectre_dir:
        index_pseudo_mode = ps_nume_mode[ps_noeud_cmp.index("ACCE    Y")]
    elif "OZ" in spectre_dir:
        index_pseudo_mode = ps_nume_mode[ps_noeud_cmp.index("ACCE    Z")]
    else:
        raise Exception(
            "Direction '{spectre_dir}' ne correspond pas à aucun axe global".format(
                spectre_dir=spectre_dir
            )
        )
    if index_pseudo_mode is None:
        raise Exception(
            "Direction '{spectre_dir}' n'existe pas dans la base pseudo-modale".format(
                spectre_dir=spectre_dir
            )
        )
    # get field in mode_statique
    if option in ["DEPL", "REAC_NODA", "FORC_NODA"]:
        phi_ps = pseudo_mode.getField(option, index_pseudo_mode).getValues()
        # pseudo-mode
        R_c = (phi_ps - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
    elif option in ["EFGE_ELNO", "EGRU_ELNO", "SIEF_ELGA", "SIGM_ELNO", "SIPO_ELNO", "SIEF_ELNO"]:
        phi_ps = pseudo_mode.getField(option, index_pseudo_mode).getValues()
        # pseudo-mode
        R_c = (phi_ps - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
    if option == "VITE":  # on accepte que VITE = DEPL * omega (pseudo-vitesse relative)
        phi_ps = pseudo_mode.getField("DEPL", index_pseudo_mode).getValues()
        UTMESS("F", "SEISME_10", valk=option)
        # pseudo-mode
        R_c = (phi_ps * w_r - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
    if option == "ACCE_ABSOLU":
        phi_ps = pseudo_mode.getField("DEPL", index_pseudo_mode).getValues()
        # this component is zero because this correction is automatic for mono-appui
        R_c = np.zeros(np.shape(phi_ps))
    # return
    return R_c, S_r_freq_coup


def corr_pseudo_mode_mult(
    option,
    pseudo_mode,
    l_group_no,
    mesh,
    direction,
    amors,
    freq_coup,
    pr_wr2_phi,
    w_r,
    spectre_nappe,
    spectre_coeff,
    corr_freq,
    spectre_nature,
):
    """correction par pseudo-mode/mode statique pour multi_appui
    Args:
        option      : option of field to combine
        pseudo_mode : pseudo_mode (mode_meca)
        l_group_no  : list of all group_no of support
        mesh    : mesh extracted from mode_meca by mode_meca.getMesh()
        direction   : direction
        amors       : list of damping coefficients
        freq_coup   : scalar value of cutting frequency at which ZPA is read
        pr_wr2_phi  : list of produit rho*phi/omega^2
        w_r         : list of omega = 2*pi*freq
        spectre_nappe: "nappe" objects for spectrum
        spectre_coeff: scale factor for spectrum
        corr_freq   : corr_freq option (OUI or NON)
        spectre_nature: nature of spectrum

    Returns:
        R_c (ndarray): response by correction by pseudo-mode
    """
    # Run corr_pseudo_mode_mult
    # ZPA at the cutting frequency
    S_r_freq_coup = spectre_nappe(amors[-1], freq_coup) * spectre_coeff
    # correction of spectrum by corr_freq
    if corr_freq == "OUI":
        S_r_freq_coup = S_r_freq_coup * np.sqrt(1 - amors[-1] ** 2)
    # correction of spectrum by nature of spectrum
    if spectre_nature == "DEPL":
        S_r_freq_coup = ((2 * np.pi * freq_coup) ** 2) * S_r_freq_coup
    elif spectre_nature == "VITE":
        S_r_freq_coup = 2 * np.pi * freq_coup * S_r_freq_coup
    elif spectre_nature == "ACCE":
        pass
    # searche for index (NUME_CMP) in MODE_STATIQUE / PSEUDO_MODE
    # list of node_num and node_name in group_no
    nodes_num_per_grno, nodes_name_per_grno = get_nodes(mesh, GROUP_NO=l_group_no)
    # list of all noeud_cmp
    l_noeud_cmp = []
    for i_node in range(len(nodes_num_per_grno)):
        node_name = nodes_name_per_grno[i_node]
        noeud_cmp = node_name.ljust(8) + direction.replace("O", "D")
        l_noeud_cmp.append(noeud_cmp)
    ps_nume_modes = pseudo_mode.getAccessParameters()["NUME_MODE"]
    ps_noeud_cmps = pseudo_mode.getAccessParameters()["NOEUD_CMP"]
    if all(noeud_cmp in ps_noeud_cmps for noeud_cmp in l_noeud_cmp):
        # ps_nume_mode = ps_nume_modes[ps_noeud_cmps.index(noeud_cmp)]
        l_phi_ps = []
        for noeud_cmp in l_noeud_cmp:
            ps_nume_mode = ps_nume_modes[ps_noeud_cmps.index(noeud_cmp)]
            l_phi_ps.append(pseudo_mode.getField(option, ps_nume_mode))
        phi_ps = l_phi_ps[0].copy()
        for x in l_phi_ps[1:]:
            phi_ps += x
        # get static mode for component (noeud_cmp)
        if option in ["DEPL", "REAC_NODA", "FORC_NODA"]:
            R_c_noeud = (phi_ps.getValues() - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
        elif option in [
            "EFGE_ELNO",
            "EGRU_ELNO",
            "SIEF_ELGA",
            "SIGM_ELNO",
            "SIPO_ELNO",
            "SIEF_ELNO",
        ]:
            R_c_noeud = (phi_ps.getValues() - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
        if option == "VITE":  # pseudo-mode is not allowed
            UTMESS("F", "SEISME_10", valk=option)
            R_c_noeud = (phi_ps.getValues() * w_r - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
        if (
            option == "ACCE_ABSOLU"
        ):  # correction by pseudo-mode is not allowed for ACCE_ABSOLU in mutl_appui
            UTMESS("F", "SEISME_10", valk=option)
            R_c_noeud = np.zeros(phi_ps.size())
    else:  # noued_cmp is not calculated in pseudo-mode --> need to be done
        raise Exception(
            "NOEUD_CMP: {0} n'existe pas dans la base MODE_STATIQUE - mot-clé PSEUDO_MODE".format(
                noeud_cmp
            )
        )
    # return
    return R_c_noeud, S_r_freq_coup


def corr_pseudo_mode_enveloppe(
    option,
    pseudo_mode,
    amors,
    freq_coup,
    pr_wr2_phi,
    w_r,
    spectre_dir,
    S_r_freq_coup,
    corr_freq,
    spectre_nature,
):
    """correction by pseudo-mode/mode statique
    for TYPE_ANALYSE = ENVELOPPE
    same as 'corr_pseudo_mode_mono' but arguments = S_r_freq_coup instead of "spectre_nappe" and "spectre_coef"

    Args:
        option      : option of field to combine
        pseudo_mode : pseudo_mode (mode_meca)
        amors       : list of damping coefficients
        freq_coup   : scalar value of cutting frequency at which ZPA is read
        pr_wr2_phi  : list of produit rho*phi/omega^2
        w_r         : list of omega = 2*pi*freq
        spectre_dir : direction
        S_r_freq_coup : maximal value of ZPA at the cutting frequency
        corr_freq   : corr_freq option (OUI or NON)
        spectre_nature: nature of spectrum

    Returns:
        R_c (ndarray): response by correction of pseudo-mode
    """
    # Run corr_pseudo_mode_enveloppe"
    # # ZPA at cut frequency
    # S_r_freq_coup = spectre_nappe(amors[-1], freq_coup) * spectre_coeff
    # correction of spectrum if corr_freq is "OUI"
    if corr_freq == "OUI":
        S_r_freq_coup = S_r_freq_coup * np.sqrt(1 - amors[-1] ** 2)
    # nature of spectrum to ACCE spectrum
    if spectre_nature is "DEPL":
        S_r_freq_coup = ((2 * np.pi * freq_coup) ** 2) * S_r_freq_coup
    elif spectre_nature is "VITE":
        S_r_freq_coup = 2 * np.pi * freq_coup * S_r_freq_coup
    elif spectre_nature is "ACCE":
        S_r_freq_coup = S_r_freq_coup
    # search for index (NUME_CMP) in MODE_STATIQUE corresponding to direction
    ps_noeud_cmp = pseudo_mode.getAccessParameters()["NOEUD_CMP"]
    ps_nume_mode = pseudo_mode.getAccessParameters()["NUME_MODE"]

    if "OX" in spectre_dir:
        index_pseudo_mode = ps_nume_mode[ps_noeud_cmp.index("ACCE    X")]
        components = "DX"
    elif "OY" in spectre_dir:
        index_pseudo_mode = ps_nume_mode[ps_noeud_cmp.index("ACCE    Y")]
        components = "DY"
    elif "OZ" in spectre_dir:
        index_pseudo_mode = ps_nume_mode[ps_noeud_cmp.index("ACCE    Z")]
        components = "DZ"
    else:
        raise Exception(
            "Direction '{spectre_dir}' ne correspond pas à aucun axe global".format(
                spectre_dir=spectre_dir
            )
        )
    if index_pseudo_mode is None:
        raise Exception(
            "Direction '{spectre_dir}' n'existe pas dans la base pseudo-modale".format(
                spectre_dir=spectre_dir
            )
        )
    # get field in mode_statique
    if option in ["DEPL", "REAC_NODA", "FORC_NODA"]:
        phi_ps = pseudo_mode.getField(option, index_pseudo_mode).getValues()
        # pseudo-mode
        R_c = (phi_ps - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
    elif option in ["EFGE_ELNO", "EGRU_ELNO", "SIEF_ELGA", "SIGM_ELNO", "SIPO_ELNO", "SIEF_ELNO"]:
        phi_ps = pseudo_mode.getField(option, index_pseudo_mode).getValues()
        # pseudo-mode
        R_c = (phi_ps - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
    if option == "VITE":  # on accepte que VITE = DEPL * omega (pseudo-vitesse relative)
        phi_ps = pseudo_mode.getField("DEPL", index_pseudo_mode).getValues()
        UTMESS("F", "SEISME_10", valk=option)
        # pseudo-mode
        R_c = (phi_ps * w_r - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
    if option == "ACCE_ABSOLU":
        phi_ps = pseudo_mode.getField("DEPL", index_pseudo_mode).getValues()
        # this component is zero because this correction is automatic for mono-appui
        R_c = np.zeros(np.shape(phi_ps))
    # return
    return R_c


def get_phis(mode_meca, option, nume_ordres):
    """get eigen-vector
    Args:
         mode_meca  : modale basis
         option     : fields to get, except for (VITE, ACCE_ABSOLU)
         nume_ordres: nume_ordre of feild to be gotten
    Returns:
        phis (ndarray)
    """
    # Run get_phis
    # check if field exists in modal basis (mode_meca)
    if option not in mode_meca.getFieldsNames():
        raise Exception(
            "Le champs '{option}' n'existe pas dans la base modale".format(option=option)
        )
    # preparing to get vector phi
    phis = []
    if all(nume_ordre in mode_meca.getIndexes() for nume_ordre in nume_ordres):
        for imode in nume_ordres:
            if option in ["DEPL", "REAC_NODA", "FORC_NODA"]:
                phis.append(mode_meca.getField(option, imode).getValues())
            elif option in [
                "EFGE_ELNO",
                "EGRU_ELNO",
                "SIEF_ELGA",
                "SIGM_ELNO",
                "SIPO_ELNO",
                "SIEF_ELNO",
            ]:
                phis.append(mode_meca.getField(option, imode).getValues())
    # eigen-vector
    phis = np.array(phis)
    # return
    return phis


def get_depl_mult_appui(depl_mult_appui):
    """get input form key_factor DEPL_MULT_APPUI
    Args:
         input of key_factor DEPL_MULT_APPUI by users
    Returns:
        [mode_stat, group_no_refe, D_e, dirs] for each support
        where
        mode_stat{nom_appui}    : static mode of structure (mode_stat)
        group_no_refe{nom_appui}: group_no of the unique node of reference (not used)
        D_e[nom_appui]          : displacement values
        dirs[nom_appui]         : directions of displacements
    """
    # run get_depl_mult_appui
    if depl_mult_appui is not None:
        mode_stat = {}
        group_no_refe = {}
        D_e = {}
        dirs = {}
        dir = []
        # get information on nom_appui
        for i in range(len(depl_mult_appui)):
            # get nom_appui
            nom_appui = depl_mult_appui[i].get("NOM_APPUI")
            D_e[nom_appui] = {}
        # for each nom_appui, get information about DDS
        for i in range(len(depl_mult_appui)):
            # get nom_appui
            nom_appui = depl_mult_appui[i].get("NOM_APPUI")
            # add mode_stat corresponding to nom_appui
            mode_stat[nom_appui] = depl_mult_appui[i].get("MODE_STAT")
            # add group_no_refe corresponding to nom_appui
            group_no_refe[nom_appui] = depl_mult_appui[i].get("GROUP_NO_REFE")
            # add D_e corresponding to nom_appui
            D_e[nom_appui] = {}
            if depl_mult_appui[i].get("DX") is not None:
                D_e[nom_appui]["OX"] = depl_mult_appui[i].get("DX")
                dir.append("OX")
            # axe Y
            if depl_mult_appui[i].get("DY") is not None:
                D_e[nom_appui]["OY"] = depl_mult_appui[i].get("DY")
                dir.append("OY")
            # axe Z
            if depl_mult_appui[i].get("DZ") is not None:
                D_e[nom_appui]["OZ"] = depl_mult_appui[i].get("DZ")
                dir.append("OZ")
            # if all DX, DY and DZ are not present
            if all(depl_mult_appui[i].get(xx) is None for xx in ["DX", "DY", "DZ"]):
                raise Exception("Aucun DDS n'est trouvé pour {0}".format(nom_appui))
            # add direction to nom_appui
            dirs[nom_appui] = dir
            # reset parameter dir
            dir = []
        # saving
        depl_mult_appuis = [mode_stat, group_no_refe, D_e, dirs]
        # get all directions of displacements input by users
        dir_standard = ["OX", "OY", "OZ"]
        D_e_dirs = []
        for i_dir_standard in range(len(dir_standard)):
            if any(dir_standard[i_dir_standard] in dir for dir in dirs):
                D_e_dirs.append(dir_standard[i_dir_standard])
    else:
        depl_mult_appuis = [[], [], [], []]
        D_e_dirs = []
    # return
    return depl_mult_appuis, D_e_dirs

class Resu:
    """
    Args:
        type_resu       : input of TYPE_RESU
        mode_meca       : input of MODE_MECA
        mesh            : the mesh
    """
    _response_by_type = {"VALE_QS": "PART_QS{}",
                         "VALE_DIRE": "DIRE{}",
                         "VALE_DDS": "PART_DDS{}",
                         "VALE_INER": "PART_INER{}",
                         "VALE_DYNA": "PART_DYNA{}",
                         "VALE_TOTA": "TOTA{}"}
    _newmark_by_type = {"VALE_TOTA": "NEWM_{}",
                        "VALE_DYNA": "DYNA_NEWN_{}",
                        "VALE_INER": "INER_NEWM_{}",
                        "VALE_DDS": "DDS_NEWM_{}",
                        "VALE_QS": "QS_NEWM_{}"}
    def __init__(self, type_resu, mode_meca, mesh):
        self._type_resu = type_resu
        self._mode_meca = mode_meca
        self._mesh = mesh
        self._result = MultipleElasticResult()
        self._nbIndexes = 0
        self._field_by_option= {}

    def _incr_index(self):
        self._nbIndexes +=1
        if self._nbIndexes==1:
            self._result.allocate(self._nbIndexes)
        else:
            self._result.resize(self._nbIndexes)

    def _setField(self, option, values, name):
        field = self._copy_field(option)
        field.setValues(values)
        self._incr_index()
        self._result.setField(field, option, value=name, para="NOM_CAS")

    def _copy_field(self, option):
        """Create identical field of option field"""
        field = self._field_by_option.get(option)
        if not field:
            if option in ("VITE", "ACCE_ABSOLU"):
                option2 = "DEPL"
            else:
                option2 = option
            field = self._mode_meca.getField(option2, 1).copy()
            self._field_by_option[option] = field
        return field

    def _get_type(self, vale_type):
        """return type_resu corresponding to vale_type"""
        for type_resu in self._type_resu:
            if type_resu.get("TYPE")==vale_type:
                return type_resu

    def add_dire_response(self, option, direction, R, vale_type):
        """print out directional response
        Args:
            option         : field to be printed
            direction      : direction
            R              : directional response
            vale_type      : type of response
        """
        type_resu = self._get_type(vale_type)
        if type_resu is None:
            return

        list_axe = type_resu.get("LIST_AXE")
        direction = direction.replace("O", "")
        if list_axe and direction in list_axe:
            self._setField(option, R, self._response_by_type[vale_type].format("_" + direction))

    def add_response(self, option, R, R_newmark_all, vale_type):
        """print out response
        Args:
            option          : field to be printed
            R               : total response
            R_newmark_all   : response by NEWMARK
            vale_type       : type of response
        """
        type_resu = self._get_type(vale_type)
        if type_resu is None:
            return

        self._setField(option, R, self._response_by_type[vale_type].format(""))
        if type_resu.get("NEWMARK") == "OUI":
            for i, r_newmark in enumerate(R_newmark_all):
                self._setField(option, r_newmark, self._newmark_by_type[vale_type].format(i + 1))

    def add_spectral_response(self, option, direction, R, nume_ordres_resu, appui=""):
        """print out spectral response (mono or multi appuis)
        Args:
            option          : field to be printed
            direction       : direction
            R               : spectral response
            nume_ordres_resu: list of nume ordre in R
            appui           : support to be given for multi_appui
        """
        type_resu = self._get_type("VALE_SPEC")
        if type_resu is None:
            return

        _, nume_ordres, _, _ = filter_ordre_freq(self._mode_meca.LIST_PARA(), type_resu)
        list_axe = type_resu.get("LIST_AXE")
        list_appui = type_resu.get("LIST_APPUI")
        tout_appui = type_resu.get("TOUT_APPUI")
        direction = direction.replace("O", "")
        if direction in list_axe and (not appui or tout_appui == "OUI" or appui in list_appui):
            for nume_ordre in nume_ordres:
                i_nume_ordre = (nume_ordres_resu.tolist()).index(nume_ordre)
                self._setField(option, R[i_nume_ordre], f"SPEC_{nume_ordre}{direction}{appui}")

    def get(self):
        model = self._mode_meca.getModel()
        if model:
            self._result.setModel(model)
        else:
            self._result.setMesh(self._mesh)
        cara_elem = self._mode_meca.getElementaryCharacteristics()
        if cara_elem:
            self._result.setElementaryCharacteristics(cara_elem)
        return self._result


def comb_sism_modal_ops(self, **args):
    """Execute the command COMB_SISM_MODAL.

    Arguments:
        **args (dict): User's keywords.

    Returns:
        resu: result for linear modal combinaison
    """
    mode_meca = args.get("MODE_MECA")
    amor_reduit = args.get("AMOR_REDUIT")
    list_amor = args.get("LIST_AMOR")
    amor_gene = args.get("AMOR_GENE")
    mode_corr = args.get("MODE_CORR")
    pseudo_mode = args.get("PSEUDO_MODE")
    freq_coup_in = args.get("FREQ_COUP")
    type_analyse = args.get("TYPE_ANALYSE")
    spectre_in = args.get("SPECTRE")
    appuis_in = args.get("APPUIS")
    comb_mode = args.get("COMB_MODE")
    comb_direction = args.get("COMB_DIRECTION")
    group_appui_correle = args.get("GROUP_APPUI_CORRELE")
    # CUMUL_INTRA : Combination of the contributions of each response of support inside a group
    cumul_intra = args.get("CUMUL_INTRA")
    # CUMUL_INTER : Combination of the contributions of each group of supports
    cumul_inter = args.get("CUMUL_INTER")
    comb_dds_correle = args.get("COMB_DDS_CORRELE")
    depl_mult_appui = args.get("DEPL_MULT_APPUI")
    type_resu = args.get("TYPE_RESU")
    verbosity = args["INFO"]
    setFortranLoggingLevel(verbosity)
    # --------------------------------------------------------------------------
    # exploring mode_meca
    # --------------------------------------------------------------------------
    # mesh in mode_meca
    mesh = mode_meca.getMesh()
    # list of parameters in mode_meca
    list_para = mode_meca.LIST_PARA()
    # get numeber of orders of modes to be combined
    freqs, nume_ordres, nume_modes, gene_masses = filter_ordre_freq(list_para, args)

    # Get damping coefficients
    amor_reduit = get_amor_reduit(mode_meca, nume_ordres, amor_reduit, list_amor, amor_gene)
    # if number of damping coefficient is less than modes to combine
    # add last value of damping coefficient to complete the damping vector
    nb_modes = len(nume_modes)
    amors = np.array(list(amor_reduit) + [amor_reduit[-1]] * (nb_modes - len(amor_reduit)))
    #  participation factor in mono_appui
    l_fact_partici = []
    for key_fact in ["FACT_PARTICI_DX", "FACT_PARTICI_DY", "FACT_PARTICI_DZ"]:
        l_fact_partici_par_dir = []
        for i_nume_ordre in range(len(nume_ordres)):
            # search for index of number of order of considered mode
            index_nume_ordre = list_para["NUME_ORDRE"].index(nume_ordres[i_nume_ordre])
            l_fact_partici_par_dir.append(list_para[key_fact][index_nume_ordre])
        l_fact_partici.append(np.array(l_fact_partici_par_dir))
    # Masse effective modale
    l_masse_effe = []
    for key_fact in ["MASS_EFFE_DX", "MASS_EFFE_DY", "MASS_EFFE_DZ"]:
        l_masse_effe_par_dir = []
        for i_nume_ordre in range(len(nume_ordres)):
            # search for index of number of order of considered mode
            index_nume_ordre = list_para["NUME_ORDRE"].index(nume_ordres[i_nume_ordre])
            l_masse_effe_par_dir.append(list_para[key_fact][index_nume_ordre])
        l_masse_effe.append(np.array(l_masse_effe_par_dir))
    # Masse effective modale
    l_masse_effe_un = []
    for key_fact in ["MASS_EFFE_UN_DX", "MASS_EFFE_UN_DY", "MASS_EFFE_UN_DZ"]:
        l_masse_effe_un_par_dir = []
        for i_nume_ordre in range(len(nume_ordres)):
            # search for index of number of order of considered mode
            index_nume_ordre = list_para["NUME_ORDRE"].index(nume_ordres[i_nume_ordre])
            l_masse_effe_un_par_dir.append(list_para[key_fact][index_nume_ordre])
        l_masse_effe_un.append(np.array(l_masse_effe_un_par_dir))
    # Get spectres
    if type_analyse == "MONO_APPUI":
        spectres = get_spectres_mono_appui(spectre_in)
    elif type_analyse == "MULT_APPUI":
        spectres = get_spectres_mult_appui(spectre_in, mesh)
        # search for all directions presented by users
        dir_all_old = []
        for i_appui in range(len(spectres)):
            dir_all_old += spectres[i_appui][0]
        dir_all_old = list(set(dir_all_old))
        # compared to standard directions
        dir_standard = ["OX", "OY", "OZ"]
        dir_all = []
        for i_dir_standard in range(len(dir_standard)):
            if dir_standard[i_dir_standard] in dir_all_old:
                dir_all.append(dir_standard[i_dir_standard])
    elif type_analyse == "ENVELOPPE":
        spectres = get_spectres_mult_appui(spectre_in, mesh)
    # Get support (APPUI)
    if type_analyse == "MULT_APPUI":
        # Get support (APPUI)
        nom_appuis_all, l_nodes_num_all, l_nodes_name_all, l_group_no_all = get_appuis(
            appuis_in, mesh
        )
    # Get all groupes of support (group_appui)
    if type_analyse == "MULT_APPUI":
        group_appui_all, nom_group_appui_all = get_group_appuis(spectres, group_appui_correle)
    # Get imposed displacement DDS
    if type_analyse == "MULT_APPUI":
        depl_mult_appuis, D_e_dirs_all = get_depl_mult_appui(depl_mult_appui)
    # --------------------------------------------------------------------------
    # Output result preparing
    # --------------------------------------------------------------------------

    resu = Resu(type_resu, mode_meca, mesh)

    # --------------------------------------------------------------------------
    # Combinaison
    # --------------------------------------------------------------------------
    # Iteration on option result to combine
    i_option = 0

    comb_modal_response = CombModalResponse(comb_mode, type_analyse, amors, freqs)
    for option in args["OPTION"]:
        i_option += 1
        if type_analyse == "MONO_APPUI":
            # step 1: Get eigen-vector to combine from mode_meca
            if option not in ["VITE", "ACCE_ABSOLU"]:
                phis = get_phis(mode_meca, option, nume_ordres)
            elif option == "VITE" or option == "ACCE_ABSOLU":
                phis = get_phis(mode_meca, "DEPL", nume_ordres)
            else:
                raise Exception("OPTION '{option}' n'est pas pris en compte".format(option=option))
            # step 2: spectral value
            l_R_x = []  # list of directional total result
            l_part_d = []  # list of directional result part dynamique
            l_part_s = []  # list of directional result part pseudo-statique
            l_R_prim = []  # list of RCCM part primaire
            # preparing for printing out to INFO
            l_SA = {}  # list of info about read spectra at eigen-frequencies by direction
            l_pseudo = {}  # list of info about correction by pseudo mode by direction
            # iteration on direction
            for i_dir in range(len(spectres[0])):
                # eigen-pulsation before corr_freq
                w_r = 2 * np.pi * freqs
                # Spectrum interpolation
                [spectre_dir, spectre_nappe, spectre_coeff, spectre_corr_freq, spectre_nature] = [
                    spectres[i][i_dir] for i in range(5)
                ]
                # Correction for frequency by corr_freq
                if spectre_corr_freq == "OUI":
                    correct = np.sqrt(1 - amors**2)
                else:
                    correct = 1
                # Pulsation afeter corr_freq
                w_r *= correct
                # Spectrale values at eigen-frequencies
                S_r_freq = []
                for i_freq in range(len(freqs)):
                    S_r_freq.append(spectre_nappe(amors[i_freq], freqs[i_freq]) * spectre_coeff)
                # Correction by corr_freq for spectrum
                S_r_freq *= correct
                # Cutting frequency
                if freq_coup_in is not None:
                    freq_coup = freq_coup_in
                else:
                    freq_coup = freqs[-1]
                #  Correction of spectrum by nature of spectrum
                if spectre_nature == "DEPL":
                    S_r_freq = ((2 * np.pi * freqs) ** 2) * S_r_freq
                elif spectre_nature == "VITE":
                    S_r_freq = 2 * np.pi * freqs * S_r_freq
                elif spectre_nature == "ACCE":
                    pass
                # Participation factor by direction
                if "OX" in spectre_dir:
                    fact_partici = l_fact_partici[0]
                    components = "DX"
                elif "OY" in spectre_dir:
                    fact_partici = l_fact_partici[1]
                    components = "DY"
                elif "OZ" in spectre_dir:
                    fact_partici = l_fact_partici[2]
                    components = "DZ"
                else:
                    raise Exception(
                        "Direction '{spectre_dir}' ne correspond pas à aucun axe global".format(
                            spectre_dir=spectre_dir
                        )
                    )
                # Spectral response
                if option not in ["VITE", "ACCE_ABSOLU"]:
                    R_mi_all = (S_r_freq * fact_partici / w_r**2)[:, None] * phis
                    pr_wr2_phi_all = (fact_partici / w_r**2)[:, None] * phis
                    pr_wr2_phi_c_all = (fact_partici / (2 * np.pi * freqs) ** 2)[
                        :, None
                    ] * phis  # interrogration ??? pq ne pas utiliser omega corrige?
                elif option == "VITE":  # ici: phis correspond à DEPL
                    R_mi_all = (S_r_freq * fact_partici / w_r)[:, None] * phis
                    pr_wr2_phi_all = (fact_partici / w_r)[:, None] * phis
                    pr_wr2_phi_c_all = (fact_partici / w_r)[:, None] * phis
                elif option == "ACCE_ABSOLU":  # ici: phis correspond à DEPL
                    R_mi_all = (S_r_freq * fact_partici)[:, None] * phis
                    pr_wr2_phi_all = (fact_partici)[:, None] * phis
                    pr_wr2_phi_c_all = (fact_partici)[:, None] * phis
                # in case where the first mode is bigger than cutting frequency
                if freq_coup is not None and freq_coup >= freqs[0]:
                    R_mi = R_mi_all
                    pr_wr2_phi = pr_wr2_phi_all
                    pr_wr2_phi_c = pr_wr2_phi_c_all
                elif freq_coup < freqs[0]:
                    R_mi = np.zeros(np.shape(phis))
                    pr_wr2_phi = np.zeros(np.shape(phis))
                    pr_wr2_phi_c = np.zeros(np.shape(phis))
                    # Raise alarm for zero mode to be considered before cutting frequency
                    UTMESS("A", "SEISME_96", valr=freq_coup)
                # Print output for spectral value for each mode and direction
                resu.add_spectral_response(option, spectre_dir, R_mi_all, nume_ordres)
                # step 3: modal combinaison
                # Get input COMB_MODE
                R_m2, R_qs = comb_modal_response.get(R_mi)
                # Automatic correction for ACCE_ABSOLU in mono-appui
                if option == "ACCE_ABSOLU":
                    # field of unit value for acce_absolu
                    acce_unitaire = mode_meca.getField("DEPL", 1).copy()
                    acce_unitaire.setValues({components: 1.0}, [])
                    S_r_freq_coup = spectre_nappe(amors[-1], freq_coup) * spectre_coeff
                    R_tt = (acce_unitaire.getValues() - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
                    # add to combined modale responses in square
                    R_m2 += R_tt**2
                # modale response
                R_m = np.sqrt(R_m2)
                # step 4 : Entrainement zero pour mon_appui
                R_e2 = np.zeros(np.shape(R_m2))
                # step 5 : pseudo-mode response
                if mode_corr == "OUI":
                    # check if cutting frequency is present
                    if freq_coup_in is None:
                        UTMESS("A", "SEISME_95", valr=freq_coup)
                    # calculate pseudo-mode
                    R_c, S_r_freq_coup = corr_pseudo_mode_mono(
                        option,
                        pseudo_mode,
                        amors,
                        freq_coup,
                        pr_wr2_phi_c,
                        w_r,
                        spectre_dir,
                        spectre_nappe,
                        spectre_coeff,
                        spectre_corr_freq,
                        spectre_nature,
                    )
                    # save for INFO
                    l_pseudo[spectre_dir] = [freq_coup, S_r_freq_coup]
                else:
                    R_c = np.zeros(np.shape(R_m2))
                    S_r_freq_coup = None
                # step 6 : reponse by direction
                # total
                R_x = np.sqrt(R_m2 + (R_qs + R_c) ** 2 + R_e2)
                # inertial part (part primaire)
                R_prim = np.sqrt(R_m2 + (R_qs + R_c) ** 2)
                # add total directionnal responses
                l_R_x.append(R_x)
                # POST_ROCHE/ part dynamique et pseudo statique
                l_part_d.append(R_m)
                l_part_s.append(R_c)
                # RCCM part primaire
                l_R_prim.append(R_prim)
                # Print out response for mono_appui by direction
                resu.add_dire_response(option, spectre_dir, R_x, "VALE_DIRE")
                resu.add_dire_response(option, spectre_dir, R_m, "VALE_DYNA")
                resu.add_dire_response(option, spectre_dir, R_c, "VALE_QS")
                resu.add_dire_response(option, spectre_dir, R_prim, "VALE_INER")
                # save for INFO
                l_SA[spectre_dir] = S_r_freq
            # step 7 : reponse by directional combinaison
            # Get input COMB_DIRECTION
            comb_direction = args["COMB_DIRECTION"]
            R_xyz, R_newmark_all = comb_directions(comb_direction, l_R_x)
            # POST_ROCHE / part dynamique et pseudo statique
            R_d, Rd_newmark_all = comb_directions(comb_direction, l_part_d)
            R_ps, Rps_newmark_all = comb_directions(comb_direction, l_part_s)
            # RCCM
            R_prim, R_prim_newmark_all = comb_directions(comb_direction, l_R_prim)
            R_seco, R_seco_newmark_all = np.zeros(np.shape(R_prim)), np.zeros(
                np.shape(R_prim_newmark_all)
            )
            # Print output for combined response
            resu.add_response(option, R_xyz, R_newmark_all, "VALE_TOTA")
            resu.add_response(option, R_d, Rd_newmark_all, "VALE_DYNA")
            resu.add_response(option, R_ps, Rps_newmark_all, "VALE_QS")
            resu.add_response(option, R_prim, R_prim_newmark_all, "VALE_INER")
            # Print out for INFO = 1 or 2
            if verbosity and i_option == 1:
                # about mode_meca
                list_para = mode_meca.LIST_PARA()
                # shown_name
                show_name, show_type = _get_object_repr(mode_meca)
                # info for modal basis to be considered/combined
                for direction in ["OX", "OY", "OZ"]:
                    if "OX" in direction:
                        fact_partici = l_fact_partici[0]
                        masse_effe = l_masse_effe[0]
                    elif "OY" in direction:
                        fact_partici = l_fact_partici[1]
                        masse_effe = l_masse_effe[1]
                    elif "OZ" in direction:
                        fact_partici = l_fact_partici[2]
                        masse_effe = l_masse_effe[2]
                    UTMESS("I", "SEISME_48")
                # about spectra
                for i_dir in range(len(spectres[0])):
                    # Spectrum information
                    [
                        spectre_dir,
                        spectre_nappe,
                        spectre_coeff,
                        spectre_corr_freq,
                        spectre_nature,
                    ] = [spectres[i][i_dir] for i in range(5)]
                    # nature of spectra
                    UTMESS("I", "SEISME_17", valk=spectre_nature)
                    # info of read value on spectra
                    UTMESS("I", "SEISME_53")
                    for i_freq in range(len(freqs)):
                        dict_args = dict(
                            vali=nume_modes[i_freq],
                            valr=(freqs[i_freq], amors[i_freq], l_SA[spectre_dir][i_freq]),
                            valk=spectre_dir,
                        )
                        UTMESS("I", "SEISME_54", **dict_args)
                # about correction by pseudo-mode
                if mode_corr == "OUI":
                    for i_dir in range(len(spectres[0])):
                        [
                            spectre_dir,
                            spectre_nappe,
                            spectre_coeff,
                            spectre_corr_freq,
                            spectre_nature,
                        ] = [spectres[i][i_dir] for i in range(5)]
                        # cutting frequency et ZPA
                        UTMESS("I", "SEISME_56")
                        dict_args = dict(
                            valr=(l_pseudo[spectre_dir][0], l_pseudo[spectre_dir][1]),
                            valk=(spectre_dir, "mono_appui"),
                        )
                        UTMESS("I", "SEISME_57", **dict_args)
                # about combinaison of response due to DDS
                if comb_dds_correle:
                    UTMESS("I", "SEISME_19", valk=comb_dds_correle)
                # about directional combinaison
                UTMESS("I", "SEISME_18", valk=comb_direction)
            # end mono_appui
        elif type_analyse == "MULT_APPUI":
            # Run combinaison for mult_appui
            # step 1: Get eigen-vector by option to calculated
            if option not in ["VITE", "ACCE_ABSOLU"]:
                phis = get_phis(mode_meca, option, nume_ordres)
            elif option == "VITE" or option == "ACCE_ABSOLU":
                phis = get_phis(mode_meca, "DEPL", nume_ordres)
            else:
                raise Exception("OPTION '{option}' n'est pas pris en compte".format(option=option))
            # step 2: spectral value
            l_R_x = []  # list of directional total result
            l_part_d = []  # list of directional result part dyn
            l_part_s = []  # list of directional result part pseudo statique
            l_R_prim = []  # list of directional result RCCM part primaire
            l_R_seco = []  # list of directional result RCCM part secondaire
            # preparing for printing out to INFO
            l_SA = {}  # dict of info about read spectra at eigen-frequencies by direction
            l_pseudo = {}  # dict of info about correction by pseudo mode by direction
            l_parti = {}  # dict of info about participation factor in case of mult_appui
            # Iteration on directions
            for i_dir in range(max(len(dir_all), len(D_e_dirs_all))):
                direction = dir_all[i_dir]
                # ("iteration on direction:", direction)
                # save for print out as INFO
                # dict of info about read spectra at eigen-frequencies by direction
                l_SA[direction] = {}
                # dict of info about correction by pseudo mode by direction
                l_pseudo[direction] = {}
                # dict of info about participation factor in case of mult_appui
                l_parti[direction] = {}
                # Iteration on appuis
                nb_appui = len(nom_appuis_all)
                R_appuis = []  # list reponse tous appuis
                for i_appui in range(nb_appui):
                    nom_appui = spectres[i_appui][5]
                    # ("iteration on support:", nom_appui)
                    # Reponse oscillator et pseudo-mode
                    if direction in spectres[i_appui][0]:
                        # "check of direction %s found in all spectra of support %s"
                        index_dir = spectres[i_appui][0].index(direction)
                        # Extraction info at appui
                        nom_appui = spectres[i_appui][5]
                        [
                            spectre_dir,
                            spectre_nappe,
                            spectre_coeff,
                            spectre_corr_freq,
                            spectre_nature,
                        ] = [spectres[i_appui][i][index_dir] for i in range(5)]
                        # Pulsation before correction by corr_freq
                        w_r = 2 * np.pi * freqs
                        # Correction of spectrum by corr_freq
                        if spectre_corr_freq == "OUI":
                            # pulsation propre amortie
                            w_r *= np.sqrt(1 - amors**2)
                            correct = np.sqrt(1 - amors**2)
                        else:
                            w_r *= 1
                            correct = 1
                        # Spectrum interpolation
                        S_r_freq = []
                        for i_freq in range(len(freqs)):
                            S_r_freq.append(
                                spectre_nappe(amors[i_freq], freqs[i_freq]) * spectre_coeff
                            )
                        # correction by corr_freq
                        S_r_freq *= correct
                        # cutting frequency
                        if freq_coup_in is not None:
                            freq_coup = freq_coup_in
                        else:
                            freq_coup = freqs[-1]
                        # Correction of spectrum by nature
                        if spectre_nature == "DEPL":
                            S_r_freq *= (2 * np.pi * freqs) ** 2
                        elif spectre_nature == "VITE":
                            S_r_freq *= 2 * np.pi * freqs
                        elif spectre_nature == "ACCE":
                            S_r_freq *= 1
                        # iteration on node in appui
                        l_nodes_num = l_nodes_num_all[nom_appuis_all.index(nom_appui)]
                        l_nodes_name = l_nodes_name_all[nom_appuis_all.index(nom_appui)]
                        l_group_no = l_group_no_all[nom_appuis_all.index(nom_appui)]
                        R_m_noeuds = []  # stockage reponse oscillateur de tous noeuds
                        R_c_noeuds = []  # stockage reponse pseudo-mode de tous noeuds
                        # save for INFO
                        l_SA[spectre_dir][nom_appui] = S_r_freq
                        # modal spectral response
                        # preparing a list of reac_noda for all modes and all group_no
                        reac_noda = []
                        dofs = (
                            mode_meca.getDOFNumbering()
                            .getEquationNumbering()
                            .getDOFs(list_cmp=[direction.replace("O", "D")], list_grpno=l_group_no)
                        )
                        for imode in nume_ordres:
                            # all values of reac_node for mode
                            reac_noda_all = mode_meca.getField("REAC_NODA", imode)
                            # values for all nodes in the same group_no without knowing the order
                            reac_noda_by_mode = reac_noda_all.getValues(dofs)
                            reac_noda.append(reac_noda_by_mode)
                        # Reac_node for all modes at 1 support
                        reac_noda = np.sum(reac_noda, axis=1)
                        # participation factor for all modes at 1 support
                        fact_partici = -1.0 * reac_noda / (gene_masses * w_r**2)
                        # Spectral response at node, mode, direction
                        if option not in ["VITE", "ACCE_ABSOLU"]:
                            R_mi_all = (S_r_freq * fact_partici / w_r**2)[:, None] * phis
                            pr_wr2_phi_c_all = (fact_partici / (2 * np.pi * freqs) ** 2)[
                                :, None
                            ] * phis
                        elif option == "VITE":  # ici: phis correspond à DEPL
                            R_mi_all = (S_r_freq * fact_partici / w_r)[:, None] * phis
                            pr_wr2_phi_c_all = (fact_partici / (2 * np.pi * freqs))[:, None] * phis
                        elif option == "ACCE_ABSOLU":  # ici: phis correspond à DEPL
                            R_mi_all = (S_r_freq * fact_partici)[:, None] * phis
                            pr_wr2_phi_c_all = (fact_partici)[:, None] * phis
                        # Check if cutting frequency is smaller than firt mode
                        if freq_coup is not None and freq_coup >= freqs[0]:
                            R_mi = R_mi_all
                            pr_wr2_phi_c = pr_wr2_phi_c_all
                        elif freq_coup < freqs[0]:
                            R_mi = np.zeros(np.shape(phis))
                            pr_wr2_phi_c = np.zeros(np.shape(phis))
                            # Raise alarm message if first mode bigger than cutting frequency
                            UTMESS("A", "SEISME_96", valr=freq_coup)
                        # saving all responses at all nodes of considerd appui
                        R_m_noeuds.append(R_mi)
                        # --- pseudo mode response
                        if mode_corr == "OUI":
                            # check if freq_coup is present
                            if freq_coup_in is not None:
                                UTMESS("A", "SEISME_95", valr=freq_coup)
                            # correction by pseudo-mode
                            R_c_noeud, S_r_freq_coup = corr_pseudo_mode_mult(
                                option,
                                pseudo_mode,
                                l_group_no,
                                mesh,
                                direction,
                                amors,
                                freq_coup,
                                pr_wr2_phi_c,
                                w_r,
                                spectre_nappe,
                                spectre_coeff,
                                spectre_corr_freq,
                                spectre_nature,
                            )
                            # save for INFO
                            l_pseudo[spectre_dir][nom_appui] = [freq_coup, S_r_freq_coup]
                        else:
                            R_c_noeud = np.zeros(np.shape(phis[0]))
                        # saving all pseudo-modes at all nodes of considerd appui
                        R_c_noeuds.append(R_c_noeud)
                        # step 3: LINE combinaison intra-appui: considered as CORRELATED
                        R_m_appui = np.sum(R_m_noeuds, axis=0)
                        R_c_appui = np.sum(R_c_noeuds, axis=0)
                    else:  # if direction is not found in SPECTRE then nul
                        R_m_appui = np.zeros((len(w_r), len(phis[0])))
                        R_c_appui = np.zeros(np.shape(phis[0]))
                    # save for INFO
                    l_parti[spectre_dir][nom_appui] = fact_partici
                    # Print out spectral response at each mode, direction and appui
                    # VALE_SPEC
                    resu.add_spectral_response(option, direction, R_m_appui, nume_ordres, nom_appui)
                    # response due to DDS
                    if (
                        depl_mult_appui is not None
                        and nom_appui in depl_mult_appuis[3]
                        and direction in depl_mult_appuis[3][nom_appui]
                    ):
                        # Get information from DEPL_MULT_APPUI
                        mode_stat = depl_mult_appuis[0][nom_appui]
                        # iteration on node in APPUI
                        l_nodes_num = l_nodes_num_all[nom_appuis_all.index(nom_appui)]
                        l_nodes_name = l_nodes_name_all[nom_appuis_all.index(nom_appui)]
                        R_e_noeuds = []  # stock of DDS response at all nodes in considered APPUI
                        for i_node in range(len(l_nodes_num)):
                            # ("Iteration sur les nodes dans un appui")
                            node_name = l_nodes_name[i_node]
                            noeud_cmp = node_name.ljust(8) + direction.replace("O", "D")
                            stat_nume_modes = mode_stat.getAccessParameters()["NUME_MODE"]
                            stat_noeud_cmps = mode_stat.getAccessParameters()["NOEUD_CMP"]
                            if noeud_cmp in stat_noeud_cmps:
                                stat_nume_mode = stat_nume_modes[stat_noeud_cmps.index(noeud_cmp)]
                                # Get static mode
                                if option in ["DEPL", "REAC_NODA", "FORC_NODA"]:
                                    phi_stat = mode_stat.getField(
                                        option, stat_nume_mode
                                    ).getValues()
                                    # Reponse due to DDS at node
                                    R_e_noeud = (
                                        np.array(phi_stat)
                                        * depl_mult_appuis[2][nom_appui][direction]
                                    )
                                elif option in [
                                    "EFGE_ELNO",
                                    "EGRU_ELNO",
                                    "SIEF_ELGA",
                                    "SIGM_ELNO",
                                    "SIPO_ELNO",
                                    "SIEF_ELNO",
                                ]:
                                    phi_stat = mode_stat.getField(
                                        option, stat_nume_mode
                                    ).getValues()
                                    # Reponse due to DDS at node
                                    R_e_noeud = (
                                        np.array(phi_stat)
                                        * depl_mult_appuis[2][nom_appui][direction]
                                    )

                                elif option in ["VITE", "ACCE_ABSOLU"]:
                                    raise Exception(
                                        "Il est inutile de calculer la réponse d'entrainement pour {option}".format(
                                            option=option
                                        )
                                    )
                            else:  # noued_cmp is not found --> raise fatal erreur
                                raise Exception(
                                    "NOUED_CMP: {0} n'existe pas dans la base mode_statique".format(
                                        noeud_cmp
                                    )
                                )
                            # save responses DDS at all nodes in considered APPUI
                            R_e_noeuds.append(R_e_noeud)
                        # save en array
                        R_e_noeuds = np.array(R_e_noeuds)
                        # calculate response at APPUI
                        # step 3: LINE combinaison intra-appui: considered as CORRELATED
                        R_e_appui = np.sum(R_e_noeuds, axis=0)
                    else:
                        R_e_appui = np.zeros(np.shape(phis[0]))
                    # all responses at one APPUI
                    R_appuis.append([nom_appui, R_m_appui, R_e_appui, R_c_appui])
                # step 4: iteration on group_appui
                l_R_x_j = []  # list des reponse directionnelle de tous les group_appui
                l_part_d_j = []  # list des reponse part dyn  de tous les group_appui
                l_part_s_j = []  # list des reponse part sta  de tous les group_appui
                l_R_prim_j = []  # list des reponse RCCM part primaire de tous les group_appui
                l_R_seco_j = []  # list des reponse RCCM part primaire de tous les group_appui
                for j_group_appui in range(len(group_appui_all)):
                    # ("iteration on nom_group_appui:", nom_group_appui)
                    appui_in_groupe = group_appui_all[j_group_appui]
                    # ("nom_group_appui contient appui_in_groupe:", appui_in_groupe)
                    R_m_j, R_e_j, R_c_j = [], [], []
                    for [appui, R_m, R_e, R_c] in R_appuis:
                        if appui in appui_in_groupe:
                            R_m_j.append(R_m), R_e_j.append(R_e), R_c_j.append(R_c)
                    R_m_j = np.array(R_m_j)
                    R_e_j = np.array(R_e_j)
                    R_c_j = np.array(R_c_j)
                    # step 5: combinaison of all appuis in a group_appui
                    # rules are different for different components
                    # pour oscillator and pseudo-mode: CUMUL_INTRA
                    # pour dds: rule definied in COMB_DDS_CORRELE
                    # response of oscillator by group_appui

                    R_m_group_appui = comb_appui_corr(cumul_intra, R_m_j)
                    # response of pseudo-mode by group_appui
                    R_c_group_appui = comb_appui_corr(cumul_intra, R_c_j)
                    # response of DDS by group_appui)
                    R_e_group_appui = comb_appui_corr(comb_dds_correle, R_e_j)
                    # step 6: modal combinaison pour R_m_group_appui
                    R_m2, R_qs = comb_modal_response.get(R_m_group_appui)
                    # Automatic correction for ACCE_ABSOLU: not used in mult-appui
                    if option == "ACCE_ABSOLU":
                        # raise fatal error message to stop
                        UTMESS("F", "SEISME_10", valk=option)
                        # participation factor
                        if "OX" in spectre_dir:
                            fact_partici = l_fact_partici[0]
                            components = "DX"
                        elif "OY" in spectre_dir:
                            fact_partici = l_fact_partici[1]
                            components = "DY"
                        elif "OZ" in spectre_dir:
                            fact_partici = l_fact_partici[2]
                            components = "DZ"
                        # unit field
                        acce_unitaire = mode_meca.getField("DEPL", 1).copy()
                        acce_unitaire.setValues({components: 1.0}, [])
                        # recalculate pr_wr2_phi for acce_absolu
                        pr_wr2_phi_all = (fact_partici)[:, None] * phis
                        pr_wr2_phi = pr_wr2_phi_all[freqs <= freq_coup]
                        # i_appui=0 : only first support to be considered
                        index_dir = spectres[0][0].index(direction)
                        [
                            spectre_dir,
                            spectre_nappe,
                            spectre_coeff,
                            spectre_corr_freq,
                            spectre_nature,
                        ] = [spectres[0][i][index_dir] for i in range(5)]
                        S_r_freq_coup = spectre_nappe(amors[-1], freq_coup) * spectre_coeff
                        R_tt = (
                            acce_unitaire.getValues() - np.sum(pr_wr2_phi, axis=0)
                        ) * S_r_freq_coup
                        # reponse oscillator by adding correction for ACCE_ABSOLU:Not used
                        R_m2 += R_tt**2
                    # reponse oscillator
                    R_m = np.sqrt(R_m2)
                    # step 7 : reponse by direction for group_appui
                    # ("reponse directionnelle : sqrt(Rm**2 + Rc**2 + Re**2)")
                    R_x_j = np.sqrt(R_m2 + (R_qs + R_c_group_appui) ** 2 + R_e_group_appui**2)
                    l_R_x_j.append(R_x_j)
                    # RCCM part primaire
                    R_prim_group_appui = np.sqrt(R_m2 + (R_qs + R_c_group_appui) ** 2)
                    # POST_ROCHE/ part dynamique et pseudo statique
                    l_part_d_j.append(R_m)
                    l_part_s_j.append(R_c_group_appui)
                    # RCCM
                    l_R_prim_j.append(R_prim_group_appui)
                    l_R_seco_j.append(R_e_group_appui)
                # step 8: combinaison of all group_appui
                # While nb of group_appui > 1: rule = QUAD (considered as DECORRELATED)
                if len(l_R_x_j) > 1:
                    # combi all group_appui
                    if cumul_inter == "QUAD":
                        R_x = np.sqrt(np.sum(np.array(l_R_x_j) ** 2, axis=0))
                        part_d_x = np.sqrt(np.sum(np.array(l_part_d_j) ** 2, axis=0))
                        part_s_x = np.sqrt(np.sum(np.array(l_part_s_j) ** 2, axis=0))
                        R_prim_x = np.sqrt(np.sum(np.array(l_R_prim_j) ** 2, axis=0))
                        R_seco_x = np.sqrt(np.sum(np.array(l_R_seco_j) ** 2, axis=0))
                    elif cumul_inter == "ABS":  # HB: NEW METHOD
                        R_x = np.sum(np.abs(np.array(l_R_x_j)), axis=0)
                        part_d_x = np.sum(np.abs(np.array(l_part_d_j)), axis=0)
                        part_s_x = np.sum(np.abs(np.array(l_part_s_j)), axis=0)
                        R_prim_x = np.sum(np.abs(np.array(l_R_prim_j)), axis=0)
                        R_seco_x = np.sum(np.abs(np.array(l_R_seco_j)), axis=0)
                    elif cumul_inter == "LINE":  # HB: NEW METHOD
                        R_x = np.sum(np.array(l_R_x_j), axis=0)
                        part_d_x = np.sum(np.array(l_part_d_j), axis=0)
                        part_s_x = np.sum(np.array(l_part_s_j), axis=0)
                        R_prim_x = np.sum(np.array(l_R_prim_j), axis=0)
                        R_seco_x = np.sum(np.array(l_R_seco_j), axis=0)

                elif len(l_R_x_j) == 1:  # un seul group appui = multi appui correle
                    # combi group_appui not done if one group_appui
                    R_x = l_R_x_j[0]
                    part_d_x = l_part_d_j[0]
                    part_s_x = l_part_s_j[0]
                    R_prim_x = l_R_prim_j[0]
                    R_seco_x = l_R_seco_j[0]
                l_R_x.append(R_x)
                # POST_ROCHE / part dynamique et pseudo statique
                l_part_d.append(part_d_x)
                l_part_s.append(part_s_x)
                # RCCM part primaire
                l_R_prim.append(R_prim_x)
                l_R_seco.append(R_seco_x)
                # Print out reponse of each direction
                resu.add_dire_response(option, direction, R_x, "VALE_DIRE")
                resu.add_dire_response(option, direction, part_d_x, "VALE_DYNA")
                resu.add_dire_response(option, direction, part_s_x, "VALE_QS")
                resu.add_dire_response(option, direction, R_seco_x, "VALE_DDS")
                resu.add_dire_response(option, direction, R_prim_x, "VALE_INER")
            # step 9 : combinaison of all directions
            # Get input COMB_DIRECTION
            comb_direction = args["COMB_DIRECTION"]
            R_xyz, R_newmark_all = comb_directions(comb_direction, l_R_x)
            # POST_ROCHE/ part dynamique et pseudo statique
            R_d, Rd_newmark_all = comb_directions(comb_direction, l_part_d)
            R_ps, Rps_newmark_all = comb_directions(comb_direction, l_part_s)
            # RCCM
            R_prim, R_prim_newmark_all = comb_directions(comb_direction, l_R_prim)
            R_seco, R_seco_newmark_all = comb_directions(comb_direction, l_R_seco)
            # ----------------------------------------------------------
            # Print out total response
            resu.add_response(option, R_xyz, R_newmark_all, "VALE_TOTA")
            resu.add_response(option, R_d, Rd_newmark_all, "VALE_DYNA")
            resu.add_response(option, R_ps, Rps_newmark_all, "VALE_QS")
            resu.add_response(option, R_prim, R_prim_newmark_all, "VALE_INER")
            resu.add_response(option, R_seco, R_seco_newmark_all, "VALE_DDS")
            if verbosity and i_option == 1:
                # about mode_meca
                list_para = mode_meca.LIST_PARA()
                # shown_name
                show_name, show_type = _get_object_repr(mode_meca)
                # info for modal basis to be considered/combined
                for direction in ["OX", "OY", "OZ"]:
                    if "OX" in direction:
                        fact_partici = l_fact_partici[0]
                        masse_effe = l_masse_effe[0]
                    elif "OY" in direction:
                        fact_partici = l_fact_partici[1]
                        masse_effe = l_masse_effe[1]
                    elif "OZ" in direction:
                        fact_partici = l_fact_partici[2]
                        masse_effe = l_masse_effe[2]
                    UTMESS("I", "SEISME_48")
                # about spectra
                for i_dir in range(max(len(dir_all), len(D_e_dirs_all))):
                    direction = dir_all[i_dir]
                    for i_appui in range(nb_appui):
                        nom_appui = spectres[i_appui][5]
                        if direction in spectres[i_appui][0]:
                            index_dir = spectres[i_appui][0].index(direction)
                            # Extraction info at appui
                            [
                                spectre_dir,
                                spectre_nappe,
                                spectre_coeff,
                                spectre_corr_freq,
                                spectre_nature,
                            ] = [spectres[i_appui][i][index_dir] for i in range(5)]
                            # nature of spectra
                            UTMESS("I", "SEISME_17", valk=spectre_nature)
                            # info of read value on spectra
                            UTMESS("I", "SEISME_97")
                            for i_freq in range(len(freqs)):
                                valr = (
                                    freqs[i_freq],
                                    amors[i_freq],
                                    l_SA[spectre_dir][nom_appui][i_freq],
                                    l_parti[spectre_dir][nom_appui][i_freq],
                                )
                                karg = dict(
                                    vali=nume_modes[i_freq],
                                    valr=valr,
                                    valk=(spectre_dir, nom_appui),
                                )
                                UTMESS("I", "SEISME_98", **karg)
                            # about correction by pseudo-mode
                            if mode_corr == "OUI":
                                # cutting frequency et ZPA
                                UTMESS("I", "SEISME_56")
                                dict_args = dict(
                                    valr=(
                                        l_pseudo[spectre_dir][nom_appui][0],
                                        l_pseudo[spectre_dir][nom_appui][1],
                                    ),
                                    valk=(spectre_dir, nom_appui),
                                )
                                UTMESS("I", "SEISME_57", **dict_args)
                # about combinaison of response due to DDS
                if comb_dds_correle:
                    UTMESS("I", "SEISME_19", valk=comb_dds_correle)
                    # about directional combinaison
                    UTMESS("I", "SEISME_18", valk=comb_direction)
        # end mult_appui
        if type_analyse == "ENVELOPPE":
            # step 1: Get eigen-vector to combine from mode_meca
            if option not in ["VITE", "ACCE_ABSOLU"]:
                phis = get_phis(mode_meca, option, nume_ordres)
            elif option == "VITE" or option == "ACCE_ABSOLU":
                phis = get_phis(mode_meca, "DEPL", nume_ordres)
            else:
                raise Exception("OPTION '{option}' n'est pas pris en compte".format(option=option))
            # step 2: spectral value
            l_R_x = []  # list of directional total result
            l_part_d = []  # list of directional result part dynamique
            l_part_s = []  # list of directional result part pseudo-statique
            l_R_prim = []  # list of RCCM part primaire
            # preparing for printing out to INFO
            l_SA = {}  # list of info about read spectra at eigen-frequencies by direction
            l_pseudo = {}  # list of info about correction by pseudo mode by direction
            # ----------------------------------------
            # Preparing for type analysis = enveloppe :
            l_axis = []
            l_corr_freq = []
            l_spec_nat = []
            for i_spec in range(len(spectres)):
                l_axis += spectres[i_spec][0]
                l_corr_freq += spectres[i_spec][3]
                l_spec_nat += spectres[i_spec][4]
            axes_retenu = sorted(set(l_axis))
            is_corr_freq = sorted(set(l_corr_freq))
            spec_nat = sorted(set(l_spec_nat))  # spectre nature

            corr_freq_overall = "NON"
            # correction in frequencies is taken into account if correction is required for all supports:
            if len(is_corr_freq) == 1:
                corr_freq_overall = is_corr_freq[0]

            # iteration on direction
            for i_dir in range(len(axes_retenu)):
                # eigen-pulsation before corr_freq
                w_r = 2 * np.pi * freqs
                # ------------------------------------
                # ENVELOPE : Selection of maximal value of spectrum :
                val_spec_allgr = []
                axis = axes_retenu[i_dir]
                for i_spec, spec in enumerate(spectres):
                    if axis in spec[0]:
                        i_axis = spec[0].index(axis)
                        # corr = np.ones(len(freqs))  # = unité - correction for frequencies
                        # if spec[3][i_axis] == 'OUI':
                        #     corr = np.sqrt(1 - amors**2)
                        val_spec = []  # value of spectre for each support group
                        for i_freq, freq in enumerate(freqs):
                            val_spec.append(spec[1][i_axis](amors[i_freq], freq) * spec[2][i_axis])
                        val_spec_allgr.append(val_spec)
                val_spec_allgr = np.array(val_spec_allgr)
                S_r_freq = np.max(val_spec_allgr, axis=0)
                # If correction for frequencies for all axis :
                if corr_freq_overall == "OUI":
                    correct = np.sqrt(1 - amors**2)
                else:
                    correct = 1

                # # Spectrum interpolation
                # [spectre_dir, spectre_nappe, spectre_coeff, spectre_corr_freq, spectre_nature] = [
                #     spectres[i][i_dir] for i in range(5)
                # ]
                # # Correction for frequency by corr_freq
                # if spectre_corr_freq == "OUI":
                #     correct = np.sqrt(1 - amors**2)
                # else:
                #     correct = 1
                # Pulsation afeter corr_freq
                w_r *= correct
                # # Spectrale values at eigen-frequencies
                # S_r_freq = []
                # for i_freq in range(len(freqs)):
                #     S_r_freq.append(spectre_nappe(amors[i_freq], freqs[i_freq]) * spectre_coeff)
                # Correction by corr_freq for spectrum
                S_r_freq *= correct
                # Cutting frequency
                if freq_coup_in is not None:
                    freq_coup = freq_coup_in
                else:
                    freq_coup = freqs[-1]
                #  Correction of spectrum by nature of spectrum
                spectre_nature = spec_nat[0]
                if spectre_nature is "DEPL":
                    S_r_freq = ((2 * np.pi * freqs) ** 2) * S_r_freq
                elif spectre_nature is "VITE":
                    S_r_freq = 2 * np.pi * freqs * S_r_freq
                elif spectre_nature is "ACCE":
                    S_r_freq = S_r_freq
                # Participation factor by direction
                spectre_dir = axes_retenu[i_dir]
                if "OX" in spectre_dir:
                    fact_partici = l_fact_partici[0]
                    components = "DX"
                elif "OY" in spectre_dir:
                    fact_partici = l_fact_partici[1]
                    components = "DY"
                elif "OZ" in spectre_dir:
                    fact_partici = l_fact_partici[2]
                    components = "DZ"
                else:
                    raise Exception(
                        "Direction '{spectre_dir}' ne correspond pas à aucun axe global".format(
                            spectre_dir=spectre_dir
                        )
                    )
                # Spectral response
                if option not in ["VITE", "ACCE_ABSOLU"]:
                    R_mi_all = (S_r_freq * fact_partici / w_r**2)[:, None] * phis
                    pr_wr2_phi_all = (fact_partici / w_r**2)[:, None] * phis
                    pr_wr2_phi_c_all = (fact_partici / (2 * np.pi * freqs) ** 2)[
                        :, None
                    ] * phis  # interrogration ??? pq ne pas utiliser omega corrige?
                elif option == "VITE":  # ici: phis correspond à DEPL
                    R_mi_all = (S_r_freq * fact_partici / w_r)[:, None] * phis
                    pr_wr2_phi_all = (fact_partici / w_r)[:, None] * phis
                    pr_wr2_phi_c_all = (fact_partici / w_r)[:, None] * phis
                elif option == "ACCE_ABSOLU":  # ici: phis correspond à DEPL
                    R_mi_all = (S_r_freq * fact_partici)[:, None] * phis
                    pr_wr2_phi_all = (fact_partici)[:, None] * phis
                    pr_wr2_phi_c_all = (fact_partici)[:, None] * phis
                # in case where the first mode is bigger than cutting frequency
                if freq_coup is not None and freq_coup >= freqs[0]:
                    R_mi = R_mi_all
                    pr_wr2_phi = pr_wr2_phi_all
                    pr_wr2_phi_c = pr_wr2_phi_c_all
                elif freq_coup < freqs[0]:
                    R_mi = np.zeros(np.shape(phis))
                    pr_wr2_phi = np.zeros(np.shape(phis))
                    pr_wr2_phi_c = np.zeros(np.shape(phis))
                    # Raise alarm for zero mode to be considered before cutting frequency
                    UTMESS("A", "SEISME_96", valr=freq_coup)
                # Print output for spectral value for each mode and direction
                resu.add_spectral_response(option, spectre_dir, R_mi_all, nume_ordres)
                # step 3: modal combinaison
                # Get input COMB_MODE
                R_m2, R_qs = comb_modal_response.get(R_mi)
                # Automatic correction for ACCE_ABSOLU in mono-appui
                if option == "ACCE_ABSOLU":
                    # field of unit value for acce_absolu
                    acce_unitaire = mode_meca.getField("DEPL", 1).copy()
                    acce_unitaire.setValues({components: 1.0}, [])
                    # S_r_freq_coup = spectre_nappe(amors[-1], freq_coup) * spectre_coeff
                    # --------------------------------------------------------------------
                    # ENVELOPE :
                    # Searching for the envelope of the ZPA value at the cutting frequency
                    spec_freq_coup_allgr = []
                    for i_spec, spec in enumerate(spectres):
                        if axis in spec[0]:
                            i_axis = spec[0].index(axis)
                            spec_freq_coup_allgr.append(
                                spec[1][i_axis](amors[-1], freq_coup) * spec[2][i_axis]
                            )
                    # maximal values of ZPA at the cutting frequency
                    S_r_freq_coup = max(spec_freq_coup_allgr)

                    R_tt = (acce_unitaire.getValues() - np.sum(pr_wr2_phi, axis=0)) * S_r_freq_coup
                    # add to combined modale responses in square
                    R_m2 += R_tt**2
                # modale response
                R_m = np.sqrt(R_m2)
                # step 4 : Entrainement zero pour mon_appui
                R_e2 = np.zeros(np.shape(R_m2))
                # step 5 : pseudo-mode response
                if mode_corr == "OUI":
                    # check if cutting frequency is present
                    if freq_coup_in is None:
                        UTMESS("A", "SEISME_95", valr=freq_coup)

                    # Searching for the envelope of the ZPA value at the cutting frequency
                    spec_freq_coup_allgr = []
                    for i_spec, spec in enumerate(spectres):
                        if axis in spec[0]:
                            i_axis = spec[0].index(axis)
                            spec_freq_coup_allgr.append(
                                spec[1][i_axis](amors[-1], freq_coup) * spec[2][i_axis]
                            )
                    # maximal values of ZPA at the cutting frequency
                    S_r_freq_coup = max(spec_freq_coup_allgr)
                    spectre_corr_freq = corr_freq_overall  # for all supports
                    R_c = corr_pseudo_mode_enveloppe(
                        option,
                        pseudo_mode,
                        amors,
                        freq_coup,
                        pr_wr2_phi_c,
                        w_r,
                        spectre_dir,
                        S_r_freq_coup,
                        spectre_corr_freq,
                        spectre_nature,
                    )
                    # # calculate pseudo-mode
                    # R_c, S_r_freq_coup = corr_pseudo_mode_mono(
                    #     option,
                    #     pseudo_mode,
                    #     amors,
                    #     freq_coup,
                    #     pr_wr2_phi_c,
                    #     w_r,
                    #     spectre_dir,
                    #     spectre_nappe,
                    #     spectre_coeff,
                    #     spectre_corr_freq,
                    #     spectre_nature,
                    # )
                    # save for INFO
                    l_pseudo[spectre_dir] = [freq_coup, S_r_freq_coup]
                else:
                    R_c = np.zeros(np.shape(R_m2))
                    S_r_freq_coup = None
                # step 6 : reponse by direction
                # total
                R_x = np.sqrt(R_m2 + (R_qs + R_c) ** 2 + R_e2)
                # inertial part (part primaire)
                R_prim = np.sqrt(R_m2 + (R_qs + R_c) ** 2)
                # add total directionnal responses
                l_R_x.append(R_x)
                # POST_ROCHE/ part dynamique et pseudo statique
                l_part_d.append(R_m)
                l_part_s.append(R_c)
                # RCCM part primaire
                l_R_prim.append(R_prim)
                # Print out response for mono_appui by direction
                resu.add_dire_response(option, spectre_dir, R_x, "VALE_DIRE")
                resu.add_dire_response(option, spectre_dir, R_m, "VALE_DYNA" )
                resu.add_dire_response(option, spectre_dir, R_c, "VALE_QS")
                resu.add_dire_response(option, axis, R_prim, "VALE_INER")
                # VALE_DDS
                if any(type_resu[i].get("TYPE") == "VALE_DDS" for i in range(len(type_resu))):
                    print("Pas de VALE_DDS calculée pour la méthode ENVELOPPE")
                # save for INFO
                l_SA[spectre_dir] = S_r_freq
            # step 7 : reponse by directional combinaison
            # Get input COMB_DIRECTION
            comb_direction = args["COMB_DIRECTION"]
            R_xyz, R_newmark_all = comb_directions(comb_direction, l_R_x)
            # POST_ROCHE / part dynamique et pseudo statique
            R_d, Rd_newmark_all = comb_directions(comb_direction, l_part_d)
            R_ps, Rps_newmark_all = comb_directions(comb_direction, l_part_s)
            # RCCM
            R_prim, R_prim_newmark_all = comb_directions(comb_direction, l_R_prim)
            R_seco, R_seco_newmark_all = np.zeros(np.shape(R_prim)), np.zeros(
                np.shape(R_prim_newmark_all)
            )
            # Print output for combined response
            resu.add_response(option, R_xyz, R_newmark_all, "VALE_TOTA")
            resu.add_response(option, R_d, Rd_newmark_all, "VALE_DYNA")
            resu.add_response(option, R_ps, Rps_newmark_all,"VALE_QS")
            resu.add_response(option, R_prim, R_prim_newmark_all, "VALE_INER")
            # Print out for INFO = 1 or 2
            if verbosity and i_option == 1:
                # about mode_meca
                list_para = mode_meca.LIST_PARA()
                # shown_name
                show_name, show_type = _get_object_repr(mode_meca)
                # info for modal basis to be considered/combined
                for direction in ["OX", "OY", "OZ"]:
                    if "OX" in direction:
                        fact_partici = l_fact_partici[0]
                        masse_effe = l_masse_effe[0]
                        masse_effe_un = l_masse_effe_un[0]
                    elif "OY" in direction:
                        fact_partici = l_fact_partici[1]
                        masse_effe = l_masse_effe[1]
                        masse_effe_un = l_masse_effe_un[1]
                    elif "OZ" in direction:
                        fact_partici = l_fact_partici[2]
                        masse_effe = l_masse_effe[2]
                        masse_effe_un = l_masse_effe_un[2]
                    UTMESS("I", "SEISME_48")
                # about spectra
                for i_dir in range(len(axes_retenu)):
                    # Spectrum information
                    # [
                    #     spectre_dir,
                    #     spectre_nappe,
                    #     spectre_coeff,
                    #     spectre_corr_freq,
                    #     spectre_nature,
                    # ] = [spectres[i][i_dir] for i in range(5)]
                    spectre_dir = axes_retenu[i_dir]
                    # nature of spectra
                    UTMESS("I", "SEISME_17", valk=spectre_nature)
                    # info of read value on spectra
                    UTMESS("I", "SEISME_53")
                    for i_freq in range(len(freqs)):
                        vali = nume_modes[i_freq]
                        valr = (freqs[i_freq], amors[i_freq], l_SA[spectre_dir][i_freq])
                        UTMESS("I", "SEISME_54", vali=vali, valr=valr, valk=spectre_dir)
                # about correction by pseudo-mode
                if mode_corr == "OUI":
                    # for i_dir in range(len(spectres[0])):
                    #     [
                    #         spectre_dir,
                    #         spectre_nappe,
                    #         spectre_coeff,
                    #         spectre_corr_freq,
                    #         spectre_nature,
                    #     ] = [spectres[i][i_dir] for i in range(5)]
                    spectre_dir = axes_retenu[i_dir]
                    # cutting frequency et ZPA
                    UTMESS("I", "SEISME_56")
                    valr = (l_pseudo[spectre_dir][0], l_pseudo[spectre_dir][1])
                    UTMESS("I", "SEISME_57", valr=valr, valk=(spectre_dir, "ENVELOPPE"))
                # about combinaison of response due to DDS
                if comb_dds_correle:
                    UTMESS("I", "SEISME_19", valk=comb_dds_correle)
                # about directional combinaison
                UTMESS("I", "SEISME_18", valk=comb_direction)
            # end mono_appui
    # end
    return resu.get()
