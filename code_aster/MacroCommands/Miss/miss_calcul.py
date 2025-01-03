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

# person_in_charge: mathieu.courtois at edf.fr

"""Module permettant le lancement d'un calcul MISS3D

Les classes définies sont :
    CalculMiss      : classe générique d'un appel à MISS3D
    CalculMissImpe  : spécialisation au calcul d'impédance
    CalculMissFichierTemps : calcul dans le domaine de Laplace
    CalculMissIssf  : calcul avec prise en compte de l'interaction
                      sol-structure-fluide
    CalculMissPost  : étapes limitées au post-traitement
"""

import os
import os.path as osp
import re
import shutil
import socket
import traceback
from math import cos, pi, sin

import aster
import numpy as NP

from ...Cata.Syntax import _F
from ...CodeCommands import IMPR_MACR_ELEM, MACR_ELEM_DYNA
from ...Helpers.LogicalUnit import LogicalUnitFile
from ...Messages import UTMESS
from ...Utilities import ExecutionParameter, Timer, config
from ...Utilities.misc import _print, decode_str, encode_str, send_file, set_debug
from ...Utilities.mpi_utils import MPI
from ...Utilities.System import ExecCommand
from .miss_fichier_cmde import MissCmdeGen
from .miss_fichier_interf import fichier_chp, fichier_ext, fichier_mvol, fichier_sign
from .miss_fichier_option import fichier_option
from .miss_fichier_sol import fichier_sol
from .miss_post import PostMissFactory
from .miss_resu_aster import ResuAsterReader
from .miss_utils import copie_fichier

MODE = "wb"


class CalculMiss:

    """Définition d'un calcul MISS3D."""

    option_calcul = None

    @staticmethod
    def factory(parent, param):
        """Factory that returns the CalculMiss object"""
        if param["TYPE_RESU"] == "FICHIER_TEMPS":
            return CalculMissFichierTemps(parent, param)
        elif param["TYPE_RESU"] != "FICHIER" and not param["_exec_Miss"]:
            return CalculMissPost(parent, param)
        elif param["ISSF"] == "OUI":
            return CalculMissIssf(parent, param)
        else:
            return CalculMissImpe(parent, param)

    def __init__(self, parent, parameters):
        """Initialisations"""
        if not self.option_calcul:
            raise NotImplementedError("option_calcul non défini")
        self.parent = parent
        self.param = parameters
        self.data = None
        self.verbose = parameters["INFO"] >= 2
        self.debug = self.verbose
        self.resu_aster_reader = None
        if self.verbose:
            _print("Paramètres du calcul", self.param)
        if self.debug:
            self.timer = Timer()
            set_debug(True)

    def run(self):
        """Enchaine les tâches élémentaires"""
        self.prepare_donnees()
        self.execute()
        return self.post_traitement()

    def init_reader(self):
        """Initialise le lecteur de fichier Aster"""
        mcgrp = []
        for grp in self.param:
            if grp.startswith("GROUP_MA_") and self.param[grp] is not None:
                mcgrp.append(grp)
        self.resu_aster_reader = ResuAsterReader(len(mcgrp))

    def prepare_donnees(self):
        """Préparation des données"""
        self.cree_reptrav()
        self.init_reader()
        self.cree_resultat_aster()
        self.cree_fichier_sol()
        self.cree_fichier_mvol()
        self.cree_fichier_pc()
        self.cree_fichier_chp()
        self.cree_commande_miss()
        self.cree_fichier_option()
        # libérer la structure contenant les données numériques
        self.init_data()

    def init_data(self):
        """Libérer la structure contenant les données numériques"""
        self.data = None

    def execute(self):
        """Exécute MISS3D."""
        self._dbg_trace("Start")

        copie_fichier(self._fichier_tmp("in"), osp.join(self.param["_WRKDIR"], "MISS.IN"))
        cmd = ExecutionParameter().get_option("prog:run_miss3d") + " " + self.param["VERSION"]
        iret = 4
        try:
            os.chdir(self.param["_WRKDIR"])
            if self.verbose:
                _print("Exécution de MISS3D dans :", self.param["_WRKDIR"])
                os.system("ls -la")
            comment = "Lancement de la commande :\n%s" % cmd
            if self.verbose:
                aster.affiche("MESSAGE", comment)
            iret, output, error = ExecCommand(
                cmd, alt_comment=comment, verbose=False, separated_stderr=True
            )
            if self.verbose:
                _print("Contenu du répertoire après l'exécution de MISS3D :")
                os.system("ls -la")
        finally:
            os.chdir(self.param["_INIDIR"])
        UTMESS("I", "EXECLOGICIEL0_9", valk=output)
        miss_out = ""
        if osp.exists(self._fichier_tmp("OUT")):
            with open(self._fichier_tmp("OUT"), "r") as f:
                miss_out = f.read()
        is_ok = iret == 0 and (miss_out.find("INSUFFISAN") < 0 and miss_out.find("*** STOP") < 0)
        if not is_ok or self.verbose:
            aster.affiche("MESSAGE", miss_out)
        if not is_ok:
            UTMESS("I", "EXECLOGICIEL0_10", valk=error, print_as="E")
            UTMESS("F", "EXECLOGICIEL0_3", vali=[0, iret])
        self._dbg_trace("Stop")

    def post_traitement(self):
        """Opérations de post-traitement."""
        self._dbg_trace("Start")
        self.fichier_resultat()

        post = PostMissFactory(self.param["TYPE_RESU"], self.parent, self.param)
        post.set_filename_callback(self._fichier_tmp)
        toReturn = post.run()
        # nécessaire s'il y a deux exécutions Miss dans le même calcul Aster
        self.menage()
        self._dbg_trace("Stop")
        if self.debug:
            print(self.timer)
        return toReturn

    def fichier_resultat(self):
        """Copie les fichiers résultats dans les unités logiques."""
        if self.param["UNITE_IMPR_ASTER"] and osp.exists(self._fichier_tmp("aster")):
            copie_fichier(
                self._fichier_tmp("aster"), self._fichier_aster(self.param["UNITE_IMPR_ASTER"])
            )
        if osp.exists(self._fichier_tmp("resu_impe")):
            copie_fichier(
                self._fichier_tmp("resu_impe"), self._fichier_aster(self.param["UNITE_RESU_IMPE"])
            )
        if osp.exists(self._fichier_tmp("resu_forc")):
            copie_fichier(
                self._fichier_tmp("resu_forc"), self._fichier_aster(self.param["UNITE_RESU_FORC"])
            )

    def cree_reptrav(self):
        """Création du répertoire d'exécution de MISS."""
        if not osp.exists(self.param["_WRKDIR"]):
            os.makedirs(self.param["_WRKDIR"])

    def menage(self):
        """Suppression des fichiers/répertoires de travail."""
        if osp.exists(self.param["_WRKDIR"]) and self.param["REPERTOIRE"] is None:
            shutil.rmtree(self.param["_WRKDIR"])

    def cree_resultat_aster(self):
        """Produit le(s) fichier(s) issu(s) d'Aster."""
        self._dbg_trace("Start")
        ulaster = LogicalUnitFile.new_free()
        mael = self.param["MACR_ELEM_DYNA"]
        if mael is None:
            opts = {}
            if self.param["MATR_RIGI"]:
                opts["MATR_RIGI"] = self.param["MATR_RIGI"]
            if self.param["MATR_MASS"]:
                opts["MATR_MASS"] = self.param["MATR_MASS"]
            __mael = MACR_ELEM_DYNA(BASE_MODALE=self.param["BASE_MODALE"], **opts)
            mael = __mael
        if self.param["_hasPC"]:
            grma = self.param["GROUP_MA_CONTROL"]
            self.param.set("_nbPC", get_number_PC(self.parent, mael, grma))
        other_groups = {}
        if self.param["ISSF"] == "OUI":
            other_groups = _F(
                GROUP_MA_FLU_STR=self.param["GROUP_MA_FLU_STR"],
                GROUP_MA_FLU_SOL=self.param["GROUP_MA_FLU_SOL"],
                GROUP_MA_SOL_SOL=self.param["GROUP_MA_SOL_SOL"],
            )
        else:
            if "GROUP_MA_SOL_SOL" in self.param:
                other_groups = _F(GROUP_MA_SOL_SOL=self.param["GROUP_MA_SOL_SOL"])
        IMPR_MACR_ELEM(
            MACR_ELEM_DYNA=mael,
            FORMAT="MISS_3D",
            GROUP_MA_INTERF=self.param["GROUP_MA_INTERF"],
            GROUP_MA_CONTROL=self.param.get("GROUP_MA_CONTROL"),
            FORMAT_R="1PE16.9",
            SOUS_TITRE="PRODUIT PAR CALC_MISS",
            UNITE=ulaster.unit,
            **other_groups
        )
        ulaster.release()
        copie_fichier(ulaster.filename, self._fichier_tmp("aster"))
        self.data = self.resu_aster_reader.read(self._fichier_tmp("aster"))
        self._dbg_trace("Stop")

    def cree_fichier_mvol(self):
        """Produit le fichier de maillage."""
        self._dbg_trace("Start")
        content = fichier_mvol(self.data)
        with open(self._fichier_tmp("mvol"), MODE) as f:
            f.write(content.encode())
        self._dbg_trace("Stop")

    def cree_fichier_pc(self):
        """Produit les fichiers pour les points de contrôle"""
        if not self.param["_hasPC"]:
            return
        self._dbg_trace("Start")
        content = fichier_ext(self.data)
        with open(self._fichier_tmp("ext"), MODE) as f:
            f.write(content.encode())
        content = fichier_sign(self.param)
        with open(self._fichier_tmp("01.sign"), MODE) as f:
            f.write(content.encode())
        self._dbg_trace("Stop")

    def cree_fichier_chp(self):
        """Produit le fichier chp (modes statiques, dynamiques...)."""
        self._dbg_trace("Start")
        content = fichier_chp(self.param, self.data)
        with open(self._fichier_tmp("chp"), MODE) as f:
            f.write(content.encode())
        self._dbg_trace("Stop")

    def cree_fichier_sol(self):
        """Ecrit le fichier de sol."""
        self._dbg_trace("Stop")
        if self.param["TABLE_SOL"] is not None:
            tabsol = self.param["TABLE_SOL"].EXTR_TABLE()
            sol_content = fichier_sol(tabsol, self.data, self.param)
            if self.verbose:
                _print("Fichier de sol", sol_content)
            with open(self._fichier_tmp("sol"), MODE) as f:
                f.write(sol_content.encode())
        self._dbg_trace("Stop")

    def cree_commande_miss(self, ext="in", lapl_temps=False):
        """Produit le fichier de commandes Miss"""
        self._dbg_trace("Start")
        # Execute méthode classique
        generator = MissCmdeGen(
            self.param, self.data, self._fichier_tmp_local, lapl_temps=lapl_temps
        )
        content = generator.build()
        with open(self._fichier_tmp(ext), MODE) as f:
            f.write(content.encode())
        self._dbg_trace("Stop")

    def cree_fichier_option(self):
        """Ecrit le fichier OPTMIS."""
        self._dbg_trace("Start")
        option_content = fichier_option(self.param)
        if self.verbose:
            _print("Fichier d'option", option_content)
        with open(self._fichier_tmp("optmis"), MODE) as f:
            f.write(option_content.encode())
        self._dbg_trace("Stop")

    # --- utilitaires internes
    def _fichier_tmp(self, ext):
        """Retourne le nom d'un fichier MISS dans WRKDIR"""
        fich = "%s.%s" % (self.param["PROJET"], ext)
        return osp.join(self.param["_WRKDIR"], fich)

    def _fichier_tmp_local(self, ext):
        """Retourne le nom d'un fichier MISS en relatif par rapport au
        répertoire d'exécution de Miss"""
        return osp.join("./", osp.basename(self._fichier_tmp(ext)))

    def _fichier_aster(self, unite):
        """Nom du fichier d'unité logique unite dans le répertoire d'exécution
        de Code_Aster"""
        filename = LogicalUnitFile.filename_from_unit(unite)
        filename = osp.join(self.param["_INIDIR"], filename)
        return filename

    def _dbg_trace(self, on_off):
        if not self.debug:
            return
        stack = traceback.format_stack(limit=5)[-2]
        mat = re.search("File ['\"]*(.*?)['\"]*, *line ([0-9]+), *in (.*)", stack)
        getattr(self.timer, on_off)("CALC_MISS." + mat.group(3))


class CalculMissImpe(CalculMiss):

    """Définition d'un calcul MISS3D de type MISS_IMPE."""

    option_calcul = "MISS_IMPE"


class CalculMissIssf(CalculMiss):

    """Définition d'un calcul MISS3D de type MISS_IMPE avec ISSF."""

    option_calcul = "MISS_IMPE"


class CalculMissPost(CalculMiss):

    """Définition d'une exécution de CALC_MISS où seul le
    post-traitement est demandé."""

    option_calcul = "POST-TRAITEMENT"

    def prepare_donnees(self):
        """Préparation des données."""

    def execute(self):
        """Exécute MISS3D."""

    def fichier_resultat(self):
        """Copie les fichiers résultats dans les unités logiques."""


class CalculMissFichierTemps(CalculMiss):

    """Définition d'une exécution de CALC_MISS dans le domaine de Laplace"""

    option_calcul = "IMPE_LAPL"

    def init_attr(self):
        """Initialisations"""
        rank = MPI.ASTER_COMM_WORLD.Get_rank()
        size = MPI.ASTER_COMM_WORLD.Get_size()
        cwd = osp.join(os.getcwd(), self.param["_WRKDIR"])
        host = socket.gethostname()
        path = "{}:{}".format(host, cwd)
        if size > 1:
            UTMESS("I", "PARALLEL_1", vali=(rank, size), valk=path)
        ipath = encode_str(path)
        # send proc #0 directory to others
        self._results_path = []
        buffer = MPI.ASTER_COMM_WORLD.bcast(ipath, 0)
        self._results_path.append(decode_str(buffer))
        # proc #0 gathers all working directories
        alldirs = MPI.ASTER_COMM_WORLD.gather(path, 0)
        if rank == 0:
            self._results_path.extend(alldirs[1:])

        self.dt = self.param["PAS_INST"]
        N_inst = int(self.param["INST_FIN"] / self.param["PAS_INST"])
        self.reducFactor = self.param["FACTEUR_INTERPOL"]
        self.cutOffValue = self.param["PCENT_FREQ_CALCUL"] / 100.0
        if N_inst % 2 != 0:
            UTMESS("F", "MISS0_18")
        eps = self.param["PRECISION"]
        self.rho = eps ** (1.0 / (2.0 * N_inst))
        factor = self.param["COEF_SURECH"]
        self.L_points = int(factor * N_inst)
        self.nbr_freq = self.L_points // 2 + 1

        # Variables à rajouter dans 'param'
        self.param.set("LIST_FREQ", None)
        self.param.set("FREQ_IMAG", None)
        # Nombre de modes (dynamiques) d'interface
        modes_nb = self.data.mode_stat_nb
        self.param.set("NB_MODE", modes_nb)

        # Creation du fichier sol
        self.param.set("FICHIER_SOL_INCI", self._fichier_tmp_local("sol.inci"))

    def execute(self):
        """Exécute MISS3D : calcul du champ incident + boucle sur les fréquences complexes"""
        self.init_attr()

        if self.param["EXCIT_SOL"]:
            copie_fichier(self._fichier_tmp("inci"), self._fichier_tmp("in"))
            aster.affiche("MESSAGE", "LANCEMENT DU CALCUL DE LA FORCE SISMIQUE EN FREQUENCE")
            CalculMiss.execute(self)
            copie_fichier(self._fichier_tmp("OUT"), self._fichier_tmp("OUT.inci"))
            copie_fichier(self._fichier_tmp("sol"), self._fichier_tmp("sol.inci"))
            copie_fichier(self._fichier_tmp("resu_forc"), self._fichier_tmp("forc_Laplace"))
            copie_fichier(
                self._fichier_tmp("resu_forc"),
                self._fichier_aster(self.param["EXCIT_SOL"]["UNITE_RESU_FORC"]),
            )

        CalculMiss.cree_fichier_sol(self)
        aster.affiche("MESSAGE", "BOUCLE SUR LES FREQUENCES COMPLEXES")
        self._exec_boucle_lapl()
        self.build_global_file()

    def _exec_boucle_lapl(self):
        """Exécute MISS3D dans une boucle sur toutes les fréquences complexes"""
        L_points = self.L_points
        dt = self.dt
        rho = self.rho
        reduc_factor = self.reducFactor
        # TODO la règle UN_PARMI(FREQ_MIN, LIST_FREQ, FREQ_IMAG) ne convient pas!
        # au moins mettre FREQ_MIN à None
        self.param.set("FREQ_MIN", None)
        if self.cutOffValue == 1:
            fc = self.cutOffValue * self.nbr_freq / float(reduc_factor)
        else:
            fc = NP.int_(NP.ceil(self.cutOffValue * self.nbr_freq / float(reduc_factor)))
        freq_list1 = NP.arange(0, fc * reduc_factor)
        freq_list2 = NP.arange(fc * reduc_factor, self.nbr_freq, reduc_factor)
        self.freq_list = freq_list1.tolist() + freq_list2.tolist()
        rank = MPI.ASTER_COMM_WORLD.Get_rank()
        size = MPI.ASTER_COMM_WORLD.Get_size()

        for k in self.freq_list:
            # round-robin partition
            if k % size != rank:
                UTMESS("I", "MISS0_24", vali=(k, k % size))
                continue
            UTMESS("I", "MISS0_25", vali=(k, rank))

            if (k == 0) or (k == self.nbr_freq - 1):
                self.param.set("LIST_FREQ", (0.1e-4,))
                self.param.set(
                    "FREQ_IMAG",
                    (
                        1.5
                        - 2.0 * rho * cos(2 * pi * k / L_points)
                        + 0.5 * rho * rho * cos(4 * pi * k / L_points)
                    )
                    / dt,
                )
            else:
                self.param.set(
                    "LIST_FREQ",
                    (
                        -(
                            -2.0 * rho * sin(2 * pi * k / L_points)
                            + 0.5 * rho * rho * sin(4 * pi * k / L_points)
                        )
                        / dt
                        / (2 * pi),
                    ),
                )
                self.param.set(
                    "FREQ_IMAG",
                    (
                        1.5
                        - 2.0 * rho * cos(2 * pi * k / L_points)
                        + 0.5 * rho * rho * cos(4 * pi * k / L_points)
                    )
                    / dt,
                )

            self.param.set("_calc_impe", True)
            CalculMiss.cree_commande_miss(self)
            str00 = (
                str(self.param["FREQ_IMAG"])
                + " + i . "
                + str(self.param["LIST_FREQ"][0])
                + " ("
                + str(k + 1)
                + "/"
                + str(self.nbr_freq)
                + ")"
            )
            aster.affiche("MESSAGE", "FREQUENCE COMPLEXE COURANTE =  " + str00)
            CalculMiss.execute(self)
            resname = self._fichier_tmp("resu_impe_%04d" % k)
            copie_fichier(self._fichier_tmp("resu_impe"), resname)

            if rank != 0:
                UTMESS("I", "PARALLEL_2", valk=(resname, self._results_path[0]))
                send_file(resname, self._results_path[0])

        # libérer la structure contenant les données numériques
        CalculMiss.init_data(self)
        MPI.ASTER_COMM_WORLD.Barrier()

    def build_global_file(self):
        """Build the file by concatenating those on each frequency"""
        rank = MPI.ASTER_COMM_WORLD.Get_rank()
        size = MPI.ASTER_COMM_WORLD.Get_size()
        if rank == 0:
            fimpe = self._fichier_tmp("impe_Laplace")
            fd = open(fimpe, "w")
            for k in self.freq_list:
                resname = self._fichier_tmp("resu_impe_%04d" % k)
                with open(resname, "r") as f:
                    fd.write(f.read())
            fd.close()
            for k in range(1, size):
                UTMESS("I", "PARALLEL_2", valk=(fimpe, self._results_path[k]))
                send_file(fimpe, self._results_path[k])
        MPI.ASTER_COMM_WORLD.Barrier()

    def cree_commande_miss(self):
        """Produit le fichier de commandes Miss du champ incident (.inci)"""
        CalculMiss.cree_commande_miss(self, ext="inci", lapl_temps=True)

    def fichier_resultat(self):
        """Libérer la structure contenant les données numériques"""
        if self.param["UNITE_IMPR_ASTER"] and osp.exists(self._fichier_tmp("aster")):
            copie_fichier(
                self._fichier_tmp("aster"), self._fichier_aster(self.param["UNITE_IMPR_ASTER"])
            )

    def init_data(self):
        """Libérer la structure contenant les données numériques"""
        pass


def get_number_PC(parent, macr_elem, lgrpc):
    """Retourne le nombre de points de contrôle"""
    mail = macr_elem.getDOFNumbering().getMesh()
    assert mail is not None, "impossible de récupérer le maillage du macro-élément"
    lgrpma = mail.LIST_GROUP_MA()
    result = sum([nbel for name, nbel in lgrpma if name in lgrpc])
    return result
