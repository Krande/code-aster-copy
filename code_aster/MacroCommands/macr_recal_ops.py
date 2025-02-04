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

# person_in_charge: mathieu.courtois@edf.fr

import copy
import glob
import os
import sys

import numpy as NP
from asrun.common.sysutils import on_64bits
from asrun.profil import AsterProfil
from run_aster.config import CFG

from libaster import onFatalError

from ..CodeCommands import DEFI_LIST_REEL
from ..Messages import UTMESS
from .Recal import reca_algo, reca_calcul_aster, reca_interp, reca_message, reca_utilitaires, recal
from .Recal.reca_controles import gestion
from .Recal.reca_evol import evolutivo
from .Utils.optimize import fmin, fminBFGS, fminNCG

try:
    import Gnuplot

    HAS_GNUPLOT = True
except ImportError:
    HAS_GNUPLOT = False

debug = False

INFO = 1
NOMPRO = "MACR_RECAL"


# ------------------------------------------------------------------------
def Sortie(LIST_NOM_PARA, LIST_PARA, val, CALCUL_ASTER, Mess):
    """Sortie de la macro, on renvoie les parametres obtenus"""

    UTMESS("I", "RECAL0_39", valk=str(CALCUL_ASTER.evaluation_fonction), files=Mess.get_filename())

    LIST_NOM_PARA_ALPHA = [para[0] for para in LIST_PARA]
    LIST_NOM_PARA_ALPHA.sort()
    lival = []
    for i in LIST_NOM_PARA:
        lival.append(val[LIST_NOM_PARA_ALPHA.index(i)])
    nomres = DEFI_LIST_REEL(VALE=lival)

    return nomres


# ------------------------------------------------------------------------
def force_list(obj, typref=list):
    """Retourne 'obj' sous forme d'une liste de 'typref'."""
    if type(obj) not in (list, tuple):
        assert type(obj) == typref, "%s != %s" % (type(obj), typref)
        obj = [obj]
    elif len(obj) > 0:
        elt = obj[0]
        if type(elt) != typref:
            obj = [obj]
    return obj


# ------------------------------------------------------------------------
def macr_recal_ops(
    self,
    UNITE_ESCL,
    RESU_EXP=None,
    LIST_POIDS=None,
    LIST_PARA=None,
    RESU_CALC=None,
    ITER_MAXI=None,
    ITER_FONC_MAXI=None,
    RESI_GLOB_RELA=None,
    UNITE_RESU=None,
    PARA_DIFF_FINI=None,
    GRAPHIQUE=None,
    PARA_OPTI=None,
    COURBE=None,
    METHODE=None,
    INFO=None,
    **args
):
    """Macro commande realisant le recalage de modeles Aster"""
    # Gestion des Exceptions
    prev_onFatalError = onFatalError()
    onFatalError("EXCEPTION")

    # Reecriture des mots-cles
    if COURBE:
        RESU_EXP = []
        RESU_CALC = []
        LIST_POIDS = []
        LIST_POIDS_DEFAUT = []
        LIST_POIDS_TMP = []
        for m in COURBE:
            lpar, lval = m["FONC_EXP"].Valeurs()
            LIST_EXP = []
            LIST_CALC_TMP = []
            for j in range(len(lpar)):
                LIST_EXP_TMP = []
                LIST_EXP_TMP.append(lpar[j])
                LIST_EXP_TMP.append(lval[j])
                LIST_EXP.append(LIST_EXP_TMP)
            RESU_EXP.append(NP.array(LIST_EXP))
            LIST_CALC_TMP.append(m["NOM_FONC_CALC"])
            LIST_CALC_TMP.append(m["PARA_X"])
            LIST_CALC_TMP.append(m["PARA_Y"])
            RESU_CALC.append(LIST_CALC_TMP)
            if m["POIDS"] is not None:
                LIST_POIDS_TMP.append(m["POIDS"])
            LIST_POIDS_DEFAUT.append(1.0)
        if len(LIST_POIDS_TMP) == 0:
            LIST_POIDS = NP.array(LIST_POIDS_DEFAUT)
        else:
            LIST_POIDS = NP.array(LIST_POIDS_TMP)

    if PARA_OPTI:
        j = 0
        LIST_PARA = []
        for m in PARA_OPTI:
            LIST_TMP = []
            LIST_TMP.append(m["NOM_PARA"])
            LIST_TMP.append(m["VALE_INI"])
            LIST_TMP.append(m["VALE_MIN"])
            LIST_TMP.append(m["VALE_MAX"])
            LIST_PARA.append(LIST_TMP)

    nomres = macr_recal(
        self,
        UNITE_ESCL,
        force_list(RESU_EXP, NP.ndarray),
        LIST_POIDS,
        force_list(LIST_PARA),
        force_list(RESU_CALC),
        ITER_MAXI,
        ITER_FONC_MAXI,
        RESI_GLOB_RELA,
        UNITE_RESU,
        PARA_DIFF_FINI,
        GRAPHIQUE,
        METHODE,
        INFO,
        **args
    )

    onFatalError(prev_onFatalError)
    return nomres


# ------------------------------------------------------------------------
def macr_recal(
    self,
    UNITE_ESCL,
    RESU_EXP,
    POIDS,
    LIST_PARA,
    RESU_CALC,
    ITER_MAXI,
    ITER_FONC_MAXI,
    RESI_GLOB_RELA,
    UNITE_RESU,
    PARA_DIFF_FINI,
    GRAPHIQUE,
    METHODE,
    INFO,
    **args
):
    ASTER_ROOT = os.environ["ASTER_ROOT"]

    try:
        sys.path.append(os.path.join(ASTER_ROOT, "ASTK", "ASTK_SERV", "lib"))
        sys.path.append(
            os.path.join(
                ASTER_ROOT,
                "lib",
                "python%s.%s" % (sys.version_info[0], sys.version_info[1]),
                "site-packages",
            )
        )
    except:
        pass

    # _____________________________________________
    #
    # RECUPERATION DU PROFIL DU CALCUL MAITRE
    # _____________________________________________
    # Lecture du fichier .export dans le repertoire temporaire d'execution
    list_export = glob.glob("*.export")
    if len(list_export) == 0:
        UTMESS("F", "RECAL0_4")
    elif len(list_export) > 1:
        UTMESS("F", "RECAL0_5")
    prof = AsterProfil(list_export[0])

    # _____________________________________________
    #
    # PARAMETRES
    # _____________________________________________
    TOLE_PARA = args["TOLE_PARA"]
    TOLE_FONC = args["TOLE_FONC"]

    # Pour les calculs esclaves
    CALCUL_ESCLAVE = {}.fromkeys(
        [
            "LANCEMENT",
            "MODE",
            "UNITE_SUIVI",
            "CLASSE",
            "ACTUALISATION",
            "memjeveux",
            "memjob",
            "mem_aster",
            "tpmax",
            "tpsjob",
            "mpi_nbnoeud",
            "mpi_nbcpu",
            "NMAX_SIMULT",
        ]
    )

    dESCLAVE = args["CALCUL_ESCLAVE"][0].cree_dict_valeurs(args["CALCUL_ESCLAVE"][0].mc_liste)
    for i in list(dESCLAVE.keys()):
        if dESCLAVE[i] is None:
            del dESCLAVE[i]

    CALCUL_ESCLAVE["LANCEMENT"] = dESCLAVE["LANCEMENT"]
    if "UNITE_SUIVI" in dESCLAVE:
        CALCUL_ESCLAVE["UNITE_SUIVI"] = dESCLAVE["UNITE_SUIVI"]
    else:
        CALCUL_ESCLAVE["UNITE_SUIVI"] = None
    if "MODE" in dESCLAVE:
        CALCUL_ESCLAVE["MODE"] = dESCLAVE["MODE"]
    else:
        CALCUL_ESCLAVE["MODE"] = prof["mode"][0].upper()

    LANCEMENT = CALCUL_ESCLAVE["LANCEMENT"]

    # Parametres de l'algorithme genetique
    if "NB_PARENTS" in args:
        NB_PARENTS = args["NB_PARENTS"]
    if "NB_FILS" in args:
        NB_FILS = args["NB_FILS"]
    if "ECART_TYPE" in args:
        ECART_TYPE = args["ECART_TYPE"]
    if "ITER_ALGO_GENE" in args:
        ITER_ALGO_GENE = args["ITER_ALGO_GENE"]
    if "RESI_ALGO_GENE" in args:
        RESI_ALGO_GENE = args["RESI_ALGO_GENE"]

    if "GRAINE" in args:
        UTMESS("A", "RECAL0_43")
        GRAINE = args["GRAINE"]
    else:
        GRAINE = None

    # Parametres concernant le recalage d'un modele dynamique
    if "DYNAMIQUE" in args:
        DYNAMIQUE = args["DYNAMIQUE"]
    else:
        DYNAMIQUE = None

    # _____________________________________________
    #
    # VERIFICATION PREALABLE SUR GNUPLOT
    # _____________________________________________

    if GRAPHIQUE:
        dGRAPHIQUE = GRAPHIQUE[0].cree_dict_valeurs(GRAPHIQUE[0].mc_liste)
        if "FORMAT" in dGRAPHIQUE and dGRAPHIQUE["FORMAT"] == "GNUPLOT":
            # On essaie d'importer Gnuplot -> PAS DE GRAPHIQUE
            if not HAS_GNUPLOT:
                GRAPHIQUE = None
                UTMESS("A", "RECAL0_3")

    # _____________________________________________
    #
    # PARAMETRES DU MODE DISTRIBUTION
    # _____________________________________________
    if LANCEMENT == "DISTRIBUTION":
        addmem = CFG.get("addmem", 0)
        if debug:
            print("tpsjob:", prof.param["tpsjob"][0])
            print("tpmax:", prof.args["tpmax"])
            print("mem_aster:", prof.param.get("mem_aster", [100])[0])
            print("memjeveux:", prof.args["memjeveux"])
            print("memjob:", prof.param["memjob"][0])
            print("addmem:", addmem)
            print("dESCLAVE:", dESCLAVE)

        # Pour la conversion mega-mots / mega-octets
        if on_64bits():
            facw = 8
        else:
            facw = 4

        # Recuperation du parametre mem_aster
        try:
            mem_aster = int(prof["mem_aster"][0])
        except ValueError:
            mem_aster = 100
        if mem_aster in (0, 100):
            if CALCUL_ESCLAVE["MODE"] == "INTERACTIF":
                UTMESS("A", "RECAL0_6")
            mem_aster = 100
        CALCUL_ESCLAVE["mem_aster"] = mem_aster

        # Utilisation du mot-cle TEMPS
        if "TEMPS" in dESCLAVE:
            CALCUL_ESCLAVE["tpsjob"] = int(dESCLAVE["TEMPS"] / 60)
            CALCUL_ESCLAVE["tpmax"] = int(dESCLAVE["TEMPS"])
        else:
            # Recuperation depuis le calcul maitre
            CALCUL_ESCLAVE["tpsjob"] = prof.param["tpsjob"][0]
            CALCUL_ESCLAVE["tpmax"] = prof.args["tpmax"]

        # Utilisation du mot-cle MEMOIRE
        if "MEMOIRE" in dESCLAVE:
            mem = int(dESCLAVE["MEMOIRE"]) + addmem
            CALCUL_ESCLAVE["memjob"] = int(mem * 1024)
            # Calcul du parametre memjeveux esclave
            memjeveux = int(mem / facw)
            try:
                if mem_aster == 100:
                    CALCUL_ESCLAVE["memjeveux"] = memjeveux
                else:
                    CALCUL_ESCLAVE["memjeveux"] = float(
                        int((float(mem_aster) / 100.0) * float(memjeveux))
                    )
            except:
                UTMESS("F", "RECAL0_8")
        else:
            # Recuperation depuis le calcul maitre
            CALCUL_ESCLAVE["memjob"] = int(prof.param["memjob"][0]) + addmem * 1024
            CALCUL_ESCLAVE["memjeveux"] = prof.args["memjeveux"] + int(addmem / facw)

        # Utilisation du mot-cle MPI_NBCPU
        if "MPI_NBCPU" in dESCLAVE:
            # Verifie que le calcul maitre est bien en MPI sur 1 cpu
            mpi_nbcpu = str(prof["mpi_nbcpu"][0])
            if mpi_nbcpu != "1":
                UTMESS("A", "RECAL0_7")

            CALCUL_ESCLAVE["mpi_nbcpu"] = int(dESCLAVE["MPI_NBCPU"])

        # Utilisation du mot-cle MPI_NBNOEUD
        if "MPI_NBNOEUD" in dESCLAVE:
            CALCUL_ESCLAVE["mpi_nbnoeud"] = int(dESCLAVE["MPI_NBNOEUD"])

        # Parametres batch
        if CALCUL_ESCLAVE["MODE"] == "BATCH":
            if "CLASSE" in dESCLAVE:
                CALCUL_ESCLAVE["CLASSE"] = dESCLAVE["CLASSE"]
            if "ACTUALISATION" in dESCLAVE:
                CALCUL_ESCLAVE["ACTUALISATION"] = dESCLAVE["ACTUALISATION"]

            # Affichage parametres batch
            if CALCUL_ESCLAVE["CLASSE"]:
                classe = CALCUL_ESCLAVE["CLASSE"]
            else:
                classe = " -auto- "
            valk = (
                str(CALCUL_ESCLAVE["tpmax"]),
                str(int(CALCUL_ESCLAVE["memjob"]) / 1024),
                str(int(float(CALCUL_ESCLAVE["memjeveux"]) * facw)),
                classe,
            )
            UTMESS("I", "RECAL0_69", valk=valk)
        if debug:
            print("CALCUL_ESCLAVE:", CALCUL_ESCLAVE)

    # _____________________________________________
    #
    # VERIFICATIONS
    # _____________________________________________

    if float(PARA_DIFF_FINI) > 0.1:
        UTMESS("A", "RECAL0_76", valk=(str(PARA_DIFF_FINI)))

    # _____________________________________________
    #
    # INITIALISATIONS
    # _____________________________________________
    # Stocke l'ordre initial des parametres pour restituer dans le bon ordre
    # les valeurs en sortie de la macro
    LIST_NOM_PARA = [para[0] for para in LIST_PARA]

    # On classe les parametres
    LIST_PARA.sort()

    # Pour les algorithmes d'optimize.py, on a des limitations
    if METHODE in ["FMIN", "FMINBFGS", "FMINNCG"]:
        # On ne peut tracer qu'a la derniere iteration
        if GRAPHIQUE:
            if GRAPHIQUE["AFFICHAGE"] == "TOUTE_ITERATION":
                UTMESS("I", "RECAL0_10", valk=METHODE)
        # Les bornes ne sont pas gerees
        UTMESS("I", "RECAL0_11", valk=METHODE)

    # _______________________________________________
    #
    # GESTION DE L'OPTION FACULTATIVE POUR LES POIDS
    # _______________________________________________
    if POIDS is None:
        POIDS = NP.ones(len(RESU_EXP))

    # _____________________________________________
    #
    # GESTION DES ERREURS DE SYNTAXE
    # _____________________________________________
    texte_erreur, texte_alarme = gestion(
        UNITE_ESCL, LIST_PARA, RESU_CALC, RESU_EXP, POIDS, GRAPHIQUE, UNITE_RESU, METHODE
    )
    if texte_erreur != "":
        UTMESS("F", "RECAL0_12", valk=texte_erreur)
    if texte_alarme != "":
        UTMESS("A", "RECAL0_12", valk=texte_alarme)

    # _____________________________________________
    #
    # INITIALISATIONS
    # _____________________________________________
    iter = 0
    restant, temps_iter = 0.0, 0.0
    restant, temps_iter, err = reca_utilitaires.temps_CPU(restant, temps_iter)
    para, val, borne_inf, borne_sup = reca_utilitaires.transforme_list_Num(LIST_PARA, RESU_EXP)
    val_init = copy.copy(val)

    # Fonctionnelle en sortie (vectorielle ou scalaire)
    if METHODE in ["FMIN", "FMINBFGS", "FMINNCG", "GENETIQUE", "HYBRIDE"]:
        vector_output = False
    else:
        vector_output = True

    # OBJET "CALCUL"
    CALCUL_ASTER = reca_calcul_aster.CALCUL_ASTER(
        jdc=self,
        METHODE=METHODE,
        UNITE_ESCL=UNITE_ESCL,
        UNITE_RESU=UNITE_RESU,
        para=para,
        reponses=RESU_CALC,
        PARA_DIFF_FINI=PARA_DIFF_FINI,
        vector_output=vector_output,
        DYNAMIQUE=DYNAMIQUE,
        # LANCEMENT       = LANCEMENT,
        CALCUL_ESCLAVE=CALCUL_ESCLAVE,
        INFO=INFO,
    )

    CALCUL_ASTER.RESU_EXP = RESU_EXP
    CALCUL_ASTER.RESU_CALC = RESU_CALC
    CALCUL_ASTER.LIST_PARA = LIST_PARA

    if CALCUL_ESCLAVE["UNITE_SUIVI"]:
        CALCUL_ASTER.unity_follow = CALCUL_ESCLAVE["UNITE_SUIVI"]

    # Instances des classes pour le calcul de l'erreur et le
    # dimensionnemnt/adim
    Dim = reca_algo.Dimension(copy.copy(val_init))
    CALCUL_ASTER.Simul = reca_interp.Sim_exp(RESU_EXP, POIDS)
    CALCUL_ASTER.Dim = Dim
    CALCUL_ASTER.reca_algo = reca_algo

    if GRAPHIQUE:
        CALCUL_ASTER.UNITE_GRAPHIQUE = GRAPHIQUE["UNITE"]

    # Dans le cas de la dynamique avec appariement manual des MAC, on passe la
    # flag correspondant a True
    if METHODE in [
        "HYBRIDE",
        "LEVENBERG",
        "GENETIQUE",
    ]:  # AAC --> j'ai modifie et donne la possibilite d'afficher la fenetre mac pour levenb et gene
        if DYNAMIQUE is not None and DYNAMIQUE["APPARIEMENT_MANUEL"] == "OUI":
            CALCUL_ASTER.graph_mac = True

    # Instance de la classe gérant l'affichage des resultats du calcul de
    # l'optimisation
    Mess = reca_message.Message(para, RESU_EXP, copy.copy(val_init), UNITE_RESU)
    Mess.initialise()

    # Calcul de F
    #    erreur = CALCUL_ASTER.calcul_F(val)
    # Calcul de F et G
    #    erreur, residu, A_nodim, A = CALCUL_ASTER.calcul_FG(val)
    #    sys.exit()

    # Mode INCLUDE : on doit executer les commandes PRE ici
    if LANCEMENT == "INCLUSION":
        UNITE_INCLUDE = UNITE_ESCL
        recal.make_include_files(
            UNITE_INCLUDE=UNITE_INCLUDE, calcul=RESU_CALC, parametres=LIST_PARA
        )

    # -------------------------------------------------------------------------------
    # Pas d'optimisation (juste une evaluation de la fonctionnelle pour le point courant)
    # -------------------------------------------------------------------------------
    #
    if ITER_MAXI <= 0:
        erreur = CALCUL_ASTER.calcul_F(val)
        residu = 0
        iter = 0
        L_F = CALCUL_ASTER.Lcalc[0]
        CALCUL_ASTER.evaluation_fonction = 1

    # -------------------------------------------------------------------------------
    # Algorithme FMIN (pas d'adimensionnement car n'utilise pas de gradient)
    # -------------------------------------------------------------------------------
    #
    elif METHODE == "FMIN":
        UTMESS("I", "RECAL0_13", valk=METHODE, files=Mess.get_filename())
        val, fval, warnflag = fmin(
            CALCUL_ASTER.calcul_F, val, maxiter=ITER_MAXI, maxfun=ITER_FONC_MAXI, fulloutput=1
        )

        iter_fonc = CALCUL_ASTER.evaluation_fonction
        if warnflag == 1:
            UTMESS("I", "RECAL0_54", files=Mess.get_filename())
        if warnflag == 2:
            UTMESS("I", "RECAL0_55", files=Mess.get_filename())
        Mess.affiche_etat_final_convergence(
            iter, ITER_MAXI, iter_fonc, ITER_FONC_MAXI, RESI_GLOB_RELA, residu=0, Act=[]
        )
        Mess.affiche_fonctionnelle(fval)
        Mess.affiche_valeurs(val)
        nomres = Sortie(LIST_NOM_PARA, LIST_PARA, val, CALCUL_ASTER, Mess)
        return nomres

    # -------------------------------------------------------------------------------
    # Algorithme GENETIQUE (pas d'adimensionnement car n'utilise pas de gradient)
    # -------------------------------------------------------------------------------
    #
    elif METHODE == "GENETIQUE":
        UTMESS("I", "RECAL0_13", valk=METHODE, files=Mess.get_filename())
        nb_parents = NB_PARENTS
        nb_fils = NB_FILS
        nb_iter = ITER_ALGO_GENE
        sigma = ECART_TYPE
        err_min = RESI_ALGO_GENE
        graine = GRAINE
        val = evolutivo(
            CALCUL_ASTER,
            val,
            nb_iter,
            err_min,
            nb_parents,
            nb_fils,
            sigma,
            borne_inf,
            borne_sup,
            graine,
        )
        nomres = Sortie(LIST_NOM_PARA, LIST_PARA, val, CALCUL_ASTER, Mess)
        return nomres

    # -------------------------------------------------------------------------------
    # Pour tous les autres methodes, on adimensionne
    # -------------------------------------------------------------------------------
    #
    else:
        # -------------------------------------------------------------------------------
        # Si METHODE=='HYBRIDE', on lance d'abord l'algo genetique et ensuite celui de
        # Levenberg-Marquardt qui demarre avec le jeu de parametres issu de
        # genetique
        if METHODE == "HYBRIDE":
            nb_parents = NB_PARENTS
            nb_fils = NB_FILS
            nb_iter = ITER_ALGO_GENE
            sigma = ECART_TYPE
            err_min = RESI_ALGO_GENE
            graine = GRAINE
            val_gene = evolutivo(
                CALCUL_ASTER,
                val,
                nb_iter,
                err_min,
                nb_parents,
                nb_fils,
                sigma,
                borne_inf,
                borne_sup,
                graine,
            )
            val = copy.copy(val_gene)
            val_init = copy.copy(val)
            # AA ? CALCUL_ASTER.graph_mac = True

        # Calcul de F et G
        erreur, residu, A_nodim, A = CALCUL_ASTER.calcul_FG(val)
        E = recal.CALC_ERROR(experience=RESU_EXP, X0=val, calcul=RESU_CALC, poids=POIDS)
        E.CalcError(CALCUL_ASTER.Lcalc)
        E.CalcSensibilityMatrix(CALCUL_ASTER.Lcalc, val, dX=None, pas=PARA_DIFF_FINI)

        L_init = E.L_init
        L_J_init = E.L_J_init
        J_init = E.J_init
        J = E.J
        A = E.A
        A_nodim = E.A_nodim
        erreur = E.erreur
        residu = E.residu
        gradient_init = E.gradient_init

        # Calcul du lambda_init
        l = reca_algo.lambda_init(NP.dot(NP.transpose(A), A))

        Mess.affiche_result_iter(iter, J, val, residu, NP.array([]))

        CALCUL_ASTER.L_init = L_init
        CALCUL_ASTER.L_J_init = L_J_init
        CALCUL_ASTER.J_init = J_init
        CALCUL_ASTER.A_init = A
        CALCUL_ASTER.gradient_init = gradient_init
        CALCUL_ASTER.residu_init = residu

        # On teste un manque de temps CPU
        restant, temps_iter, err = reca_utilitaires.temps_CPU(restant, temps_iter)
        if err == 1:
            ier = ier + 1
            return ier

        # -------------------------------------------------------------------------------
        # Methode FMINBFGS et FMINNCG
        # -------------------------------------------------------------------------------
        #
        if METHODE in ["FMINBFGS", "FMINNCG"]:
            UTMESS("I", "RECAL0_13", valk=METHODE, files=Mess.get_filename())

            # Derivees
            f = CALCUL_ASTER.calcul_F2
            fprime = CALCUL_ASTER.calcul_G
            warnflag = 0

            if "GRADIENT" in args and args["GRADIENT"] == "NON_CALCULE":
                f = CALCUL_ASTER.calcul_F
                fprime = None

            if fprime:
                UTMESS("I", "RECAL0_14")
            else:
                UTMESS("I", "RECAL0_15")

            # Lancement de l'optimisation
            if METHODE == "FMINBFGS":
                val, fval, func_calls, grad_calls, warnflag = fminBFGS(
                    f=f,
                    x0=val,
                    fprime=fprime,
                    maxiter=ITER_MAXI,
                    avegtol=RESI_GLOB_RELA,
                    fulloutput=1,
                )

            elif METHODE == "FMINNCG":
                val, fval, func_calls, grad_calls, hcalls, warnflag = fminNCG(
                    f=f,
                    x0=val,
                    fprime=fprime,
                    fhess_p=None,
                    fhess=None,
                    maxiter=ITER_MAXI,
                    avextol=RESI_GLOB_RELA,
                    fulloutput=1,
                )

            # Affichage des messages de sortie
            iter_fonc = CALCUL_ASTER.evaluation_fonction
            if warnflag:
                UTMESS("I", "RECAL0_55", files=Mess.get_filename())
            Mess.affiche_etat_final_convergence(
                iter, ITER_MAXI, iter_fonc, ITER_FONC_MAXI, RESI_GLOB_RELA, residu=0, Act=[]
            )
            Mess.affiche_fonctionnelle(fval)
            Mess.affiche_valeurs(val)

            # Permet d'avoir un diagnostic NOOK pour le job
            if warnflag:
                iter = ITER_MAXI

            L_F = CALCUL_ASTER.L
            residu = fval
            ecart_fonc = 0  # non calcule avec ces methodes
            ecart_para = 0  # non calcule avec ces methodes

        # -------------------------------------------------------------------------------
        # Methode Levenberg-Marquardt
        # ----------------------------------------------------------------------
        elif METHODE in ["LEVENBERG", "HYBRIDE"]:
            # ___________________________________________________________
            #
            # BOUCLE PRINCIPALE DE L'ALGORITHME de Levenberg-Marquardt
            # ___________________________________________________________

            UTMESS("I", "RECAL0_13", valk=METHODE, files=Mess.get_filename())
            while iter < ITER_MAXI:
                iter = iter + 1
                new_val, s, l, Act = reca_algo.Levenberg_bornes(
                    val, Dim, val_init, borne_inf, borne_sup, A, erreur, l, UNITE_RESU
                )

                # On teste la variation sur les parametres
                ecart_para = reca_algo.calcul_norme2(NP.array(new_val) - NP.array(val))
                if debug:
                    print(
                        "AA0/ecart para=%s\nAA0/oldpara/newpara=%s %s" % (ecart_para, val, new_val)
                    )
                if ecart_para < TOLE_PARA:
                    UTMESS("I", "RECAL0_51", valr=ecart_para, files=Mess.get_filename())
                    break

                # Calculs au point courant val et toutes les perturbations par
                # differences finies (N+1 calculs distribues ou inclus)
                CALCUL_ASTER.calcul_FG(new_val)

                # Calcul de l'erreur et de la matrice des sensibilites
                old_J = copy.copy(J)
                E.CalcError(CALCUL_ASTER.Lcalc)
                new_J = E.J

                l = reca_algo.actualise_lambda(
                    l, Dim.adim(val), Dim.adim(new_val), A, erreur, new_J, J
                )
                E.CalcSensibilityMatrix(CALCUL_ASTER.Lcalc, new_val, dX=None, pas=PARA_DIFF_FINI)

                L_F = CALCUL_ASTER.Lcalc[0]
                A = E.A_nodim
                val = copy.copy(new_val)
                erreur = copy.copy(E.erreur)
                J = E.J

                if debug:
                    print("AA0/L_F=", L_F)
                    print("AA0/l=", l)
                    print("AA0/erreur=", erreur)
                    print("AA0/J=", J)
                    print("AA0/A_nodim=", A)

                # Calcul de la matrice des sensibilites
                A = Dim.adim_sensi(A)

                # Calcul du residu
                residu = reca_algo.test_convergence(gradient_init, erreur, A, s)

                if debug:
                    print("AA0/residu=", residu)
                    print("AA0/new_val=", new_val)
                    print("AA0/A=", A)

                # On calcule la variation sur la fonctionnelle
                ecart_fonc = abs(new_J - old_J)

                # Affichage iteration
                Mess.affiche_result_iter(iter, J, val, residu, Act, ecart_para, ecart_fonc)

                # On teste la variation sur la fonctionnelle
                if ecart_fonc < TOLE_FONC:
                    UTMESS("I", "RECAL0_52", valr=ecart_fonc, files=Mess.get_filename())
                    break

                if GRAPHIQUE:
                    if GRAPHIQUE["AFFICHAGE"] == "TOUTE_ITERATION":
                        GRAPHE_UL_OUT = GRAPHIQUE["UNITE"]
                        if "FORMAT" in dGRAPHIQUE and dGRAPHIQUE["FORMAT"] == "XMGRACE":
                            pilote = GRAPHIQUE["PILOTE"]
                        else:
                            pilote = "INTERACTIF"
                        reca_utilitaires.graphique(
                            GRAPHIQUE["FORMAT"],
                            L_F,
                            RESU_EXP,
                            RESU_CALC,
                            iter,
                            GRAPHE_UL_OUT,
                            pilote,
                        )

                # On teste le residu
                if residu <= RESI_GLOB_RELA:
                    UTMESS("I", "RECAL0_50", valr=residu, files=Mess.get_filename())
                    break

                # On teste un manque de temps CPU
                restant, temps_iter, err = reca_utilitaires.temps_CPU(restant, temps_iter)
                if err == 1:
                    UTMESS("I", "RECAL0_53", files=Mess.get_filename())
                    break

            # _____________________________________________
            #
            # FIN DES ITERATIONS
            # CONVERGENCE OU ECHEC
            # _____________________________________________
            iter_fonc = CALCUL_ASTER.evaluation_fonction
            Mess.affiche_etat_final_convergence(
                iter, ITER_MAXI, iter_fonc, ITER_FONC_MAXI, RESI_GLOB_RELA, residu, Act
            )
            reca_algo.calcul_etat_final(para, A, iter, ITER_MAXI, RESI_GLOB_RELA, residu, Mess)

        # ----------------------------------------------------------------------
    # _____________________________________________
    #
    # FIN DES ITERATIONS POUR TOUS LES ALGOS
    # _____________________________________________
    if GRAPHIQUE:
        fichier = None
        # Pour les algorithmes d'optimize.py, on ne peut tracer qu'a la
        # derniere iteration
        if (
            (GRAPHIQUE["AFFICHAGE"] == "ITERATION_FINALE")
            or (METHODE in ["FMIN", "FMINBFGS", "FMINNCG"])
            or (ITER_MAXI <= 0)
        ):
            UTMESS("I", "RECAL0_17")
            GRAPHE_UL_OUT = GRAPHIQUE["UNITE"]
            pilote = GRAPHIQUE["PILOTE"]
            reca_utilitaires.graphique(
                GRAPHIQUE["FORMAT"], L_F, RESU_EXP, RESU_CALC, iter, GRAPHE_UL_OUT, pilote, fichier
            )

    # Si pas de convergence alors diagnostic NOOK_TEST_RESU
    # if (residu > RESI_GLOB_RELA) and  (ecart_fonc > TOLE_FONC) and
    # (ecart_para < TOLE_PARA):

    if debug:
        print("residu, RESI_GLOB_RELA=", residu, RESI_GLOB_RELA, (residu > RESI_GLOB_RELA))
        print("ecart_fonc, TOLE_FONC=", ecart_fonc, TOLE_FONC, (ecart_fonc > TOLE_FONC))
        print("ecart_para, TOLE_PARA=", ecart_para, TOLE_PARA, (ecart_para > TOLE_PARA))

    if residu > RESI_GLOB_RELA:
        _tmp = []
        _tmp.append({"PARA": "ITER_MAXI", "LISTE_R": 0.0})

    # _____________________________________________
    #
    # CREATIONS DE LA LISTE DE REELS CONTENANT
    # LES VALEURS DES PARAMETRES A CONVERGENCE
    # _____________________________________________
    nomres = Sortie(LIST_NOM_PARA, LIST_PARA, val, CALCUL_ASTER, Mess)
    return nomres
