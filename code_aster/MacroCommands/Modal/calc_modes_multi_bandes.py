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

import aster
from ...Messages import UTMESS

from ...Cata.Syntax import _F
from ...CodeCommands import EXTR_MODE, IMPR_CO, INFO_MODE, MODI_MODELE, NUME_DDL
from ...Objects import (
    AssemblyMatrixDisplacementReal,
    GeneralizedAssemblyMatrixReal,
    GeneralizedAssemblyMatrixComplex,
)
from ...Utilities.mpi_utils import MPI
from .mode_iter_simult import MODE_ITER_SIMULT


def calc_modes_multi_bandes(self, stop_erreur, sturm, INFO, **args):
    """
    Macro-command CALC_MODES, case of the simultaneous iterations method
    over several frequency bands, with optional parallelization.
    Can be used only in the case of vibration modes (TYPE_RESU='DYNAMIQUE')
    """
    args = _F(args)
    SOLVEUR = args.get("SOLVEUR")
    SOLVEUR_MODAL = args.get("SOLVEUR_MODAL")
    VERI_MODE = args.get("VERI_MODE")

    MATR_RIGI = args["MATR_RIGI"]
    MATR_MASS = args["MATR_MASS"]
    CALC_FREQ = args["CALC_FREQ"]
    METHODE = SOLVEUR_MODAL["METHODE"]
    STOP_BANDE_VIDE = "NON"
    CARA_ELEM = args.get("CARA_ELEM")
    CHAM_MATER = args.get("CHAM_MATER")

    # ----------------------------------------------------------------------
    #
    # 1a. INITIALISATIONS
    #
    # ----------------------------------------------------------------------
    dbg = False  # True si on souhaite faire un IMPR_CO intermediaire, False sinon

    # Recuperation parametres solveur lineaire
    dSolveur = SOLVEUR[0].cree_dict_valeurs(SOLVEUR[0].mc_liste)
    for i in list(dSolveur.keys()):
        if dSolveur[i] is None:
            del dSolveur[i]
    if "TYPE_RESU" in dSolveur:  # because TYPE_RESU is a keyword with a 'global' position
        del dSolveur["TYPE_RESU"]
    if "OPTION" in dSolveur:  # because OPTION is a keyword with a 'global' position
        del dSolveur["OPTION"]
    if "FREQ" in dSolveur:  # because FREQ can be a keyword with a 'global' position
        del dSolveur["FREQ"]
    solveur_lineaire = dSolveur.get("METHODE").strip()
    dSolveur_infomode = dSolveur
    # pour INFO_MODE, le mot-clé facteur SOLVEUR ne doit pas contenir les
    # mot-clés POSTTRAITEMENTS et RESI_RELA
    if "POSTTRAITEMENTS" in dSolveur_infomode:
        del dSolveur_infomode["POSTTRAITEMENTS"]
    if "RESI_RELA" in dSolveur_infomode:
        del dSolveur_infomode["RESI_RELA"]

    # Rang du processus MPI et taille du MPI_COMM_WORLD
    # Lorsqu'on ne veut q'un niveau de parallelisme (celui du solveur lineaire)
    # on bluffe l'algo en posant rang=0/nbproc=1 pour tous les procs.
    # Cependant si le solveur est autre que MUMPS on s'arrete en erreur.
    if CALC_FREQ["NIVEAU_PARALLELISME"] == "COMPLET":
        rang = MPI.ASTER_COMM_WORLD.Get_rank()
        nbproc = MPI.ASTER_COMM_WORLD.Get_size()
        niv_par = "COMPLET"
    elif CALC_FREQ["NIVEAU_PARALLELISME"] == "PARTIEL":
        rang = 0
        nbproc = 1
        niv_par = "PARTIEL"
        nbproc_real = MPI.ASTER_COMM_WORLD.Get_size()
        if (nbproc_real > 1) & (solveur_lineaire != "MUMPS"):
            aster.affiche("MESSAGE", 72 * "-")
            UTMESS("F", "MODAL_14", vali=nbproc_real, valk=solveur_lineaire)
            aster.affiche("MESSAGE", 72 * "-")
    else:
        assert False  # Pb parametrage NIVEAU_PARALLELISME

    # Construction de la liste de frequences
    if CALC_FREQ["FREQ"]:
        lborne = []
        nnfreq = len(CALC_FREQ["FREQ"])
        for i in range(0, nnfreq):
            lborne.append(CALC_FREQ["FREQ"][i])
    else:
        assert False  # Pb parametrage CALC_FREQ

    # ----------------------------------------------------------------------
    #
    # 1b. TRAVAUX PREPARATOIRES AU LANCEMENT DE LA BOUCLE
    #
    # ----------------------------------------------------------------------

    # Modification de la sd_partition pour rendre compatible les niveaux de
    # parallelisme: celui en frequence et celui des calculs ele/assemblage.
    # Pour l'instant on fonctionne en mode 'CENTRALISE' au sein de la macro
    # (ds les sous-communicateurs, les calculs elem/assemblages
    # ne sont pas parallelises).
    # On remettra le mode de fonctionnement initial en fin de Macro.
    if nbproc > 1:
        old_prtk1 = recup_modele_partition(MATR_RIGI)
        sd_modele = None
        if MATR_RIGI is not None:
            sd_modele = MATR_RIGI.getDOFNumbering().getModel()
        if sd_modele is None:
            assert False  # Pb, on arrive pas a recuperer le nom du modele
        if old_prtk1 is not None:
            motdimo = {}
            motdimo["reuse"] = sd_modele
            motdimo["MODELE"] = sd_modele
            motdimo["DISTRIBUTION"] = _F(METHODE="CENTRALISE")
            __modimo = MODI_MODELE(**motdimo)

    # INFO_MODE global sur la liste de frequences
    # Parametres basiques
    motfaci = {}
    motfaci["COMPTAGE"] = _F(
        METHODE="AUTO",
        SEUIL_FREQ=CALC_FREQ["SEUIL_FREQ"],
        NMAX_ITER_SHIFT=CALC_FREQ["NMAX_ITER_SHIFT"],
        PREC_SHIFT=CALC_FREQ["PREC_SHIFT"],
    )

    # Gestion des frequences
    gestion_frequence(solveur_lineaire, nnfreq, nbproc)

    # Parametrage du parallelisme pour la couche FORTRAN/MPI
    if nbproc > 1:
        motfaci["PARALLELISME_MACRO"] = _F(TYPE_COM=1)
    __nbmodi = INFO_MODE(
        MATR_RIGI=MATR_RIGI,
        MATR_MASS=MATR_MASS,
        INFO=INFO,
        FREQ=lborne,
        NIVEAU_PARALLELISME=niv_par,
        SOLVEUR=dSolveur_infomode,
        **motfaci
    )

    # Gestion des sous-bandes de frequences et construction (si //) de l'objet
    nbmodeth, nbsb_nonvide, proc_sb_nvide = gestion_sous_bande(
        solveur_lineaire, __nbmodi, nnfreq, nbproc, lborne, STOP_BANDE_VIDE == "OUI"
    )

    # ----------------------------------------------------------------------
    #
    # 2. BOUCLE DE MODE_ITER_SIMULT(OPTION=BANDE) + NORM_MODE/IMPR_RESU
    #
    # ----------------------------------------------------------------------
    # On ne traite pas les sous-bandes vides (gain de temps et d'affichage!)
    # On affiche un message a la place.
    # Tous les processeurs font la meme chose pour ces sous-bandes vides, par
    # contre, ils se repartissent bien sur les non-vides.
    #
    # RQ MPI:
    # Les procs s'attendent via le MPI_COMM_SPLIT opere en fin de l'INFO_MODE
    # precedent. Ce n'est pas la peine de commencer la boucle si un proc
    # fait defaut.
    # ----------------------------------------------------------------------
    freq_ini = 1.0e99
    freq_fin = -1.0e99
    motscles = {}
    motscles["FILTRE_MODE"] = []
    numsb_nonvide = 0
    for i in range(0, nnfreq - 1):
        # --------------------------------------------------------------------
        #
        # 2. SOUS-BANDE NON VIDE
        #
        # --------------------------------------------------------------------
        sb_vide = None
        do_loop = None
        if __nbmodi["NB_MODE", i + 1] == 0:
            sb_vide = True
            do_loop = False
        else:
            numsb_nonvide = numsb_nonvide + 1
            sb_vide = False
            if ((nbproc > 1) & (proc_sb_nvide[rang] == numsb_nonvide)) | (nbproc == 1):
                do_loop = True
            else:
                do_loop = False

        if (not sb_vide) & do_loop:
            motscit = {}
            if CARA_ELEM is not None:
                motscit["CARA_ELEM"] = CARA_ELEM
            if CHAM_MATER is not None:
                motscit["CHAM_MATER"] = CHAM_MATER
            motscfa = {}
            if SOLVEUR_MODAL["DIM_SOUS_ESPACE"]:
                motscfa["DIM_SOUS_ESPACE"] = SOLVEUR_MODAL["DIM_SOUS_ESPACE"]
            if SOLVEUR_MODAL["COEF_DIM_ESPACE"]:
                motscfa["COEF_DIM_ESPACE"] = SOLVEUR_MODAL["COEF_DIM_ESPACE"]
            motscfa["FREQ"] = (lborne[i], lborne[i + 1])
            motscfa["TABLE_FREQ"] = __nbmodi
            motscit["CALC_FREQ"] = _F(
                OPTION="BANDE",
                SEUIL_FREQ=CALC_FREQ["SEUIL_FREQ"],
                NMAX_ITER_SHIFT=CALC_FREQ["NMAX_ITER_SHIFT"],
                PREC_SHIFT=CALC_FREQ["PREC_SHIFT"],
                **motscfa
            )

            if sturm == "LOCAL":
                motveri = "OUI"
            else:
                motveri = "NON"

            motscit["VERI_MODE"] = _F(
                STOP_ERREUR=stop_erreur,
                SEUIL=VERI_MODE["SEUIL"],
                STURM=motveri,
                PREC_SHIFT=VERI_MODE["PREC_SHIFT"],
            )

            motscit["STOP_BANDE_VIDE"] = STOP_BANDE_VIDE

            OPTION = "SANS"  # option for detecting rigid body modes
            if METHODE == "TRI_DIAG":
                if "NMAX_ITER_ORTHO" in args:
                    motscit["NMAX_ITER_ORTHO"] = SOLVEUR_MODAL["NMAX_ITER_ORTHO"]
                if "PREC_ORTHO" in args:
                    motscit["PREC_ORTHO"] = SOLVEUR_MODAL["PREC_ORTHO"]
                if "PREC_LANCZOS" in args:
                    motscit["PREC_LANCZOS"] = SOLVEUR_MODAL["PREC_LANCZOS"]
                if "MAX_ITER_QR" in args:
                    motscit["NMAX_ITER_QR"] = SOLVEUR_MODAL["NMAX_ITER_QR"]
                    if SOLVEUR_MODAL["MODE_RIGIDE"] == "OUI":
                        OPTION = "MODE_RIGIDE"
            elif METHODE == "JACOBI":
                if "NMAX_ITER_BATHE" in args:
                    motscit["NMAX_ITER_BATHE"] = SOLVEUR_MODAL["NMAX_ITER_BATHE"]
                if "PREC_BATHE" in args:
                    motscit["PREC_BATHE"] = SOLVEUR_MODAL["PREC_BATHE"]
                if "NMAX_ITER_JACOBI" in args:
                    motscit["NMAX_ITER_JACOBI"] = SOLVEUR_MODAL["NMAX_ITER_JACOBI"]
                if "PREC_JACOBI" in args:
                    motscit["PREC_JACOBI"] = SOLVEUR_MODAL["PREC_JACOBI"]
            elif METHODE == "SORENSEN":
                if "NMAX_ITER_SOREN" in args:
                    motscit["NMAX_ITER_SOREN"] = SOLVEUR_MODAL["NMAX_ITER_SOREN"]
                if "PARA_ORTHO_SOREN" in args:
                    motscit["PARA_ORTHO_SOREN"] = SOLVEUR_MODAL["PARA_ORTHO_SOREN"]
                if "PREC_SOREN" in args:
                    motscit["PREC_SOREN"] = SOLVEUR_MODAL["PREC_SOREN"]
            elif METHODE == "QZ":
                motscit["TYPE_QZ"] = SOLVEUR_MODAL["TYPE_QZ"]
            else:
                assert False  # Pb parametrage METHODE

            # Parametrage du parallelisme pour la couche FORTRAN/MPI
            # afin d'operer les comm des modes propres et des parametres
            # modaux.
            if nbproc > 1:
                motscit["PARALLELISME_MACRO"] = _F(
                    TYPE_COM=1, IPARA1_COM=numsb_nonvide, IPARA2_COM=nbsb_nonvide
                )

            # -----------------------------------------------------------------
            #
            # 2a. Calcul des modes
            #
            # -----------------------------------------------------------------                                      )
            __nomre0 = MODE_ITER_SIMULT(
                MATR_RIGI=MATR_RIGI,
                MATR_MASS=MATR_MASS,
                INFO=INFO,
                METHODE=METHODE,
                OPTION=OPTION,
                SOLVEUR=dSolveur,
                **motscit
            )

            # -----------------------------------------------------------------
            #
            # 2b. Preparation du test de Sturm de l'etape 3a.
            #
            # -----------------------------------------------------------------
            if sturm in (
                "GLOBAL",
                "OUI",
            ):  # in the case of CALC_MODES on several bands, OUI is reset to GLOBAL
                dicomode = {}
                dicomode = __nomre0.LIST_VARI_ACCES()
                if len(dicomode["FREQ"]) != 0:
                    raux_ini = dicomode["FREQ"][0]
                    raux_fin = dicomode["FREQ"][-1]
                    if raux_ini < freq_ini:
                        freq_ini = raux_ini
                    if raux_fin > freq_fin:
                        freq_fin = raux_fin
                else:
                    assert False  # Ce test ne marche pas (PB LIST_VARI_ACCES)

            # -----------------------------------------------------------------
            #
            # 2c. Préparation pour la concaténation de l'étape 3b.
            #
            # -----------------------------------------------------------------
            motscles["FILTRE_MODE"].append(_F(MODE=__nomre0, TOUT_ORDRE="OUI"))

            # Pour verification
            if dbg:
                IMPR_CO(CONCEPT=_F(NOM=__nomre0))

        # --------------------------------------------------------------------
        #
        # 2. SOUS-BANDE VIDE
        #
        # --------------------------------------------------------------------
        elif sb_vide:
            if STOP_BANDE_VIDE == "OUI":
                UTMESS("F", "MODAL_6", vali=(i + 1,))
            elif STOP_BANDE_VIDE == "NON":
                aster.affiche("MESSAGE", 72 * "-")
                UTMESS("I", "MODAL_6", vali=(i + 1,))
                aster.affiche("MESSAGE", 72 * "-")
            else:
                assert False  # Pb parametrage STOP_BANDE_VIDE

    # -----------------------------------------------------------------------
    #
    # 3a. Test de sturm effectif
    #
    # -----------------------------------------------------------------------
    if sturm == "NON":
        aster.affiche("MESSAGE", 72 * "-")
        UTMESS("I", "MODAL_2")
        aster.affiche("MESSAGE", 72 * "-")

    elif sturm == "LOCAL":
        aster.affiche("MESSAGE", 72 * "-")
        UTMESS("I", "MODAL_3")
        aster.affiche("MESSAGE", 72 * "-")

    elif sturm in (
        "GLOBAL",
        "OUI",
    ):  # in the case of CALC_MODES on several bands, OUI is reset to GLOBAL
        # Construction des 2 bornes de la bande a tester
        if nbmodeth != 0:
            omecor = CALC_FREQ["SEUIL_FREQ"]
            precshift = VERI_MODE["PREC_SHIFT"]
            freq_ini = (1.0 - precshift) * freq_ini
            freq_fin = (1.0 + precshift) * freq_fin
            if abs(freq_ini) < omecor:
                freq_ini = -omecor
            if abs(freq_fin) < omecor:
                freq_fin = omecor

            # Parametrage du parallelisme pour la couche FORTRAN/MPI
            if nbproc > 1:
                motfaci["PARALLELISME_MACRO"] = _F(TYPE_COM=2)
            __nbmodf = INFO_MODE(
                MATR_RIGI=MATR_RIGI,
                MATR_MASS=MATR_MASS,
                INFO=INFO,
                FREQ=(freq_ini, freq_fin),
                NIVEAU_PARALLELISME=niv_par,
                SOLVEUR=dSolveur_infomode,
                **motfaci
            )

            # Recuperation du nbre de modes donnes par STURM global
            nbmodesg = __nbmodf["NB_MODE", 1]
            if nbmodeth == nbmodesg:
                aster.affiche("MESSAGE", 72 * "-")
                UTMESS("I", "MODAL_4", valr=(freq_ini, freq_fin), vali=(nbmodesg, nbmodeth))
                aster.affiche("MESSAGE", 72 * "-")
            else:
                # Message similaire a ALGELINE5_24 pour le FORTRAN
                aorf = "A" if stop_erreur == "NON" else "F"
                valargs = _F(valr=(freq_ini, freq_fin, precshift), vali=(nbmodesg, nbmodeth))
                UTMESS(aorf, "MODAL_5", **valargs)

        # La bande globale est vide
        else:
            aster.affiche("MESSAGE", 72 * "-")
            UTMESS("I", "MODAL_7", valr=(lborne[0], lborne[nnfreq - 1]))
            aster.affiche("MESSAGE", 72 * "-")
    else:
        assert False  # Pb parametrage STURM

    # -----------------------------------------------------------------------
    #
    # 3b. Concaténation des résultats
    #
    # -----------------------------------------------------------------------

    modes = EXTR_MODE(**motscles)

    # ----------------------------------------------------------------------
    #
    # 3c. Epilogue
    #
    # RQ MPI:
    # Si la SD_PARTITION existait avant de rentrer ds la Macro on remet ses
    # anciennes valeurs pour continuer la suite du fichier de commande.
    # ----------------------------------------------------------------------
    if nbproc > 1:
        if old_prtk1 is not None:
            motdimo = {}
            motdimo["reuse"] = sd_modele
            motdimo["MODELE"] = sd_modele
            motdimo["DISTRIBUTION"] = _F(METHODE=old_prtk1)
            __modimo = MODI_MODELE(**motdimo)

    return modes


# ----------------------
# Routines auxiliaires
# ----------------------
# Routine pour recuperer sd_modele + option de la sd_partition (si elle existe)
def recup_modele_partition(MATR_RIGI):
    if isinstance(MATR_RIGI, (GeneralizedAssemblyMatrixReal, GeneralizedAssemblyMatrixComplex)):
        UTMESS("F", "MODAL_18")

    model = MATR_RIGI.getDOFNumbering().getModel()

    if not model.existsPartition():
        return None
    else:
        return model.getPartitionMethod()


# Routine pour recuperer nbre de modes theorique (nbmodeth) determine par l'INFO_MODE
# prealable + le nbre de sous_bandes non vides (nbsd_nonvide) + msg d'erreurs associes.
# De plus, lorsque le mode parallele est active, elle gere les incompatibilites
# entre le nombre de sous-bandes, le nombre de processeurs et la methode utilisee
# pour le solveur lineaire.
# En mode parallele, elle remplit le vecteur proc_sb_nvide de taille le nbproc.
# proc_sb_nvide(i)=numero de la sbande non vide traitee par le proc de rang i.
# Regle 1: On garde les procs travaillant sur la meme sbande contigues.
# Regle 2: En cas de desequilibre de charge, on surcharge les premieres sbandes.
# RQ. Meme regle que pour la distribution frequence/proc operee ds le F77:
# cf. op0032.


def gestion_sous_bande(solveur_lineaire, __nbmodi, nnfreq, nbproc, lborne, stop):
    nbsb_nonvide = None
    proc_sb_nvide = []
    # Recuperation du nbre de modes total theorique
    nbmodeth = 0
    modemin = 1000000
    imin = 0
    modemax = -1
    imax = 0
    for i in range(0, nnfreq - 1):
        modei = __nbmodi["NB_MODE", i + 1]
        nbmodeth = nbmodeth + modei
        if modei > modemax:
            modemax = modei
            imax = i + 1
        if modei < modemin:
            modemin = modei
            imin = i + 1

    if (modemax > 3 * modemin) & (modemax > 50):
        UTMESS("I", "MODAL_19", vali=(imin, modemin, imax, modemax))
    if modemin < 10:
        UTMESS("I", "MODAL_20", vali=(imin, modemin))
    if modemax > 100:
        UTMESS("I", "MODAL_21", vali=(imax, modemax))

    # Recuperation du nbre de sous-bandes non vides
    nbsb_nonvide = 0
    for i in range(0, nnfreq - 1):
        if __nbmodi["NB_MODE", i + 1] != 0:
            nbsb_nonvide = nbsb_nonvide + 1

    if (nbmodeth == 0) | (nbsb_nonvide == 0):
        aster.affiche("MESSAGE", 72 * "-")
        if stop:
            UTMESS("F", "MODAL_8", valr=(lborne[0], lborne[nnfreq - 1]))
        else:
            UTMESS("A", "MODAL_8", valr=(lborne[0], lborne[nnfreq - 1]))
        aster.affiche("MESSAGE", 72 * "-")

    if nbproc > 1:
        if (nbproc < nbsb_nonvide) | ((nbproc > nbsb_nonvide) & (solveur_lineaire != "MUMPS")):
            aster.affiche("MESSAGE", 72 * "-")
            UTMESS("F", "MODAL_9", vali=(nbproc, nbsb_nonvide), valk=solveur_lineaire)
            aster.affiche("MESSAGE", 72 * "-")
        div = None
        reste = None
        div = nbproc // nbsb_nonvide
        reste = nbproc - nbsb_nonvide * div
        if (nbproc > nbsb_nonvide) & (reste != 0):
            aster.affiche("MESSAGE", 72 * "-")
            UTMESS("I", "MODAL_12", vali=(nbsb_nonvide, div, div + 1))
            aster.affiche("MESSAGE", 72 * "-")

        l1 = nbproc // nbsb_nonvide
        l11 = l1 + 1
        l2 = nbproc - (l1 * nbsb_nonvide)
        num_sb = 0
        for i in range(0, l2):
            num_sb = num_sb + 1
            for j in range(0, l11):
                proc_sb_nvide.append(num_sb)

        for i in range(l2, nbsb_nonvide):
            num_sb = num_sb + 1
            for j in range(0, l1):
                proc_sb_nvide.append(num_sb)

    else:
        proc_sb_nvide.append(-999)

    return nbmodeth, nbsb_nonvide, proc_sb_nvide


# Routine pour gerer, lorsque le mode parallele est active, les incompatibilites
# entre le nombre de frequences, le nombre de processeurs et la methode utilisee
# pour le solveur lineaire.


def gestion_frequence(solveur_lineaire, nnfreq, nbproc):
    if nbproc > 1:
        if (nbproc < nnfreq - 1) | ((nbproc > nnfreq - 1) & (solveur_lineaire != "MUMPS")):
            aster.affiche("MESSAGE", 72 * "-")
            UTMESS("F", "MODAL_10", vali=(nbproc, nnfreq - 1), valk=solveur_lineaire)
            aster.affiche("MESSAGE", 72 * "-")
        div = None
        reste = None
        div = nbproc // (nnfreq - 1)
        reste = nbproc - (nnfreq - 1) * div
        if (nbproc > nnfreq - 1) & (reste != 0):
            aster.affiche("MESSAGE", 72 * "-")
            UTMESS("I", "MODAL_11", vali=(nnfreq - 1, div, div + 1))
            aster.affiche("MESSAGE", 72 * "-")

    return
