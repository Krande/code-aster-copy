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

import os
import re
import sys
import traceback

import aster_core
import libaster

from ..Utilities import (
    ExecutionParameter,
    Options,
    Singleton,
    _,
    center,
    config,
    convert,
    copy_text_to,
    cut_long_lines,
    force_list,
    textbox,
    to_unicode,
    ufmt,
)
from ..Utilities.mpi_utils import MPI

DEBUG = False
CENTER = 1
DECORATED = 2
ALL_UNIT = 4

MAXLENGTH = 100
LINE_WITH = 80

contacter_assistance = _(
    """
Il y a probablement une erreur dans la programmation.
Veuillez contacter votre assistance technique."""
)

# Exceptions ids - keep consistency with astercxx.h/asterf.h
EXCNUM = {
    1: ("AsterError", libaster.AsterError),
    2: ("ConvergenceError", libaster.ConvergenceError),
    3: ("IntegrationError", libaster.IntegrationError),
    4: ("SolverError", libaster.SolverError),
    5: ("ContactError", libaster.ContactError),
    6: ("TimeLimitError", libaster.TimeLimitError),
}

# voir en fin de fin les faux appels à UTMESS pour la vérification des messages


def list_unit(code):
    """Retourne la liste des noms de fichiers (logiques) sur lesquels doit
    etre imprimé le message.
    """
    # 'D' pour afficher un diagnostic 'F' sans les effets
    # 'Z' levée d'exception
    d = {"E": ("MESSAGE",), "I": ("MESSAGE",), "A": ("MESSAGE",)}
    d["F"] = d["S"] = d["Z"] = d["D"] = d["E"]
    d["X"] = d["A"]
    return d.get(code, d["F"])


class MESSAGE_LOGGER(metaclass=Singleton):
    """Classe gérant l'impression de messages.
    On ne crée qu'une instance de ce type (singleton).
    Cette instance est accessible dans le module `aster_core` pour les appels
    depuis le fortran.
    """

    def __init__(self):
        """Initialisation"""
        self.init_buffer()

        # est-ce qu'une erreur <E> s'est produite
        self.erreur_E = False

        # compteur des alarmes émises { 'id_alarm' : count }
        self.nmax_alarm = 5
        self.count_alarm = {}  # dans la commande courante (pour arret à 5)
        self.count_alarm_tot = {}  # au total

        # alarmes à ignorer, à masquer (on ne les compte pas temporairement)
        self._ignored_alarm = {}
        self._hidden_alarm = {}

        # on prépare le dictionnaire des valeurs par défaut des arguments
        # (dicarg) :
        self.default_args = {}
        # initialisation des 50 premiers
        for i in range(1, 51):
            self.default_args["i%d" % i] = 99999999
            self.default_args["r%d" % i] = 9.9999e99
            self.default_args["k%d" % i] = "xxxxxx"
        # mettre en cache les messages 'I' (et uniquement 'I')
        self._cache_txt = {}
        # arguments mpi : ligne de commande à envoyer au proc #0
        self._mpi_rank = None
        self._mpi_nbcpu = None
        self.init_mpi_error()

    def init_mpi_error(self):
        """Stocke les informations nécessaires pour la gestion des erreurs en MPI."""
        if not config["ASTER_HAVE_MPI"]:
            return
        if not MPI.Is_initialized():
            return
        self._mpi_rank = MPI.ASTER_COMM_WORLD.Get_rank()
        self._mpi_nbcpu = MPI.ASTER_COMM_WORLD.Get_size()

    def __call__(self, *args, **kwargs):
        """Raccourci pour simplifier l'appel depuis astermodule.c et UTMESS."""
        self.print_message(*args, **kwargs)

    def print_message(
        self,
        code,
        idmess,
        valk=(),
        vali=(),
        valr=(),
        exc_typ=None,
        exception=False,
        print_as=None,
        files=None,
        cc=True,
    ):
        """Appelé par la routine fortran U2MESG ou à la fonction python UTMESS
        pour afficher un message.
        L'impression de ce message est différée si le `code` est suivi d'un "+".
        - code  : 'A', 'E', 'S', 'F', 'I'
        - idmess : identificateur du message
        - valk, vali, valr : liste des chaines, entiers ou réels.
        Si exception==True, on lève une exception en cas d'erreur, sinon
        c'est l'appelant qui devra s'en charger (dans le C a priori).
        'print_as', 'files', 'cc' : cf. print_buffer_content.
        """
        idmess = idmess.strip()
        # le '+' n'a pas de sens pour les messages 'I'.
        if code == "I+":
            code = "I"
        dictmess = {}
        if code == "I":
            cached = self._cache_txt.get(idmess)
            if cached:
                flags, msg = cached
                dicarg = self.build_dict_args(valk, vali, valr)
                try:
                    dictmess = {
                        "code": code,
                        "flags": flags,
                        "id_message": idmess,
                        "corps_message": ufmt(msg, dicarg),
                        "context_info": "",  # self.get_context(ctxt_msg, idmess, dicarg),
                    }
                except Exception:
                    # in case of formatting error, ignore the cached value
                    pass
        if code == "A" and ExecutionParameter().option & Options.WarningAsError:
            code = "F"
            self.print_message("E", "CATAMESS_3", valk=[idmess])
        if self._parent is None:
            self.set_parent(idmess)
        if not self.update_counter(code, idmess):
            return
        # récupération du texte du message
        if not dictmess:
            dictmess = self.get_message(code, idmess, valk, vali, valr, exc_typ)

        # on le met dans le buffer
        self.add_to_buffer(dictmess)

        # si on n'attend pas une suite, ...
        if is_last_message(code):
            id0 = self.get_current_id().strip()
            # vérification des compteurs
            self.check_limit()

            # on imprime le message en attente
            self.print_buffer_content(print_as, files, cc)

            if exception and code[0] in ("S", "F"):
                if self._mpi_rank is not None:
                    aster_core.MPI_Warn()
                exc_typ = dictmess.get("exc_typ")
                if exc_typ:
                    raise exc_typ(id0, valk, vali, valr)
                raise libaster.AsterError(id0, valk, vali, valr)
        return None

    def build_dict_args(self, valk, vali, valr):
        """Construit le dictionnaire de formatage du message."""
        # homogénéisation : uniquement des tuples + strip des chaines de
        # caractères
        valk, vali, valr = list(map(force_list, (valk, vali, valr)))
        valk = [k.strip() for k in valk]

        # variables passées au message
        dicarg = self.default_args.copy()
        for i in range(1, len(valk) + 1):
            dicarg["k%d" % i] = to_unicode(valk[i - 1])
        for i in range(1, len(vali) + 1):
            dicarg["i%d" % i] = vali[i - 1]
        for i in range(1, len(valr) + 1):
            dicarg["r%d" % i] = valr[i - 1]
        # valeur spéciale : ktout = concaténation de toutes les chaines
        dicarg["ktout"] = " ".join(valk)

        return dicarg

    def get_message(self, code, idmess, valk=(), vali=(), valr=(), exc_typ=None):
        """Retourne le texte du message dans un dictionnaire dont les clés sont :
        'code', 'id_message', 'corps_message'
        """
        # décodage : idmess => (catamess, numess)
        idmess = idmess.strip()
        x = idmess.split("_")
        assert len(x) > 1, idmess
        catamess = "_".join(x[0:-1]).lower()
        numess = int(x[-1])
        assert numess > 0 and numess < 100, idmess

        # import catamess => cata_msg
        try:
            mod = __import__("code_aster.Messages.%s" % catamess, globals(), locals(), [catamess])
            # si le dictionnaire n'existe pas, on alertera au moment du
            # formatage.
            cata_msg = getattr(mod, "cata_msg", {})
        except Exception as msg:
            # doit permettre d'éviter la récursivité (catamess réservé à
            # Utmess)
            if catamess != "catamess":
                code = "A"
                self.print_message(code, "CATAMESS_57", valk=(catamess, str(msg)))
            cata_msg = {}

        # corps du message
        fmt_msg = "?"
        try:
            dicarg = self.build_dict_args(valk, vali, valr)

            # cata_msg[num] = 'format'
            #              ou { 'message' : 'format',
            #                   'flags' : 'DECORATED | CENTER',
            #                   'context' : 'éléments de contexte' }
            if isinstance(cata_msg[numess], dict):
                fmt_msg = cata_msg[numess]["message"]
                flags = eval(cata_msg[numess].get("flags", 0))
            else:
                fmt_msg = cata_msg[numess]
                flags = 0

            dictmess = {
                "code": code,
                "flags": flags,
                "id_message": idmess,
                "corps_message": ufmt(fmt_msg, dicarg),
                "context_info": "",  # self.get_context(ctxt_msg, idmess, dicarg),
            }
            if code == "I":
                self._cache_txt[idmess] = flags, convert(fmt_msg)
        except Exception as msg:
            if code == "I":
                code = "A"
            dictmess = {
                "code": code,
                "flags": 0,
                "id_message": idmess,
                "corps_message": _(
                    """Erreur de programmation.
Le message %s n'a pas pu être formaté correctement.
Arguments :
    entiers : %s
    réels   : %s
    chaines : %s

    format  : %s
--------------------------------------------------------------------------
%s
Exception : %s
--------------------------------------------------------------------------

%s"""
                ),
                "context_info": "",
            }
            args = (
                idmess,
                vali,
                valr,
                valk,
                fmt_msg,
                "".join(traceback.format_exc()),
                msg,
                contacter_assistance,
            )
            # cette étape ne doit jamais faire planter !
            try:
                dictmess["corps_message"] = dictmess["corps_message"] % args
            except Exception:
                dictmess["corps_message"] = repr(args)
        # type d'exception
        if exc_typ:
            if isinstance(exc_typ, int):
                exc_args = EXCNUM.get(exc_typ, (None, RuntimeError))
            else:
                # exc_typ is an Exception
                exc_args = exc_typ.__class__.__name__, exc_typ
            dictmess["exc_name"], dictmess["exc_typ"] = exc_args
        return dictmess

    def GetText(self, *args, **kwargs):
        """Retourne le texte du message pret a etre imprime."""
        return self.format_message(self.get_message(*args, **kwargs))

    def init_buffer(self):
        """Initialise le buffer."""
        self._buffer = []
        self.set_parent(None)

    def is_buffer_empty(self):
        """Tell if the buffer is currently empty"""
        return len(self._buffer) < 1

    def set_parent(self, idmess):
        """Store the parent id of the current message"""
        self._parent = idmess

    def add_to_buffer(self, dictmess):
        """Ajoute le message décrit dans le buffer en vue d'une impression
        ultérieure.
        """
        self._buffer.append(dictmess)

    def get_current_code(self):
        """Retourne le code du message du buffer = code du message le plus grave
        (cf. dgrav)
        """
        dgrav = {"?": -9, "I": 0, "A": 1, "S": 4, "Z": 4, "E": 6, "D": 9, "F": 10}

        current = "?"
        exc_name = None
        exc_typ = None
        for dictmess in self._buffer:
            code = dictmess["code"][0]
            if dgrav.get(code, -9) > dgrav.get(current, -9):
                current = code
            exc_name = exc_name or dictmess.get("exc_name")
            exc_typ = exc_typ or dictmess.get("exc_typ")
        return current, exc_name, exc_typ

    def get_current_flags(self):
        """Retourne les flags du message du buffer = flags du premier."""
        return self._buffer[0]["flags"]

    def get_current_id(self):
        """Retourne l'id du message du buffer = id du premier message"""
        return self._buffer[0]["id_message"]

    def print_buffer_content(self, print_as=None, files=None, cc=True):
        """Extrait l'ensemble des messages du buffer dans un dictionnaire unique,
        imprime le message, et vide le buffer pour le message suivant.
        - code : celui du message le plus grave (cf. dgrav)
        - id   : celui du premier message qui est affiché
        - corps : concaténation de tous les messages.

        'print'_as permet d'imprimer un message sur des fichiers autres que les fichiers
        habituels de 'code'. Par ex, imprimer un message d'info sur 'ERREUR'.
        'files' : liste de noms de fichiers ou objets fichier dans lesquels
        écrire le message
        'cc' : si True, on écrit comme d'habitude et dans les 'files',
        si False, on n'écrit que sur les fichiers habituels (MESSAGE, RESULTAT,
        ERREUR) ou bien dans 'files' si fournit.
        """
        if isinstance(files, str):
            files = files.strip()
        if len(self._buffer) < 1:
            return None

        # construction du dictionnaire du message global
        dglob = {
            "flags": self.get_current_flags(),
            "id_message": self.get_current_id(),
            "liste_message": [],
            "liste_context": [],
        }
        args = self.get_current_code()
        dglob["code"], dglob["exc_name"], dglob["exc_typ"] = args
        for dictmess in self._buffer:
            dglob["liste_message"].append(dictmess["corps_message"])
            dglob["liste_context"].append(dictmess["context_info"])
        dglob["corps_message"] = "".join(dglob["liste_message"])
        dglob["context_info"] = "".join(dglob["liste_context"])

        # liste des unités d'impression en fonction du type de message
        if dglob["flags"] & ALL_UNIT:
            print_as = "E"
        l_unit = list_unit(print_as or dglob["code"])

        # texte final et impression
        if cc or not files:
            txt = self.format_message(dglob)
            for unite in l_unit:
                self.affiche(unite, txt)
        # "files"
        if files:
            copy_text_to(convert(txt), files)

        self.init_buffer()

    def disable_alarm(self, idmess, hide=False):
        """Ignore l'alarme "idmess"."""
        idmess = idmess.strip()
        if hide:
            self._hidden_alarm[idmess] = self._hidden_alarm.get(idmess, 0) + 1
        else:
            self._ignored_alarm[idmess] = self._ignored_alarm.get(idmess, 0) + 1

    def reset_alarm(self, idmess, hide=False):
        """Réactive l'alarme "idmess"."""
        idmess = idmess.strip()
        if hide:
            self._hidden_alarm[idmess] = min(self._hidden_alarm.get(idmess, 0) - 1, 0)
        else:
            self._ignored_alarm[idmess] = min(self._ignored_alarm.get(idmess, 0) - 1, 0)

    def is_alarm_disabled(self, idmess):
        """Doit-on ignorer l'alarme "idmess" ?"""
        return self._ignored_alarm.get(idmess, 0) + self._hidden_alarm.get(idmess, 0) > 0

    def get_info_alarm(self, only_ignored=False):
        """Retourne la liste des alarmes émises, le nombre d'occurrence
        pour chacune d'elle et un indicateur disant si elle a été masquée ou pas."""
        s_alarm = set(self._ignored_alarm.keys())
        if not only_ignored:
            s_alarm.update(list(self.count_alarm_tot.keys()))
        l_all = list(s_alarm)
        l_all.sort()
        # occurrences
        l_alarm, l_occ, l_masq = [], [], []
        for idmess in l_all:
            nb = self.count_alarm_tot.get(idmess, 0)
            if nb > 0:
                l_alarm.append(idmess)
                l_occ.append(nb)
                l_masq.append(int(self._ignored_alarm.get(idmess) is not None))
        return list(zip(l_alarm, l_occ, l_masq))

    def get_info_alarm_nb(self, only_ignored=False):
        """Retourne le nombre d'alarme émises (et non masquées)."""
        res = self.get_info_alarm(only_ignored)
        res = [item for item in res if item[2] == 0]
        return len(res)

    def info_alarm(self, only_ignored=False):
        """Fournit les infos sur les alarmes activées."""
        # on sépare des éventuels messages en attente
        self.print_buffer_content()
        # entete
        dictmess = self.get_message("I", "CATAMESS_89")
        self.add_to_buffer(dictmess)
        # occurrences
        res = self.get_info_alarm(only_ignored)
        not_seen = set(self._ignored_alarm.keys())
        for idmess, nb, masq in res:
            mark = " " + "(*)" * masq
            dictmess = self.get_message("I", "CATAMESS_90", valk=(mark, idmess), vali=nb)
            self.add_to_buffer(dictmess)
            not_seen.discard(idmess)
        if not res:
            dictmess = self.get_message("I", "CATAMESS_92")
            self.add_to_buffer(dictmess)
        self.print_buffer_content()
        if not_seen:
            code = "A"
            if self._mpi_rank == 0:
                self.print_message(code, "CATAMESS_87", valk=list(not_seen), exception=True)

    def update_counter(self, code, idmess):
        """Update the counters of alarms.
        The counter is updated only for the first message in the buffer. So
        it is important to call this method before adding the message into the
        buffer.
        Return True if everything is ok, False if the message will be skipped."""
        if code[0] == "A":
            parent = self._parent
            if (
                parent == idmess
                and self._hidden_alarm.get(idmess, 0) == 0
                and self.is_buffer_empty()
            ):
                self.count_alarm[idmess] = self.count_alarm.get(idmess, 0) + 1
                self.count_alarm_tot[idmess] = self.count_alarm_tot.get(idmess, 0) + 1
            if self.is_alarm_disabled(parent) or self.count_alarm.get(parent, 0) > self.nmax_alarm:
                # ignorer l'alarme ou count_alarm > max, on passe
                if is_last_message(code):
                    self.set_parent(None)
                return False
        return True

    def check_limit(self):
        """Vérifications des compteurs et réaction si besoin."""
        code = self.get_current_code()[0]
        idmess = self.get_current_id().strip()
        if code == "E":
            self.erreur_E = True
        elif code == "F":
            self.erreur_E = False
        elif code == "A" and self.count_alarm.get(idmess, 0) == self.nmax_alarm:
            # Pour mettre en relief le message CATAMESS_41, on le sépare
            # de la dernière alarme
            self.print_buffer_content()
            dictmess = self.get_message(code, "CATAMESS_41", valk=idmess, vali=self.nmax_alarm)
            self.add_to_buffer(dictmess)

    def check_counter(self, info_alarm=0, silent=0):
        """Méthode "jusqu'ici tout va bien" ! (Interface C : chkmsg)
        Si des erreurs <E> se sont produites, on arrete le code en <F>.
        Appelée par FIN ou directement au cours de l'exécution d'une commande.
        Retourne un entier : 0 si tout est ok.
        Si silent==1, on n'émet pas de message, on ne s'arrete pas.
        """
        iret = 0
        if sys.is_finalizing():
            return iret
        if self.erreur_E:
            iret = 4
            self.erreur_E = False
            if not silent:
                self.print_message("F", "CATAMESS_6", exception=True)
        if info_alarm:
            self.info_alarm()
        return iret

    def reset_command(self):
        """Méthode appelée entre les commandes. (Interface C : resmsg)
        On remet à zéro le compteur d'alarme,
        on vérifie les erreurs <E> en attente."""
        self.check_counter()
        # reset des alarmes
        self.count_alarm = {}
        # reset du cache, sans doute inutile car l'ensemble des messages représente
        # environ 1 Mo.
        if len(self._cache_txt) > 1000:
            self._cache_txt = {}

    def format_message(self, dictmess):
        """Formate le message décrit dans un dico :
        'code'          : A, E, S, F, I
        'id_message'    : identification du message
        'corps_message' : texte
        """
        dcomm = {
            "A": _(
                """Ceci est une alarme. Si vous ne comprenez pas le sens de cette
alarme, vous pouvez obtenir des résultats inattendus !"""
            ),
            "E": _("""Cette erreur sera suivie d'une erreur fatale."""),
            "S": _(
                """Cette erreur est fatale. Le code s'arrête. Toutes les étapes
du calcul ont été sauvées dans la base jusqu'au moment de l'arret."""
            ),
            "F": _("""Cette erreur est fatale. Le code s'arrête."""),
        }

        dmsg = dictmess.copy()
        typmess = self.get_type_message(dictmess)
        comm = dcomm.get(typmess, "")
        if re.search("^DVP", dmsg["id_message"]):
            comm += contacter_assistance
        comm = comm.strip()

        if dmsg["flags"] & DECORATED or typmess != "I":
            text = ["<{0}> <{1}>".format(typmess, dmsg["id_message"]), ""]
        else:
            text = []
        text.append(dmsg["corps_message"].lstrip())
        if comm:
            text.extend(["", comm])

        body = os.linesep.join(text)
        if dmsg["flags"] & CENTER:
            body = center(body, MAXLENGTH)
        if dmsg["flags"] & DECORATED or typmess != "I":
            body = textbox(body, MAXLENGTH)
        else:
            body = cut_long_lines(body, MAXLENGTH)
        return body

    def get_type_message(self, dictmess):
        """Retourne le type du message affiché.
        En cas d'erreur, si on lève une exception au lieu de s'arreter,
        on affiche le type de l'erreur.
        """
        code = dictmess["code"]
        typmess = code.strip()
        if self.onFatalError().startswith("EXCEPTION"):
            if typmess in ("E", "F"):
                typmess = "EXCEPTION"
        if typmess == "D":
            typmess = "F"
        # dans tous les cas, pour S et Z (exception), on affiche EXCEPTION.
        elif code == "S":
            typmess = "EXCEPTION"
        elif code == "Z":
            typmess = dictmess.get("exc_name") or "EXCEPTION"
        return typmess

    # définitions pour fonctionner sans le module aster
    def affiche(self, unite, txt):
        """Affichage du message"""
        txt = convert(txt)
        libaster.affich(unite, txt)

    def onFatalError(self):
        """Récupérer le comportement en cas d'erreur fatale."""
        return libaster.onFatalError()


def is_last_message(code):
    """Tell if a message 'code' is the last message or not."""
    return len(code) < 2 or code[1] != "+"


# unique instance du MESSAGE_LOGGER
MessageLog = MESSAGE_LOGGER()


def UTMESS(
    code, idmess, valk=(), vali=(), valr=(), exc_typ=None, print_as=None, files=None, cc=True
):
    """Utilitaire analogue à la routine fortran U2MESS/U2MESG avec les arguments
    optionnels.

    Remarques :
    - Nommer les arguments permet de ne pas tous les passer.
    - Meme fonctionnement que U2MESG :

      - appel à MessageLog
      - puis exception ou abort en fonction du niveau d'erreur.

    Arguments:
        code (str): 'A', 'E', 'S', 'F', 'I'
        idmess (str): identificateur du message
        valk, vali, valr : liste des chaines, entiers ou réels.
    """
    MessageLog(
        code,
        idmess,
        valk,
        vali,
        valr,
        exc_typ=exc_typ,
        exception=True,
        print_as=print_as,
        files=files,
        cc=cc,
    )


def ASSERT(condition, message=""):
    """Remonter un assert dans un message."""
    if condition:
        return
    stack = traceback.format_stack(limit=10)
    UTMESS("F", "DVP_9", valk=["".join(stack[:-1]), message])


def message_exception(code, exc):
    """Retourne le message associé à une exception `AsterError`
    tel qu'il aurait été imprimé par UTMESS selon la valeur de
    `code` ('I', 'A', 'S', 'F', 'Z'...)."""
    # check order of args in Exception.cxx
    return MessageLog.GetText(code, *exc.args)


def MasquerAlarme(idmess):
    """Masque une alarme : ni affichee, ni comptee.
    Utilisation dans les macros :
    - MasquerAlarme(XXX)  au debut de la macro
    - RetablirAlarme(XXX) a la fin de la macro
    Comme il s'agit d'un compteur qui est incremente puis decremente, il est
    imperatif qu'il y ait autant d'appel a MasquerAlarme qu'a RetablirAlarme.
    """
    MessageLog.disable_alarm(idmess, hide=True)


def RetablirAlarme(idmess):
    """Retablit l'etat initial pour l'alarme 'idmess'."""
    MessageLog.reset_alarm(idmess, hide=True)


# faux appels à UTMESS
def __fake__():
    UTMESS("I", "SUPERVIS_96")  # émis depuis le C (inisig)
    UTMESS("I", "SUPERVIS_99")  # émis par le logger
    UTMESS("I", "JEVEUX_44")  # émis depuis le C (iodr)
    UTMESS("I", "JEVEUX_45")  # émis depuis le C (iodr)
    UTMESS("I", "GENERIC_1")  # dans des tests pour traiter les exceptions
    UTMESS("I", "CATAMESS_1")  # pour supv002a
    UTMESS("I", "CATAMESS_2")  # pour vocab01a
    UTMESS("I", "CATAMESS_55")  # pour u2mesg.f via UTPRIN
    UTMESS("I", "CATAMESS_69")  # pour u2mesg.f via UTPRIN
    UTMESS("I", "CATAMESS_70")  # pour u2mesg.f via UTPRIN
    # utilisé ici
    UTMESS("I", "CATAMESS_3")
    UTMESS("I", "CATAMESS_6")
    UTMESS("I", "CATAMESS_41")
    UTMESS("I", "CATAMESS_57")
    UTMESS("I", "CATAMESS_87")
    UTMESS("I", "CATAMESS_89")
    UTMESS("I", "CATAMESS_90")
    UTMESS("I", "CATAMESS_92")
    # appelé par levé d'exception
    # dans TableReader.py
    UTMESS("I", "TABLE0_10")
    UTMESS("I", "TABLE0_11")
    UTMESS("I", "TABLE0_12")
    UTMESS("I", "TABLE0_13")
    UTMESS("I", "TABLE0_15")
    UTMESS("I", "TABLE0_43")
    # dans function_py.py
    UTMESS("I", "FONCT0_27")
    UTMESS("I", "FONCT0_28")
    UTMESS("I", "FONCT0_29")
    # en C++, LinearSolver
    UTMESS("I", "FACTOR_13")
