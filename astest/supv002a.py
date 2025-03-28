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

import gettext
import multiprocessing as MPR
import os
import os.path as osp
import queue
import re
import signal
import tempfile
import time
import traceback
from functools import partial
from glob import glob
from subprocess import PIPE, CalledProcessError, Popen, check_call

from code_aster import Messages
from code_aster.Messages import UTMESS, MessageLog
from code_aster.Utilities import ExecutionParameter, convert, is_int
from code_aster.Utilities import localization as LO

ENCODING = "utf-8"
VALUES = MessageLog.default_args.copy()
VALUES["ktout"] = "xxxxxx"

REI = re.compile(r"%\(i[0-9]+\)[\.0-9\-\+ ]*([a-zA-Z])", re.M)
RER = re.compile(r"%\(r[0-9]+\)[\.0-9\-\+ ]*([a-zA-Z])", re.M)
REK = re.compile(r"%\(k[0-9]+\)[\.0-9\-\+]*([a-zA-Z])", re.M)
RE1 = re.compile(r"%\((.[^0-9].*?)\)[\.0-9\-\+ ]*[a-zA-Z]", re.M)
RE2 = re.compile(r"%\(([^irk].*?)\)[\.0-9\-\+ ]*[a-zA-Z]", re.M)

RE_UNAUTH = [re.compile(r"([*#=\+\-!/\?<>&@]{4})", re.M)]

try:
    import aster

    from code_aster.Cata.Syntax import _F
    from code_aster.Commands import CREA_TABLE, TEST_TABLE

    loginfo = partial(aster.affiche, "MESSAGE")
except ImportError:

    def loginfo(msg):
        print(msg)


logdbg = None


class Checker:
    """A simple checker object."""

    def __init__(self, aspell):
        self.err = []
        self.wrn = []
        self.allwarns = []
        self.lang = None
        self.mod = None
        self.words = set()
        self.allwords = set()
        self.ignored = set()
        self.aspell = aspell

    def set_ignored(self, errors):
        self.ignored = set(errors)

    def set_current_lang(self, lang):
        self.lang = lang

    def set_current_mod(self, mod):
        self.mod = mod

    def info(self, msg):
        loginfo("<%s> %s: %s" % (self.lang, self.mod, msg))

    def error(self, error):
        self.err.append((self.lang, self.mod, error))

    def warning(self, error):
        self.wrn.append((self.lang, self.mod, error))

    def warning_spell(self, idmess, words):
        msg = "%s: unknown words %s" % (idmess, tuple(words))
        if set(words).difference(self.ignored):
            self.warning(msg)
            self.words.update(words)
        self.allwarns.append((self.lang, self.mod, msg))
        self.allwords.update(words)

    def count_errors(self):
        return len(self.err)

    def count_unknown_words(self):
        return len(self.get_unknown_words())

    def get_errors(self):
        txt = ["<%s> %s: %s" % err for err in self.err]
        return os.linesep.join(txt)

    def get_warnings(self):
        txt = ["<%s> %s: %s" % wrn for wrn in self.wrn]
        return os.linesep.join(txt)

    def get_all_warnings(self):
        txt = ["<%s> %s: %s" % wrn for wrn in self.allwarns]
        return os.linesep.join(txt)

    def get_unknown_words(self):
        lw = list(self.words)
        lw.sort()
        return lw

    def get_all_words(self):
        lw = list(self.allwords)
        lw.sort()
        return lw

    def get_modules(self):
        smod = set([mod for lang, mod, err in self.err])
        lmod = list(smod)
        lmod.sort()
        return lmod


def read_file(stream, queue):
    """Read on stream and put lines into queue"""
    while True:
        line = stream.readline().strip()
        logdbg and logdbg("aspell outputs: %r" % line)
        queue.put(line)


class AspellCall:
    """A pipe to call aspell"""

    @staticmethod
    def check_aspell():
        """Check that aspell is available.

        Returns:
            bool: *True* if aspell seems available, *False* otherwise.
        """
        try:
            check_call(["aspell", "--version"])
            isok = True
        except (CalledProcessError, FileNotFoundError):
            isok = False
        return isok

    def __init__(self, personal_dict, lang, encoding):
        """Open the pipe"""
        cmd = ["aspell", "pipe", "--encoding=%s" % encoding, "--lang=%s" % lang]
        if personal_dict and osp.exists(personal_dict):
            cmd.append("--personal=%s" % personal_dict)
        logdbg and logdbg("command: " + " ".join(cmd))
        self.lang = lang
        self.pipe = Popen(cmd, stdin=PIPE, stdout=PIPE, universal_newlines=True)
        self.inp = self.pipe.stdin
        self.out = self.pipe.stdout
        self.queue = MPR.Queue()
        self.rdr = MPR.Process(target=read_file, args=(self.out, self.queue))
        self.rdr.start()
        init = self.queue.get()
        assert "aspell" in init.lower(), "aspell probably failed to start: \n" "%s" % init
        # all characters except alphabetic ones are separators
        self._rxspl = re.compile(r"[ _0123456789\W]+", re.M | re.I | re.UNICODE)

    def __del__(self):
        """Close the pipe"""
        self.close()

    def close(self):
        """Close the pipe"""
        self.inp.close()
        self.out.close()
        self.rdr.terminate()
        self.pipe.terminate()

    def _clean(self, string):
        """clean 'string' before passing it to aspell"""
        if type(string) is not str:
            string = str(string, ENCODING)
        uwords = self._rxspl.split(string.strip())
        uwords = [i for i in uwords if i and not i[0].isdigit()]
        string = " ".join([i for i in uwords if len(i) > 1])
        return string.strip()

    def send(self, string):
        """Send 'string' and return the response of aspell"""
        string = self._clean(string)
        logdbg and logdbg("check: %r" % string)
        words = string.split()
        logdbg and logdbg("words: %r" % words)
        resp = []
        if not words:
            return words, resp
        self.inp.write(string + os.linesep)
        self.inp.flush()
        while True:
            try:
                resp.append(self.queue.get_nowait())
                # logdbg and logdbg('returns: %r' % resp[-1])
            except queue.Empty:
                if len(resp) == 0 or resp[-1] != "":
                    continue
                break
        if len(resp) != len(words) + 1:
            loginfo("warning: expected answer of aspell:\n words=%r\n resp=%r" % (words, resp))
        # logdbg and logdbg('resp: %r' % zip(words, resp))
        return words, resp

    def check(self, string):
        """Return the unknown words"""
        words, resp = self.send(string)
        unknown = [w for w, r in zip(words, resp) if r.strip() != "*"]
        return unknown


def get_cata_msg(catamess):
    """Import a messages file"""
    import importlib

    cata_msg = {}
    try:
        d = {}
        mod = __import__("code_aster.Messages.%s" % catamess, d, d, [catamess])
        importlib.reload(mod)
        cata_msg = getattr(mod, "cata_msg", {})
    except UnicodeDecodeError:
        dict_args = dict(
            valk=(
                "Encodage invalide pour le fichier de messages : '%s'" % catamess,
                traceback.format_exc(),
            )
        )
        UTMESS("F", "CATAMESS_1", **dict_args)
    except Exception:
        dict_args = dict(
            valk=("Nom du fichier de messages : '%s'" % catamess, traceback.format_exc())
        )
        UTMESS("F", "CATAMESS_1", **dict_args)
    return cata_msg


def check_format(checker, catamess, idmess, msg, typ):
    """Check the format used for the given type 'typ'."""
    dre = {
        "integer": (REI, ["d", "i"]),
        "real": (RER, ["f", "g", "e", "F", "G", "E"]),
        "string": (REK, ["s"]),
        "other1": (RE1, ["ktout"]),
        "other2": (RE2, []),
    }
    expr, l_auth = dre[typ]
    arg = expr.findall(msg)
    unauth = set(arg).difference(l_auth)
    if len(unauth) > 0:
        checker.error("%s : invalid format for type '%s' : %s" % (idmess, typ, tuple(unauth)))


def check_msg(checker, catamess, msg, key, lang):
    """Check a message."""
    idmess = "%s_%s" % (catamess, key)
    # check type : a string expected
    if type(msg) is not str:
        checker.error(lang, "%s has a wrong type" % idmess)
    if msg.strip() == "" and catamess != "vide":
        checker.error("%s is empty : use VIDE_1 for that!" % idmess)
    assert is_int(key), "unexpected key : %s" % key
    # check unauthorized formatting
    for re_unauth in RE_UNAUTH:
        mat = re_unauth.search(msg)
        if mat:
            checker.info("unrecommanded formatting in %s : %s" % (idmess, mat.groups()))
    # check formatting
    txt = None
    try:
        txt = msg % VALUES
    except Exception as exc:
        trace = repr(exc)
        checker.error("%s can not be formatted :\nmessage: %r\n%s" % (idmess, msg, trace))
    # check arguments
    for typ in ("integer", "real", "string", "other1", "other2"):
        check_format(checker, catamess, idmess, msg, typ)
    if txt and lang == "fr":
        logdbg and logdbg("idmess: %s" % idmess)
        unknown = checker.aspell.check(txt)
        if unknown:
            idmess = "%s_%s" % (catamess, key)
            unknown = [convert(word) for word in unknown]
            checker.warning_spell(idmess, unknown)


_pws = None


def get_personal_dict():
    """Build the dictionnary for Code_Aster : keywords + personal dict."""
    global _pws
    if _pws:
        return _pws
    cnt = ["personal_ws-1.1 fr 0 %s" % ENCODING]
    cata = "fort.34"
    cnt.extend(build_cata_dict(cata))

    dictdir = ExecutionParameter().get_option("rcdir")
    cawl = osp.join(dictdir, "code_aster_dict.aspell.per")
    if osp.exists(cawl):
        # ignore the first line
        with open(cawl, "r") as fper:
            cnt.extend(fper.read().splitlines()[1:])
    else:
        raise IOError(
            "no such file: {0}\nAn updated devtools repository is " "required!".format(cawl)
        )
    fd, _pws = tempfile.mkstemp(dir=os.getcwd())
    with open(fd, "w") as fobj:
        fobj.write(os.linesep.join([line for line in cnt if line.strip()]))
    return _pws


def build_cata_dict(filename):
    """Build the content of the dictionnary of the code_aster keywords.

    Arguments:
        repo (str): Path to 'src' repository.
        branches (list[str]): List of branches to be used.
        tip (bool): If *True*, uses the tip of each branch. *False* by default.

    Returns:
        str: Text of the catalog.
    """
    rkw = re.compile("([A-Z]+[A-Z_]+)", re.M | re.I)
    reg_ign = re.compile("(#.*$)", re.MULTILINE)
    allkw = set()
    with open(filename, "r") as fobj:
        txt = fobj.read()
    txt = reg_ign.sub("", txt)
    allkw.update(rkw.findall(txt))
    skw = set()
    for kw in allkw:
        skw.update(kw.split("_"))
    lkw = [kw for kw in skw if len(kw) > 1]
    lkw.sort()
    lkw.append("")
    return lkw


def check_catamess(checker, lang, l_cata):
    """Check all the messages files"""
    checker.set_current_lang(lang)
    checker.set_current_mod("-")
    loginfo("<i18n> lang=%s, domain=%s, localedir=%s" % (LO.current_lang, LO.domain, LO.localedir))
    tr = LO.translation(lang)
    if lang != "fr" and not isinstance(tr, gettext.GNUTranslations):
        checker.warning("no translation object for language '%s'" % lang)
        return
    pwl = get_personal_dict()
    if not pwl and lang == "fr":
        checker.error("Code_Aster personal dict not found: %s" % pwl)
    for catamess in l_cata:
        checker.set_current_mod(catamess)
        loginfo("<%s> checking %s..." % (lang, catamess))
        if catamess == "dvp":
            continue
        cata_msg = get_cata_msg(catamess)
        for key, msg in list(cata_msg.items()):
            if type(msg) is dict:
                msg = msg["message"]
            check_msg(checker, catamess, msg, key, lang)


def timekeeper(pid, delay):
    """Kill 'pid' if it times out"""
    time.sleep(delay)
    valk = (
        """
The process %d timed out after %d seconds.

It probably blocks reading the response of aspell on its stdout...
Try run with INFO=2 to have all the details.
"""
        % (pid, delay),
        """Interruption : kill the main process!""",
    )
    UTMESS("E", "CATAMESS_1", valk=valk)
    os.kill(pid, signal.SIGTERM)


def supv002_ops(self, ERREUR, **kwargs):
    """Fake macro-command to check messages"""
    global logdbg
    if kwargs.get("INFO") == 2:
        logdbg = loginfo
    # existing errors
    previous_errors = set(ERREUR)
    os.environ["LANG"] = "fr_FR.utf8"
    # remove all LC_xxxx variables
    keys = [k for k in list(os.environ.keys()) if k.startswith("LC_")]
    for k in keys:
        del os.environ[k]
    msgdir = osp.dirname(Messages.__file__)
    LCATA = [osp.basename(osp.splitext(cata)[0]) for cata in glob(osp.join(msgdir, "*.py"))]
    # LCATA = [osp.basename(osp.splitext(cata)[0]) for cata in glob(osp.join(msgdir, 'mecanonline9.py'))]

    # check for installation problem: http://bugs.python.org/issue3770
    do_check = True
    try:
        MPR.Queue()
        do_check = AspellCall.check_aspell()
    except ImportError as exc:
        if "sem_open implementation" in str(exc):
            do_check = False
            print("\n  <A> Problem detected ! supv002a can not run on this machine\n\n")
    if do_check:
        try:
            aspell = AspellCall(get_personal_dict(), "fr", ENCODING)
            checker = Checker(aspell)
            checker.set_ignored(ERREUR)
            # check default/native messages
            check_catamess(checker, "fr", LCATA)
            # check translated messages in english
            check_catamess(checker, "en", LCATA)
            # close aspell
            aspell.close()
        except OSError as exc:
            checker = Checker(None)
            checker.error("Can not start aspell: {0}".format(exc))

        nberr = checker.count_errors()
        errors = checker.get_errors() + """\nNumber of errors : %d""" % nberr
        allw = checker.get_all_words()
        nbwrn = len(allw)
        warns = checker.get_warnings()
        torm = list(previous_errors.difference(allw))
        torm.sort()
        new = checker.get_unknown_words()
        nbnew = len(new)
    else:
        nberr = nbnew = 0
        nbwrn = len(previous_errors)
        warns = ""
        errors = []
        torm = []

    if kwargs.get("unittest"):
        print(warns)
        print(errors)
        return 1

    if nberr > 0:
        UTMESS("A", "CATAMESS_1", valk=("%6d erreurs" % nberr, errors))
    if kwargs.get("INFO") == 2:
        warns = checker.get_all_warnings()
    if warns:
        UTMESS("A", "CATAMESS_1", valk=("Liste des alarmes et des erreurs par message", warns))
    if nbnew > 0:
        valk = ("Liste des nouvelles erreurs introduites à corriger :", str(new))
        UTMESS("A", "CATAMESS_1", valk=valk)
    if torm:
        UTMESS(
            "A",
            "CATAMESS_1",
            valk=(
                "Liste des erreurs qui n'apparaissent plus " "(à supprimer du mot-clé ERREUR) :",
                str(torm),
            ),
        )

    __tab = CREA_TABLE(LISTE=_F(PARA="NBERR", LISTE_I=nberr))

    TEST_TABLE(REFERENCE="ANALYTIQUE", VALE_CALC_I=0, VALE_REFE_I=0, NOM_PARA="NBERR", TABLE=__tab)

    __tnew = CREA_TABLE(LISTE=_F(PARA="NBNEW_ERR", LISTE_I=nbnew))

    TEST_TABLE(
        REFERENCE="ANALYTIQUE", VALE_CALC_I=0, VALE_REFE_I=0, NOM_PARA="NBNEW_ERR", TABLE=__tnew
    )

    __tabw = CREA_TABLE(LISTE=_F(PARA="NBWARN", LISTE_I=nbwrn))

    TEST_TABLE(
        CRITERE="ABSOLU",
        REFERENCE="ANALYTIQUE",
        VALE_CALC_I=len(previous_errors),
        VALE_REFE_I=0,
        PRECISION=len(previous_errors),
        NOM_PARA="NBWARN",
        TABLE=__tabw,
    )
    return


if __name__ != "__main__":
    from code_aster.Cata.Syntax import MACRO, SIMP
    from code_aster.Supervis.ExecuteCommand import UserMacro

    supv_cata = MACRO(
        nom="SUPV002",
        op=supv002_ops,
        ERREUR=SIMP(statut="o", typ="TXM", max="**"),
        INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
    )
    SUPV002 = UserMacro("SUPV002", supv_cata, supv002_ops)
else:
    # run as unittest
    # PYTHONPATH=$PYTHONPATH:/home/courtois/dev/codeaster/install/std/lib/python3.6/site-packages
    # ASTER_ROOT=/opt/aster
    # python -i astest/supv002a.33
    # logdbg = loginfo
    # aspell = AspellCall(get_personal_dict(), 'fr', 'utf-8')
    # unk = aspell.check("On ne peut pas vérifier les fotes d'orthografe.")
    # assert len(unk) == 2, unk
    # del aspell
    supv002_ops(None, [], unittest=True)
