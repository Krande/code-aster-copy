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

import hashlib
import os
import re
import sys
import tempfile
from optparse import OptionGroup, OptionParser

import aster
from ..Cata.Syntax import _F
from ..CodeCommands import CREA_TABLE, INFO_EXEC_ASTER, TEST_TABLE
from ..Messages import UTMESS


class TestFichierError(Exception):
    pass


def convert(x):
    return float(x)


def f_SOMM(somme, lx):
    return somme + sum([convert(x) for x in lx])


def f_SOMM_ABS(somme, lx):
    return somme + sum([abs(convert(x)) for x in lx])


def f_MINI(val, lx):
    return min(val, min([convert(x) for x in lx]))


def f_MAXI(val, lx):
    return max(val, max([convert(x) for x in lx]))


def f_MINI_ABS(val, lx):
    return min(val, min([abs(convert(x)) for x in lx]))


def f_MAXI_ABS(val, lx):
    return max(val, max([abs(convert(x)) for x in lx]))


dict_func_test = {
    "SOMM": f_SOMM,
    "SOMM_ABS": f_SOMM_ABS,
    "MINI": f_MINI,
    "MAXI": f_MAXI,
    "MINI_ABS": f_MINI_ABS,
    "MAXI_ABS": f_MAXI_ABS,
}

# -------------------------------------------------------------------------


def test_fichier_ops(self, **kwargs):
    """
    Macro permettant de tester la non-regression d'un fichier.
    On teste le nombre de réels présents, et, facultativement, la
    somme de ces nombres et le texte du fichier.
    """

    FICHIER = kwargs.get("FICHIER")
    NB_VALE = kwargs.get("NB_VALE")
    NB_VALE_I = kwargs.get("NB_VALE_I")
    VALE_CALC = kwargs.get("VALE_CALC")
    VALE_CALC_I = kwargs.get("VALE_CALC_I")
    VALE_CALC_K = kwargs.get("VALE_CALC_K")
    TOLE_MACHINE = kwargs.get("TOLE_MACHINE")
    CRITERE = kwargs.get("CRITERE")
    INFO = kwargs.get("INFO")

    # On importe les definitions des commandes a utiliser dans la macro
    # Le nom de la variable doit etre obligatoirement le nom de la commande
    #
    is_ok = 0
    TYPE_TEST = kwargs.get("TYPE_TEST")
    TYPE_TEST = TYPE_TEST or "SOMM"
    # # vérification non faisable dans le catalogue
    if VALE_CALC_I is not None and NB_VALE_I is None:
        UTMESS("F", "TEST0_5")
    # vérifier que le fichier a été fermé
    __tinfo = INFO_EXEC_ASTER(LISTE_INFO="ETAT_UNITE", FICHIER=FICHIER)
    if __tinfo["ETAT_UNITE", 1].find("OUVERT") > -1:
        UTMESS("S", "TEST0_2", valk=FICHIER)
    # lecture du fichier
    if not os.path.isfile(FICHIER):
        UTMESS("S", "TEST0_3", valk=FICHIER)

    with open(FICHIER, "r") as fileobj:
        # filtre par expression régulière
        try:
            fileobj = regexp_filter(fileobj, kwargs.get("EXPR_IGNORE"))
        except TestFichierError as valk:
            UTMESS("S", "TEST0_1", valk=valk)
        # calcule le nombre de valeurs et la somme ou min/max
        verbose = INFO > 1
        results = test_iter(fileobj, function=dict_func_test[TYPE_TEST], verbose=verbose)
        nbvalr, vale_r, nbvali, vale_i, chksum = results

    # produit le TEST_TABLE
    refsum = VALE_CALC_K or "not_tested"
    is_ok = int(chksum == refsum)
    mcfact = [
        _F(PARA="NBVAL_I", LISTE_I=nbvali),
        _F(PARA="VALE_I", LISTE_I=vale_i),
        _F(PARA="NBVAL", LISTE_I=nbvalr),
        _F(PARA="VALE", LISTE_R=vale_r),
        _F(PARA="TEXTE", LISTE_I=is_ok),
    ]
    __tab1 = CREA_TABLE(LISTE=mcfact)
    # message
    UTMESS("I", "TEST0_4", valk=FICHIER)
    UTMESS("I", "TEST0_13")
    if verbose or NB_VALE_I is not None:
        UTMESS("I", "TEST0_15", vali=(nbvali, NB_VALE_I or 0))
    if verbose or VALE_CALC_I is not None:
        UTMESS("I", "TEST0_16", vali=(vale_i, VALE_CALC_I or 0))
    UTMESS("I", "TEST0_14")
    UTMESS("I", "TEST0_15", vali=(nbvalr, NB_VALE or 0))
    if verbose or VALE_CALC is not None:
        UTMESS("I", "TEST0_17", valr=(vale_r, VALE_CALC or 0.0))
    if verbose or VALE_CALC_K is not None:
        UTMESS("I", "TEST0_18", valk=(chksum, refsum))
    # tests
    TEST_TABLE(
        TABLE=__tab1, NOM_PARA="NBVAL", VALE_CALC_I=NB_VALE, CRITERE="ABSOLU", TOLE_MACHINE=0
    )
    if VALE_CALC is not None:
        TEST_TABLE(
            TABLE=__tab1,
            NOM_PARA="VALE",
            CRITERE=CRITERE,
            VALE_CALC=VALE_CALC,
            TOLE_MACHINE=TOLE_MACHINE,
        )
    if NB_VALE_I is not None:
        TEST_TABLE(
            TABLE=__tab1,
            NOM_PARA="NBVAL_I",
            VALE_CALC_I=NB_VALE_I,
            CRITERE="ABSOLU",
            TOLE_MACHINE=0,
        )
    if VALE_CALC_I is not None:
        TEST_TABLE(
            TABLE=__tab1,
            NOM_PARA="VALE_I",
            CRITERE=CRITERE,
            VALE_CALC_I=VALE_CALC_I,
            TOLE_MACHINE=TOLE_MACHINE,
        )
    if VALE_CALC_K is not None:
        TEST_TABLE(
            TABLE=__tab1, NOM_PARA="TEXTE", VALE_CALC_I=int(True), TOLE_MACHINE=0, CRITERE="ABSOLU"
        )
    return


def regexp_filter(file_in, regexp_ignore, debug=False):
    """Filtre le fichier fourni (file descriptor) en utilisant les
    expressions régulières fournies.
    On retourne l'objet file vers le fichier modifié (ou non).
    """
    if not regexp_ignore:  # None or []
        return file_in
    # vérification des expressions régulières
    if type(regexp_ignore) not in (list, tuple):
        regexp_ignore = [regexp_ignore]
    l_regexp = []
    for exp in regexp_ignore:
        try:
            obj = re.compile(exp)
        except re.error as s:
            raise TestFichierError(s, str(exp))
        else:
            l_regexp.append(obj)
    # filtre du fichier
    file_out = tempfile.NamedTemporaryFile()
    file_in.seek(0)
    for i, line in enumerate(file_in):
        if debug:
            print("LIGNE", i, end=" ")
        keep = True
        for exp in l_regexp:
            if exp.search(line):
                keep = False
                if debug:
                    print(" >>>>>>>>>> IGNOREE <<<<<<<<<<")
                break
        if keep:
            file_out.write(line.encode())
            if debug:
                print()
    file_out.seek(0)
    return file_out


RE_FLOAT_EXPO = re.compile(r"[-+]?[0-9]{1}[0-9\.]+[eED][\-\+]{0,1}[0-9]+")
RE_FLOAT = re.compile(r"[-+]?[0-9]+?\.[0-9]*")
RE_INT = re.compile("[0-9]+")


def test_iter(obj, function, verbose=False):
    """
    Cette fonction compte le nombre de réels dans le fichier et une grandeur
    à partir des valeurs (somme, sommes des valeurs absolues, min/max...).
    IN :
       obj      : objet 'file' ou 'string' sur le lequel on peut itérer
       function : fonction de test   val = func_test(val, [xi, ...])
       verbose  : on affiche le résumé si info>0
    OUT :
       nombre de valeurs réelles, résultat de la fonction sur les réels,
       nombre de valeurs entières, résultat de la fonction sur les entiers,
       somme de contrôle du texte restant.
    Le résultat entier est systématiquement retourné modulo 2147483647
    pour être homogène avec les plate-formes 32 bits.
    Le nombre de valeurs lui est supposé être inférieur à cette valeur.
    """
    max_buff_size = 1000
    nbvalr = 0
    nbvali = 0
    valr = 0.0
    vali = 0.0
    hfile = hashlib.md5()
    # Si on lit tout le fichier d'un coup, on va environ 3 fois plus vite
    # que si on le lit ligne à ligne, mais on consomme en mémoire environ
    # 5 fois la taille du fichier...
    # En lisant par paquet de 1000 (ou 10000), on va quasiment aussi vite
    # en consommant très peu de mémoire.
    #    fichier     tout   ligne/ligne   1000 lignes
    #     10 Mo       3 s      10 s       3 s
    #     50 Mo      17 s      48 s      17 s
    #    100 Mo      34 s      96 s      35 s
    # l'itérateur est l'objet file lui-même ou on le crée sur la liste
    obj.seek(0)
    iterator = iter(obj)
    ok = True
    buff = []
    while ok:
        try:
            text = next(iterator)
            if type(text) is bytes:
                text = text.decode()
        except StopIteration:
            ok = False
            text = ""
        buff.append(text)
        if ok and len(buff) < max_buff_size:
            continue
        else:
            text = "".join(buff)
            buff = []
        # extract floats
        l_float = RE_FLOAT_EXPO.findall(text)
        l_float = [s.replace("D", "E") for s in l_float]
        text = RE_FLOAT_EXPO.sub("", text)
        l_float.extend(RE_FLOAT.findall(text))
        text = RE_FLOAT.sub("", text)
        nbvalr += len(l_float)
        valr = function(valr, l_float)
        # extract integers
        l_int = RE_INT.findall(text)
        text = RE_INT.sub("", text)
        nbvali += len(l_int)
        vali = function(vali, l_int)
        # add text
        text = "".join([s.strip() for s in text.split()])
        hfile.update(text.encode())
        if verbose:
            print("Nombres réels :", nbvalr)
            print(l_float)
            print("Nombres entiers :", nbvali)
            print(l_int)
            print("Texte :")
            print(text)
    chksum = hfile.hexdigest()
    return nbvalr, valr, nbvali, int(vali) % 2147483647, chksum


if __name__ == "__main__":
    p = OptionParser(usage="usage: %s fichier [options]" % sys.argv[0])
    p.add_option(
        "--type_test",
        action="store",
        dest="type_test",
        default="SOMM",
        help="type du test : SOMM, SOMM_ABS, MIN, MAX",
    )
    p.add_option(
        "--expr_ignore",
        action="store",
        dest="exp",
        type="string",
        help="expression régulière à ignorer",
    )
    p.add_option(
        "-v", "--verbose", action="store_true", dest="verbose", default=False, help="mode bavard"
    )
    opts, args = p.parse_args()

    if len(args) == 0:
        p.error("fichier à tester ?")

    if opts.exp is None:
        exp = []
    else:
        exp = [opts.exp]

    fileobj = open(args[0], "r")
    fileobj = regexp_filter(fileobj, exp)
    results = test_iter(fileobj, function=dict_func_test[opts.type_test], verbose=opts.verbose)
    print("%6d réels, vale_r = %f, %6d entiers, vale_i = %d, texte : %s" % results)
