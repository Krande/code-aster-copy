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

# aslint: disable=C4015

"""
   Description des OJB jeveux
"""
import pydoc
import sys
import traceback

import aster

from .ascheckers import CheckLog
from .asnom import SDNom
from .basetype import Type


class AsBase(Type):
    nomj = SDNom()
    optional = False

    def __init__(self, nomj=None, *args, **kwargs):
        super(AsBase, self).__init__(nomj, *args, **kwargs)
        assert self.nomj is not self.__class__.nomj
        if isinstance(nomj, str):
            self.nomj.nomj = nomj
        elif isinstance(nomj, SDNom):
            self.nomj.update(nomj.__getstate__())

    def set_name(self, nomj):
        """Positionne le nomj de self"""
        assert isinstance(self.nomj.nomj, str), "uniquement pour les concepts"
        self.nomj.nomj = nomj

    def check(self, checker=None):
        if checker is None:
            checker = CheckLog()

        # vérif déjà faite ? (en tenant compte du type)
        if checker.checkedAsBase(self):
            return checker
        checker.visitAsBase(self)

        # vérifie les enfants :
        optional = checker.optional
        checker.optional = checker.optional or self.optional
        for name in self._subtypes:
            v = getattr(self, name)
            if isinstance(v, (OJB, AsBase)):
                v.check(checker)
        for name in dir(self):
            if name.startswith("check_"):
                v = getattr(self, name)
                if callable(v):
                    try:
                        v(checker)
                    except:
                        mess = [
                            60 * "-",
                            "Erreur SDVERI (Attention : vérification " "incomplète)",
                            traceback.format_exc(),
                        ]
                        checker.err(self, "\n".join(mess))

        checker.optional = optional
        return checker

    def dump(self, indent=""):
        l = []
        checkers = []
        nomj = self.nomj()
        if self.optional:
            f = "(f)"
        else:
            f = "(o)"
        l.append(f + " " + nomj)
        for name in self._subtypes:
            obj = getattr(self, name)
            if isinstance(obj, (AsBase, OJB)):
                l.append(obj.dump(indent))
        for name in dir(self):
            if name.startswith("check_"):
                obj = getattr(self, name)
                if callable(obj) and name.startswith("check_"):
                    checkers.append(obj)

        indent = " " * len(nomj)
        for checker in checkers:
            doc = pydoc.text.document(checker)
            for line in doc.splitlines():
                l.append(indent + line)
        return "\n".join(l)

    def short_repr(self):
        return "<%s(%x,%r)>" % (self.__class__.__name__, id(self), self.nomj())

    def long_repr(self):
        if not hasattr(self, "accessible") or not self.accessible():
            # hors Aster ou en par_lot='oui'
            return self.short_repr()
        else:
            self.debugPrint()
            return ""

    def __repr__(self):
        # par défaut, on fait court !
        return self.short_repr()


# -----------------------------------------------------------------------------
class JeveuxAttr:

    """Un attribut jeveux"""

    def __init__(self, name):
        self.name = name

    def __get__(self, obj, klass):
        raise NotImplementedError

    def check(self, attrname, obj, log):
        checker = getattr(obj, "_" + attrname, None)
        if checker is None:
            return True
        val = self.__get__(obj, obj.__class__)
        if callable(checker):
            return checker(obj, attrname, val, log)
        elif val == checker:
            return True
        else:
            log.err(obj, "Attribut incorrect %s %r != %r" % (self.name, val, checker))
            return False


# -----------------------------------------------------------------------------


class JeveuxExists(JeveuxAttr):
    def __init__(self):
        pass

    def __get__(self, obj, klass):
        if obj is None:
            return self
        nomj = obj.nomj()
        if len(nomj) != 24:
            raise AssertionError(repr(nomj))
        return aster.jeveux_exists(nomj.ljust(24))


# -----------------------------------------------------------------------------


class JeveuxIntAttr(JeveuxAttr):
    def __get__(self, obj, klass):
        if obj is None:
            return self
        nomj = obj.nomj()
        if aster.jeveux_exists(nomj):
            return aster.jeveux_getattr(nomj, self.name)[0]
        else:
            return None


# -----------------------------------------------------------------------------


class JeveuxStrAttr(JeveuxAttr):
    def __get__(self, obj, klass):
        if obj is None:
            return self
        nomj = obj.nomj()
        if aster.jeveux_exists(nomj):
            return aster.jeveux_getattr(nomj, self.name)[1].strip()
        else:
            return None


# -----------------------------------------------------------------------------


class OJB(AsBase):
    _clas = None
    _genr = None
    _type = None
    _ltyp = None
    _xous = None
    _docu = None
    _exists = True

    clas = JeveuxStrAttr("CLAS")
    genr = JeveuxStrAttr("GENR")
    type = JeveuxStrAttr("TYPE")
    ltyp = JeveuxIntAttr("LTYP")
    xous = JeveuxStrAttr("XOUS")
    docu = JeveuxStrAttr("DOCU")
    exists = JeveuxExists()
    nomj = SDNom()

    def __init__(self, nomj=None, **attrs):
        super(OJB, self).__init__(nomj, **attrs)
        self.foreachattr(self.setattribute, attrs)
        self.optional = attrs.get("optional", False)

    def setattribute(self, name, prop, attrs):
        _name = "_" + name
        if name in attrs:
            setattr(self, _name, attrs[name])

    def get(self):
        nomj = self.nomj()
        if aster.jeveux_exists(nomj):
            obj_simple = aster.jeveux_getattr(nomj, "XOUS")[1].strip() == "S"
            if obj_simple:
                return aster.getvectjev(nomj)
            else:
                return aster.getcolljev(nomj)
        else:
            return None

    def get_stripped(self):
        """Fonction utilitaire, renvoie une liste de chaines 'strippées'"""
        data = self.get()
        if data is not None:
            return [x.strip() for x in data]
        else:
            return []

    def foreachattr(self, callback, *args, **kwargs):
        klass = self.__class__
        for k in dir(klass):
            v = getattr(klass, k)
            if isinstance(v, JeveuxAttr):
                callback(k, v, *args, **kwargs)

    def check(self, checker=None):
        if checker is None:
            checker = CheckLog()
        # l'objet a déjà été vérifié, on ne fait rien
        if checker.checkedOJB(self):
            return checker
        checker.visitOJB(self)
        if self.exists:
            self.foreachattr(lambda k, v, obj, c: v.check(k, obj, c), self, checker)
        else:
            if not self.optional and not checker.optional:
                checker.err(self, "n'existe pas (%r)" % self._parent)
        return checker

    def dump(self, indent=""):
        if self.optional:
            f = "(f)"
        else:
            f = "(o)"
        return f + " " + self.nomj() + " " + str(self.exists)


# -----------------------------------------------------------------------------


def Facultatif(ojb):
    ojb.optional = True
    return ojb


# -----------------------------------------------------------------------------


class OJBVect(OJB):
    lonmax = JeveuxIntAttr("LONMAX")
    lonuti = JeveuxIntAttr("LONUTI")
    _xous = "S"
    _genr = "V"


# -----------------------------------------------------------------------------


class OJBPtnom(OJB):
    nommax = JeveuxIntAttr("NOMMAX")
    nomuti = JeveuxIntAttr("NOMUTI")
    _xous = "S"
    _genr = "N"
    _type = "K"


# -----------------------------------------------------------------------------


class OJBCollec(OJB):
    stockage = JeveuxStrAttr("STOCKAGE")
    nutioc = JeveuxIntAttr("NUTIOC")
    acces = JeveuxStrAttr("ACCES")
    modelong = JeveuxStrAttr("MODELONG")
    nmaxoc = JeveuxIntAttr("NMAXOC")


# -----------------------------------------------------------------------------


class AsVI(OJBVect):
    _type = "I"


# -----------------------------------------------------------------------------


class AsVS(OJBVect):
    _type = "S"


# -----------------------------------------------------------------------------


class AsVR(OJBVect):
    _type = "R"


# -----------------------------------------------------------------------------


class AsVC(OJBVect):
    _type = "C"


# -----------------------------------------------------------------------------


class AsVL(OJBVect):
    _type = "L"


# -----------------------------------------------------------------------------


class AsVK8(OJBVect):
    _type = "K"
    _ltyp = 8


# -----------------------------------------------------------------------------


class AsVK16(OJBVect):
    _type = "K"
    _ltyp = 16


# -----------------------------------------------------------------------------


class AsVK24(OJBVect):
    _type = "K"
    _ltyp = 24


# -----------------------------------------------------------------------------


class AsVK32(OJBVect):
    _type = "K"
    _ltyp = 32


# -----------------------------------------------------------------------------


class AsVK80(OJBVect):
    _type = "K"
    _ltyp = 80


# Pour compatibilite
AsObject = OJB
AsColl = OJBCollec
AsPn = OJBPtnom
AsVect = OJBVect
