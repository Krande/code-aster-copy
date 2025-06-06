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

# TODO solve dependency TablePy/Graph
# aslint: disable=C4008

import os
import shutil
import re
import sys
from copy import copy

import numpy

from ..Messages import UTMESS
from ..Utilities import cut_long_lines, is_float, is_number, is_sequence, is_str, transpose
from ..Utilities.misc import fmtF2PY

# formats de base (identiques à ceux du module Graph)
DicForm = {
    "chead": None,  # entête avant la table
    "cfoot": "",  # lignes après la table
    "csep": " ",  # séparateur
    "ccom": "#",  # commentaire
    "ccpara": "",  # commentaire des labels
    "cdeb": "",  # début de ligne
    "cfin": "\n",  # fin de ligne
    "sepch": ";",  # séparateur entre deux lignes d'une cellule
    "formK": "%-12s",  # chaines
    "formR": "%12.5E",  # réels
    "formI": "%12d",  # entiers
}
# type par défaut des chaines de caractères
Kdef = "K24"


class TableBase:

    """Classe pour partager les méthodes d'impression entre Table et Colonne
    (c'est surtout utile pour vérifier que l'extraction et les filtres sur les
    colonnes sont corrects).
    """

    def __init__(self):
        """Constructeur."""
        self.rows = None
        self.para = None
        self.type = None
        self.titr = None

    def __repr__(self):
        return self.ReprTable()

    def Croise(self, **kargs):
        raise NotImplementedError("Must be defined in a derived class")

    def __len__(self):
        """Retourne le nombre de ligne dans la Table/Colonne."""
        return len(self.rows)

    def Impr(self, FICHIER=None, FORMAT="TABLEAU", dform=None, **opts):
        """Impresssion de la Table selon le format spécifié.
        FICHIER : nom du(des) fichier(s). Si None, on dirige vers stdout
        dform : dictionnaire de formats d'impression (format des réels,
           commentaires, saut de ligne...)
        opts  : selon FORMAT.
        """
        para = {
            "TABLEAU": {"mode": "a", "driver": self.ImprTableau},
            "ASTER": {"mode": "a", "driver": self.ImprTableau},
            "NUMPY": {"mode": "a", "driver": self.ImprNumpy},
            "XMGRACE": {"mode": "a", "driver": self.ImprGraph},
            "AGRAF": {"mode": "a", "driver": self.ImprTableau},
            "TABLEAU_CROISE": {"mode": "a", "driver": self.ImprTabCroise},
        }
        kargs = {
            "FICHIER": FICHIER,
            "FORMAT": FORMAT,
            "dform": DicForm.copy(),
            "mode": para[FORMAT]["mode"],
        }
        if dform is not None and type(dform) is dict:
            kargs["dform"].update(dform)
        # ajout des options
        kargs.update(opts)

        if not kargs.get("PAGINATION"):
            # call the associated driver
            para[FORMAT]["driver"](**kargs)

        else:
            if not is_sequence(kargs["PAGINATION"]):
                ppag = [kargs["PAGINATION"]]
            else:
                ppag = list(kargs["PAGINATION"])
            del kargs["PAGINATION"]
            npag = len(ppag)
            # paramètres hors ceux de la pagination
            lkeep = [p for p in self.para if ppag.count(p) == 0]
            # création des listes des valeurs distinctes
            lvd = []
            for p in ppag:
                lvp = list(getattr(self, p).values())
                lvn = []
                for it in lvp:
                    if it is not None and lvn.count(it) == 0:
                        lvn.append(it)
                lvn.sort()
                lvd.append(lvn)
            # création des n-uplets
            s = "[[" + ",".join(["x" + str(i) for i in range(npag)]) + "] "
            s += " ".join(["for x" + str(i) + " in lvd[" + str(i) + "]" for i in range(npag)]) + "]"
            try:
                lnup = eval(s)
            except SyntaxError as s:
                UTMESS("F", "TABLE0_20")
            # pour chaque n-uplet, on imprime la sous-table
            for nup in lnup:
                tab = self
                for i in range(npag):
                    tab = tab & (getattr(tab, ppag[i]) == nup[i])
                    sl = ""
                    if tab.titr:
                        sl = "\n"
                    tab.titr += sl + ppag[i] + ": " + str(nup[i])
                tab[lkeep].Impr(**kargs)

    def ImprTableau(self, **kargs):
        """Impression au format TABLEAU ou ASTER"""
        # fichier ou stdout
        if kargs.get("FICHIER") is not None:
            f = open(kargs["FICHIER"], kargs["mode"])
        else:
            f = sys.stdout
        # ecriture
        f.write(self.ReprTable(**kargs) + "\n")
        f.flush()
        # fermeture
        if kargs.get("FICHIER") is not None:
            f.close()

    def ImprNumpy(self, **kargs):
        """Impression au format numpy."""
        para = [i for i, j in zip(self.para, self.type) if j not in ("R", "I")]
        if para:
            UTMESS("F", "TABLE0_5", valk=repr(para))
        values = []
        for i in self.para:
            valp = list(getattr(self, i).values())
            # replace None by 0.
            valp = [x or 0.0 for x in valp]
            values.extend(valp)

        tab = numpy.array(values)
        tab.shape = (len(self.para), len(self))
        tab = tab.T
        filename = kargs.get("FICHIER")
        if filename is None:
            print(tab)
        else:
            numpy.save(filename, tab)
            if not filename.endswith(".npy"):
                shutil.move(filename + ".npy", filename)

    def ReprTable(self, FORMAT="TABLEAU", dform=None, **ignore):
        """Représentation d'une Table ou d'une Colonne sous forme d'un tableau."""
        rows = self.rows
        para = self.para
        typ = self.type
        if not is_sequence(para):
            para = [self.para]
            typ = [self.type]
        if dform is None:
            dform = DicForm.copy()
        if dform["chead"] is None:
            dform["chead"] = os.linesep.join(
                [dform["ccom"], dform["ccom"] + "-" * 80, dform["ccom"]]
            )
        # est-ce que l'attribut .type est renseigné ?
        typdef = typ != [None] * len(typ)
        ASTER = FORMAT == "ASTER"
        lspa = []
        if ASTER:
            # ['']+ pour ajouter un séparateur en début de ligne
            lspa.append("")
        # lmax : largeur max des colonnes = max(form{K,R,I},len(parametre))
        lmax = []
        # formats
        strfmt, strfmt_none = {}, {}
        for t, p in zip(typ, para):
            larg_max = max(
                [len(str(p))] + [len(FMT(dform, k, t) % 0) for k in ("formK", "formR", "formI")]
            )
            lspa.append(FMT(dform, "formK", t, larg_max, str(p)) % p)
            lmax.append(larg_max)
            assert t is not None, "Type de la colonne '%s' non défini" % p
            strfmt[p] = FMT(dform, "form" + t[0], t, larg_max)
            strfmt_none[p] = FMT(dform, "formK", t, larg_max)
        if typdef:
            stype = dform["csep"].join(
                [""] + [FMT(dform, "formK", typ[i], lmax[i]) % typ[i] for i in range(len(para))]
            )
        txt = []
        if ASTER:
            txt.append("#DEBUT_TABLE")
        if self.titr:
            if ASTER:
                txt.extend(["#TITRE " + lig for lig in self.titr.split("\n")])
            else:
                txt.extend([dform["ccom"] + lig for lig in self.titr.split("\n")])
        if " " not in dform["csep"]:
            lspa = [i.strip() for i in lspa]
        txt.append(dform["ccpara"] + dform["csep"].join(lspa))
        if ASTER and typdef:
            txt.append(stype)
        for r in rows:
            lig = []
            if ASTER:
                lig.append("")
            empty = True
            for t, p, lmax_i in zip(typ, para, lmax):
                val = r.get(p)
                if val is not None:
                    empty = False
                    lig.append(strfmt[p] % val)
                else:
                    s = strfmt_none[p] % "-"
                    # format AGRAF = TABLEAU + '\' devant les chaines de
                    # caractères !
                    if FORMAT == "AGRAF":
                        s = "\\" + s
                    lig.append(s)
            if not empty:
                lig2 = [dform["sepch"].join(ch.splitlines()) for ch in lig]
                txt.append(dform["csep"].join(lig2))
        if ASTER:
            txt.append("#FIN_TABLE")
        # ajout des debut et fin de ligne
        txt = [dform["cdeb"] + t + dform["cfin"] for t in txt]
        txt.insert(0, dform["chead"])
        txt.append(dform["cfoot"])

        return "".join(txt)

    def ImprTabCroise(self, **kargs):
        """Impression au format TABLEAU_CROISE d'une table ayant 3 paramètres."""
        # création du tableau croisé et impression au format TABLEAU
        tabc = self.Croise()
        kargs["FORMAT"] = "TABLEAU"
        tabc.Impr(**kargs)

    def ImprGraph(self, **kargs):
        """Impression au format XMGRACE : via le module Graph"""
        from .table_graph import AjoutParaCourbe, Graph

        args = kargs.copy()
        if len(self.para) != 2:
            UTMESS("A", "TABLE0_21")
            return
        # suppression des lignes contenant une cellule vide
        tnv = getattr(self, self.para[0]).NON_VIDE() & getattr(self, self.para[1]).NON_VIDE()
        # objet Graph
        graph = Graph()
        dicC = {
            "Val": [
                list(getattr(tnv, tnv.para[0]).values()),
                list(getattr(tnv, tnv.para[1]).values()),
            ],
            "Lab": tnv.para,
        }
        AjoutParaCourbe(dicC, args)
        graph.AjoutCourbe(**dicC)

        # Surcharge des propriétés du graphique et des axes
        # (bloc quasiment identique dans impr_fonction_ops)
        if args.get("TITRE"):
            graph.Titre = args["TITRE"]
        if args.get("BORNE_X"):
            graph.Min_X = args["BORNE_X"][0]
            graph.Max_X = args["BORNE_X"][1]
        if args.get("BORNE_Y"):
            graph.Min_Y = args["BORNE_Y"][0]
            graph.Max_Y = args["BORNE_Y"][1]
        if args.get("LEGENDE_X"):
            graph.Legende_X = args["LEGENDE_X"]
        if args.get("LEGENDE_Y"):
            graph.Legende_Y = args["LEGENDE_Y"]
        if args.get("ECHELLE_X"):
            graph.Echelle_X = args["ECHELLE_X"]
        if args.get("ECHELLE_Y"):
            graph.Echelle_Y = args["ECHELLE_Y"]
        if args.get("GRILLE_X"):
            graph.Grille_X = args["GRILLE_X"]
        if args.get("GRILLE_Y"):
            graph.Grille_Y = args["GRILLE_Y"]

        try:
            graph.Trace(**args)
        except TypeError:
            UTMESS("A", "TABLE0_22")


class Table(TableBase):

    """Une table est construite comme une liste de lignes, chaque ligne est
    un dictionnaire.
    On crée puis on ajoute les lignes avec la méthode append :
       t=Table()
       t.append(dict(a=1,b=2))
       t.append(dict(a=3,b=4))
    La méthode __iter__ définit un itérateur sur les lignes de la table,
    __repr__ retourne une représentation de la table, utilisée par "print t".
    Grace à la classe Colonne et à sa méthode _extract, il est possible
    de construire une sous-table qui satisfait un critère donné.
    Le critère est donné par une fonction Python qui retourne vrai
    ou faux si la valeur d'une colonne respecte le critère ou non.
    Exemple:
      def critere(valeur):
          return valeur < 10
      soustable = t.a._extract(critere)
    t.a retourne un objet intermédiaire de la classe Colonne qui mémorise
    le nom de la colonne demandée (a, ici).
    """

    def __init__(self, rows=[], para=[], typ=[], titr="", nom=""):
        """Constructeur de la Table :
        rows : liste des lignes (dict)
        para : liste des paramètres
        type : liste des types des paramètres
        titr : titre de la table
        nom : nom du concept table_sdaster dont est issue la table
        """
        self.rows = [r for r in rows if list(r.values()) != [None] * len(list(r.values()))]
        self.para = list(para)
        for i in self.para:
            if self.para.count(i) != 1:
                UTMESS("F", "TABLE0_23", valk=i)
        if len(typ) == len(self.para):
            self.type = list(typ)
        else:
            self.type = [None] * len(self.para)
        self.titr = titr
        self.nom = nom
        self.referenceToDataStructure = []

    def __getstate__(self):
        """Return internal state.

        Returns:
            tuple: Internal state.
        """
        return (self.rows, self.para, self.type, self.titr, self.nom)

    def __setstate__(self, state):
        """Restore internal state.

        Arguments:
            state (tuple): Internal state.
        """
        if isinstance(state, tuple) and len(state) == 5:
            self.rows, self.para, self.type, self.titr, self.nom = state

    def getName(self):
        """Return name of the Table"""
        return self.nom

    def copy(self):
        """Retourne une copie de la table."""
        rows = []
        for r in self.rows:
            rows.append(copy(r))
        return Table(rows, self.para[:], self.type[:], self.titr, self.nom)

    def add_para(self, para, typ):
        """Ajoute un nouveau paramètre."""
        if not is_sequence(para):
            para = [para]
        if not is_sequence(typ):
            typ = [typ]
        if len(typ) != len(para):
            typ = [typ[0]] * len(para)
        for p, t in zip(para, typ):
            if not p in self.para:
                self.para.append(p)
                self.type.append(t)

    def append(self, obj):
        """Ajoute une ligne (type dict) qui peut éventuellement définir un
        nouveau paramètre."""
        para = list(obj.keys())
        for p in para:
            if not p in self.para:
                self.add_para(p, typaster(obj[p]))
            else:
                ip = self.para.index(p)
                self.type[ip] = typaster(obj[p], self.type[ip])
        self.rows.append(obj)

    def extend(self, objs):
        """Ajoute plusieurs lignes (list of dict)."""
        for row in objs:
            self.append(row)

    def SansColonneVide(self, l_para=None):
        """Retourne une copie de la table dans laquelle on a supprimé les colonnes
        vides (les lignes vides sont automatiquement supprimées).
        """
        # ptest : colonnes potentiellement vides
        pkeep = l_para or self.para
        ptest = pkeep[:]
        for row in self:
            notNone = [p for p in ptest if row.get(p) is not None]
            ptest = [p for p in ptest if not p in notNone]
            if len(ptest) == 0:
                break
        # pkeep : on conserve les colonnes non vides
        pkeep = [p for p in pkeep if not p in ptest]
        return self[pkeep]

    def __setitem__(self, k_para, k_value):
        """Ajoute une colonne k_para dont les valeurs sont dans k_value"""
        if len(k_value) == 0:
            return
        if k_para in self.para:
            UTMESS("F", "TABLE0_24", valk=k_para)
        self.add_para(k_para, typ=typaster(k_value[0]))
        i = 0
        for row in self:
            if i < len(k_value):
                row[k_para] = k_value[i]
                self.type[-1] = typaster(k_value[i], self.type[-1])
            else:
                row[k_para] = None
            i += 1
        for j in range(i, len(k_value)):
            self.append({k_para: k_value[j]})

    def fromfunction(self, nom_para, funct, l_para=None, const=None):
        """Ajoute une colonne `nom_para` en évaluant la fonction `funct` sur
        la valeur des paramètres `l_para` (qui doivent exister dans la table).
        Si `l_para` n'est pas fourni, on prend `funct`.getVariables() (FORMULE Aster).
        On peut passer un dictionnaire de constantes dans `const`. Quand on
        utilise une FORMULE Aster, les constantes sont prises dans le contexte
        global.
        """
        # vérif préalables
        if not hasattr(funct, "__call__"):
            UTMESS("F", "TABLE0_25", valk=(funct.getName(), "__call__"))
        if nom_para in self.para:
            UTMESS("F", "TABLE0_24", valk=nom_para)
        if l_para is None:
            if not hasattr(funct, "getVariables"):
                UTMESS("F", "TABLE0_25", valk=(funct.getName(), "getVariables"))
            l_para = funct.getVariables()
        if not is_sequence(l_para):
            l_para = [l_para]
        not_found = ", ".join([p for p in l_para if not p in self.para])
        if not_found != "":
            UTMESS("F", "TABLE0_27", valk=not_found)
        if const is None:
            const = {}
        if type(const) is not dict:
            UTMESS("F", "TABLE0_28", valk=("const", "dict"))
        # liste des valeurs des paramètres
        tabpar = []
        for para in l_para:
            vals = list(getattr(self, para).values())
            tabpar.append(vals)
        tabpar = transpose(tabpar)
        # évaluation de la fonction sur ces paramètres
        vectval = []
        for lpar in tabpar:
            # si un paramètre est absent, on ne peut pas évaluer la formule
            if None in lpar:
                vectval.append(None)
            else:
                vectval.append(funct(*lpar, **const))
        # ajout de la colonne
        self[nom_para] = vectval

    def __iter__(self):
        """Itère sur les lignes de la Table"""
        return iter(self.rows)

    def __getattr__(self, column):
        """Construit un objet intermediaire (couple table, colonne)"""
        typ = None
        if not column in self.para:
            column = ""
        else:
            typ = self.type[self.para.index(column)]
        return Colonne(self, column, typ)

    def sort(self, CLES=None, ORDRE="CROISSANT"):
        """Tri de la table.
        CLES  : liste des clés de tri
        ORDRE : CROISSANT ou DECROISSANT
        """
        # par défaut, on prend tous les paramètres
        if CLES is None:
            CLES = self.para[:]
        # vérification des arguments
        if not is_sequence(CLES):
            CLES = [CLES]
        else:
            CLES = list(CLES)
        not_found = ", ".join([p for p in CLES if not p in self.para])
        if not_found != "":
            UTMESS("F", "TABLE0_27", valk=not_found)
        if not ORDRE in ("CROISSANT", "DECROISSANT"):
            UTMESS("F", "TABLE0_29", valk=ORDRE)
        # tri
        self.rows = sort_table(self.rows, self.para, self.type, CLES, (ORDRE == "DECROISSANT"))

    def __delitem__(self, args):
        """Supprime les colonnes correspondantes aux éléments de args"""
        if not is_sequence(args):
            args = [args]
        for item in args:
            try:
                ind_item = self.para.index(item)
            except ValueError:
                UTMESS("F", "TABLE0_27", valk=item)
            del self.type[ind_item]
            self.para.remove(item)
            for line in self.rows:
                if item in line:
                    del line[item]

    def __getitem__(self, args):
        """Extrait la sous table composée des colonnes dont les paramètres sont dans args"""
        if not is_sequence(args):
            args = [args]
        else:
            args = list(args)
        new_rows = []
        new_para = args
        new_type = []
        for item in new_para:
            if not item in self.para:
                return Table()
            new_type.append(self.type[self.para.index(item)])
        for line in self:
            new_line = {}
            for item in new_para:
                v = line.get(item)
                if v is not None:
                    new_line[item] = v
            new_rows.append(new_line)
        return Table(new_rows, new_para, new_type, self.titr, self.nom)

    def OrdreColonne(self, cols):
        """Réordonne les colonnes en mettant en premier 'cols'.
        Ignore les colonnes qui ne seraient pas dans 'self.para'."""
        if type(cols) not in (list, tuple):
            cols = [cols]
        new_para = [p for p in cols if p in self.para]
        others = [p for p in self.para if not p in cols]
        new_para.extend(others)
        new_type = []
        for item in new_para:
            new_type.append(self.type[self.para.index(item)])
        self.para = new_para
        self.type = new_type

    def _tuplevalues(self, dico):
        """Retourne la liste des valeurs d'une ligne dans l'ordre self.para
        ("hashable" pour en faire une clé de dict.)
        """
        return tuple(map(dico.get, self.para))

    def __and__(self, other):
        """Intersection de deux tables (opérateur &)"""
        if other.para != self.para:
            UTMESS("A", "TABLE0_30")
            return Table()
        else:
            dval_other = dict.fromkeys([self._tuplevalues(r) for r in other], 1)
            tmp = [r for r in self if dval_other.get(self._tuplevalues(r)) is not None]
            return Table(tmp, self.para, self.type, self.titr)

    def __or__(self, other):
        """Union de deux tables (opérateur |)"""
        if other.para != self.para:
            UTMESS("A", "TABLE0_30")
            return Table()
        else:
            tmp = self.rows[:]
            dval_self = dict.fromkeys([self._tuplevalues(r) for r in self], 1)
            tmp.extend([r for r in other if dval_self.get(self._tuplevalues(r)) is None])
            return Table(tmp, self.para, self.type[:], self.titr)

    def difference(self, other):
        """Différence de deux tables: retire les lignes qui sont dans 'other'"""
        if other.para != self.para:
            UTMESS("A", "TABLE0_30")
            return Table()
        else:
            dval_other = dict.fromkeys([self._tuplevalues(r) for r in other], 1)
            tmp = [r for r in self if dval_other.get(self._tuplevalues(r)) is None]
            return Table(tmp, self.para, self.type[:], self.titr)

    def values(self):
        """Renvoie la table sous la forme d'un dictionnaire de listes dont les
        clés sont les paramètres.
        """
        dico = {}
        for column in self.para:
            dico[column] = list(Colonne(self, column).values())
        return dico

    def dict_CREA_TABLE(self):
        """Renvoie le dictionnaire des mots-clés à fournir à la commande CREA_TABLE
        pour produire une table_sdaster.
        """
        self.titr = cut_long_lines(self.titr, 80)
        # il y a eu limite à 50 titres dans le fortran autant le limiter
        # maintenant
        dico = {"TITRE": ["%-80s" % lig for lig in self.titr.split("\n")][:1], "LISTE": []}
        # remplissage de chaque occurrence (pour chaque paramètre) du mot-clé
        # facteur LISTE
        for i in range(len(self.para)):
            # nom du paramètre et type si K*
            d = {"PARA": self.para[i]}
            typ = self.type[i]
            if typ[0] == "K":
                mc = "LISTE_K"
                if not typ in ("K8", "K16", "K24"):
                    UTMESS("A", "TABLE0_32", valk=(self.para[i], Kdef))
                    typ = Kdef
                d["TYPE_K"] = typ
            elif typ == "I":
                mc = "LISTE_I"
            elif typ == "R":
                mc = "LISTE_R"
            else:
                UTMESS("F", "TABLE0_31", valk=self.para[i])
            # valeurs sans trou / avec trou
            vals = list(getattr(self, self.para[i]).values())
            if typ == "R":
                try:
                    check_nan(vals)
                    check_inf(vals)
                except ValueError as err:
                    UTMESS("F", "TABLE0_33", valk=(self.para[i], str(err)))
            if vals.count(None) == 0:
                d[mc] = vals
            else:
                d["NUME_LIGN"] = [j + 1 for j in range(len(vals)) if vals[j] is not None]
                d[mc] = [v for v in vals if v is not None]
            if len(d[mc]) == 0:
                UTMESS("I", "TABLE0_34", valk=self.para[i])
            else:
                dico["LISTE"].append(d)
        if len(dico["LISTE"]) == 0:
            UTMESS("A", "TABLE0_35")
        return dico

    def Array(self, Para, Champ):
        """Renvoie sous forme de NumArray le résultat d'une extraction dans une table
        méthode utile à macr_recal
        """
        __Rep = self[Para, Champ].values()
        F = numpy.zeros((len(__Rep[Para]), 2))
        for i in range(len(__Rep[Para])):
            F[i][0] = __Rep[Para][i]
            F[i][1] = __Rep[Champ][i]
        del __Rep
        return F

    def Croise(self):
        """Retourne un tableau croisé P3(P1,P2) à partir d'une table ayant
        trois paramètres (P1, P2, P3).
        """
        if len(self.para) != 3:
            UTMESS("A", "TABLE0_36")
            return Table()
        py, px, pz = self.para
        ly, lx, lz = [list(getattr(self, p).values()) for p in self.para]
        new_rows = []
        # lpz='%s=f(%s,%s)' % (pz,px,py)
        lpz = "%s/%s" % (px, py)
        # attention aux doublons dans lx et ly
        new_para = set(ly)
        new_para.discard(None)
        new_para = list(new_para)
        new_para.sort()
        new_para.insert(0, lpz)
        # attention aux doublons dans lx et ly
        newx = set(lx)
        newx.discard(None)
        newx = list(newx)
        newx.sort()
        for x in newx:
            if x is not None:
                d = {lpz: x}
                taux = getattr(self, px) == x
                for dz in taux.rows:
                    d[dz[py]] = dz[pz]
                new_rows.append(d)
        new_type = [self.type[1]] + [self.type[2]] * (len(new_para) - 1)
        new_titr = self.titr
        if new_titr != "":
            new_titr += "\n"
        new_titr += pz + " FONCTION DE " + px + " ET " + py
        return Table(new_rows, new_para, new_type, new_titr)

    def Renomme(self, pold, pnew):
        """Renomme le paramètre `pold` en `pnew`."""
        if not pold in self.para:
            raise KeyError("Paramètre %s inexistant dans cette table" % pold)
        elif self.para.count(pnew) > 0:
            raise KeyError("Le paramètre %s existe déjà dans la table" % pnew)
        else:
            self.para[self.para.index(pold)] = pnew
            for lig in self:
                if lig.get(pold) is not None:
                    lig[pnew] = lig[pold]
                    del lig[pold]


class Colonne(TableBase):

    """Classe intermédiaire pour mémoriser un couple (table, nom de colonne)
    et exprimer les critères d'extraction sous une forme naturelle en python
    en surchargeant les operateurs <, >, != et =.
    Alors on peut écrire la requete simple :
      soustable=t.a<10
    Ainsi que des requetes plus complexes :
      soustable=t.a<10 and t.b <4
    ou
      soustable=t.a<10 or t.b <4
    Les "alias" EQ, NE, LE, LT, GE, GT permettent à la macro IMPR_TABLE
    d'utiliser directement le mot-clé utilisateur CRIT_COMP défini dans le
    catalogue : getattr(Table,CRIT_COMP).
    """

    def __init__(self, table, column, typ=None):
        """Constructeur (objet Table associé, paramètre de la colonne, type du
        paramètre).
        """
        self.Table = table
        self.rows = self.Table.rows
        self.para = column
        self.type = typ
        self.titr = ""

    def __getstate__(self):
        """Return internal state.

        Returns:
            tuple: Internal state.
        """
        return self.Table, self.para, self.type, self.titr

    def __setstate__(self, state):
        """Restore internal state.

        Arguments:
            state (tuple): Internal state.
        """
        if isinstance(state, tuple) and len(state) == 4:
            self.Table, self.para, self.type, self.titr = state
            self.rows = self.Table.rows

    def _extract(self, fun):
        """Construit une table avec les lignes de self.Table
        dont l'élément de nom self.para satisfait le critère fun,
        fun est une fonction qui retourne vrai ou faux
        """
        return Table(
            [row for row in self.Table if fun(row.get(self.para))],
            self.Table.para,
            self.Table.type,
            self.Table.titr,
        )

    def __le__(self, VALE):
        if is_sequence(VALE):
            crit = max(VALE)
        else:
            crit = VALE
        return self._extract(lambda v: v is not None and v <= crit)

    def __lt__(self, VALE):
        if is_sequence(VALE):
            crit = max(VALE)
        else:
            crit = VALE
        return self._extract(lambda v: v is not None and v < crit)

    def __ge__(self, VALE):
        if is_sequence(VALE):
            crit = min(VALE)
        else:
            crit = VALE
        return self._extract(lambda v: v is not None and v >= crit)

    def __gt__(self, VALE):
        if is_sequence(VALE):
            crit = min(VALE)
        else:
            crit = VALE
        return self._extract(lambda v: v is not None and v > crit)

    def __eq__(self, VALE, CRITERE="RELATIF", PRECISION=0.0):
        if not is_sequence(VALE):
            VALE = [VALE]
        if is_str(VALE[0]):
            stripVALE = [value.strip() for value in VALE]
            return self._extract(lambda v: str(v).strip() in stripVALE)
        else:
            if PRECISION == 0.0:
                return self._extract(lambda v: v in VALE)
            elif CRITERE == "ABSOLU":
                return self._extract(lambda v: _func_test_abs(v, VALE, PRECISION))
            else:
                return self._extract(lambda v: _func_test_rela(v, VALE, PRECISION))

    def REGEXP(self, regexp):
        """Retient les lignes dont le paramètre satisfait l'expression
        régulière `regexp`.
        """
        if not is_str(regexp):
            return self._extract(lambda v: False)
        return self._extract(lambda v: v is not None and re.search(regexp, v) is not None)

    def __ne__(self, VALE, CRITERE="RELATIF", PRECISION=0.0):
        if not is_sequence(VALE):
            VALE = [VALE]
        if is_str(VALE[0]):
            stripVALE = [value.strip() for value in VALE]
            return self._extract(lambda v: str(v).strip() not in stripVALE)
        else:
            if PRECISION == 0.0:
                return self._extract(lambda v: v not in VALE)
            elif CRITERE == "ABSOLU":
                return self._extract(lambda v: not (_func_test_abs(v, VALE, PRECISION)))
            else:
                return self._extract(lambda v: not (_func_test_rela(v, VALE, PRECISION)))

    def MAXI(self):
        # important pour les performances de récupérer le max une fois pour
        # toutes
        maxi = max(self)
        return self._extract(lambda v: v == maxi)

    def MINI(self):
        # important pour les performances de récupérer le min une fois pour
        # toutes
        mini = min(self)
        return self._extract(lambda v: v == mini)

    def MAXI_ABS(self):
        # important pour les performances de récupérer le max une fois pour
        # toutes
        maxi_abs = max([abs(v) for v in list(self.values()) if is_number(v)])
        return self._extract(lambda v: v == maxi_abs or v == -maxi_abs)

    def MINI_ABS(self):
        # important pour les performances de récupérer le min une fois pour
        # toutes
        mini_abs = min([abs(v) for v in list(self.values()) if is_number(v)])
        # tester le type de v est trop long donc pas de abs(v)
        return self._extract(lambda v: v == mini_abs or v == -mini_abs)

    def __iter__(self):
        """Itère sur les éléments de la colonne"""
        for row in self.Table:
            # si l'élément n'est pas présent on retourne None
            yield row.get(self.para)
            # yield row[self.para]

    def __getitem__(self, i):
        """Retourne la ième valeur d'une colonne"""
        return list(self.values())[i]

    def values(self):
        """Renvoie la liste des valeurs"""
        return [r.get(self.para, None) for r in self.Table]

    def not_none_values(self):
        """Renvoie la liste des valeurs non 'None'"""
        return [val for val in list(self.values()) if val is not None]

    # équivalences avec les opérateurs dans Aster
    LE = __le__
    LT = __lt__
    GE = __ge__
    GT = __gt__
    EQ = __eq__
    NE = __ne__

    def VIDE(self):
        return self.__eq__(None)

    def NON_VIDE(self):
        return self.__ne__(None)


def sort_table(rows, l_para, l_type, w_para, reverse=False):
    """Sort list of dict.
    rows     : list of dict
    l_para   : list of the keys of dict
    t_type   : list of the types of each parameter
    w_para   : keys of the sort
    """
    orows = sorted([OrderableRow(row, l_para, l_type, w_para) for row in rows], reverse=reverse)
    new_rows = [orow.data for orow in orows]
    return new_rows


class OrderableRow:
    """A row that can be compared to another.

    Arguments:
        data (dict): dict of the values.
        parameters (list[str]): parameters names of the Table.
        types (list[str]): types of the parameters of the Table.
        keys (list[str]): sub-list of the parameters to be used for ordering.
    """

    defval = {"I": 0, "R": 0.0}

    def __init__(self, data, parameters, types, keys):
        self.data = data
        self.keys = keys
        self.defaults = dict()
        for typ, para in zip(types, parameters):
            self.defaults[para] = OrderableRow.defval.get(typ, "")

    def __lt__(self, other):
        """Tell if the row must be upper than the 'other' one."""
        # 'keys' attr should be the same one
        for col in self.keys:
            default = self.defaults[col]
            left = self.data.get(col, default)
            left = left if left is not None else default
            right = other.data.get(col, default)
            right = right if right is not None else default
            if left < right:
                return True
            if left > right:
                return False
        return True


def FMT(dform, nform, typAster=None, larg=0, val=""):
    """Retourne un format d'impression Python à partir d'un type Aster ('R','I',
    'K8', 'K16'...). Si typAster is None, retourne dform[nform].
       larg : largeur minimale du format (val permet de ne pas ajouter des blancs
       si la chaine à afficher est plus longue que le format, on prend le partie
       de ne pas tronquer les chaines)
    """
    if typAster is None:
        fmt = dform[nform]
    elif typAster in ("I", "R"):
        if nform == "formK":
            # convertit %12.5E en %-12s
            fmt = re.sub(r"([0-9]+)[\.0-9]*[diueEfFgG]+", r"-\g<1>s", dform["form" + typAster])
        else:
            fmt = dform[nform]
    else:
        # typAster = Kn
        fmt = "%-" + typAster[1:] + "s"
    # on ajoute éventuellement des blancs pour atteindre la largeur demandée
    if larg != 0:
        fmt = " " * max(min(larg - len(val), larg - len(fmt % 0)), 0) + fmt
    return fmt


def merge(tab1, tab2, labels=[], restrict=False, format_r=None):
    """Assemble les deux tables tb1 et tb2 selon une liste de labels communs.
    Si labels est vide:
     - les lignes de tb2 sont ajoutés à celles de tb1,
    sinon :
     - si on trouve les valeurs de tb2 sur les labels dans tb1 (et une seule fois),
       on surcharge tb1 avec les lignes de tb2 ;
     - sinon on ajoute la ligne de tb2 à la fin de tb1.
    Si format_r est fourni, on "arrondi" les réels selon ce format lors de la
    comparaison.
    """
    _reformat = _reformat_func(format_r)
    tb1 = tab1.copy()
    tb2 = tab2.copy()
    if not is_sequence(labels):
        labels = (labels,)
    for key in labels:
        if key not in tb1.para:
            UTMESS("F", "TABLE0_27", valk=key)
        if key not in tb2.para:
            UTMESS("F", "TABLE0_27", valk=key)
    # ensemble des paramètres et des types
    n_para = tb1.para[:]
    n_type = tb1.type[:]
    for i in tb2.para:
        if i not in tb1.para:
            n_para.append(i)
            n_type.append(tb2.type[tb2.para.index(i)])
    # restriction des lignes aux labels communs (peu cher en cpu)
    dlab1 = _unique_values(tb1, labels, format_r)
    dlab2 = _unique_values(tb2, labels, format_r)
    # creation de dic1 : dictionnaire de correspondance entre les
    # lignes a merger dans les deux tableaux
    dic1 = {}
    for cle in list(dlab1.keys()):
        if dlab1[cle] is None or cle == ():
            bid = dlab1.pop(cle)
    for cle in list(dlab2.keys()):
        if dlab2[cle] is None or cle == ():
            bid = dlab2.pop(cle)
    for cle in list(dlab2.keys()):
        if cle in dlab1:
            dic1[dlab2[cle]] = dlab1[cle]
    # insertion des valeurs de tb2 dans tb1 quand les labels sont communs
    # (et uniques dans chaque table)
    # OU ajout de la ligne de tb2 dans tb1 (si restrict == False)
    rows1 = tab1.rows
    rows2 = tab2.rows
    if restrict:

        def func_append_r2(row):
            pass

    else:

        def func_append_r2(row):
            rows1.append(row)

    i2 = -1
    for r2 in rows2:
        i2 += 1
        try:
            rows1[dic1[i2]].update(r2)
        except KeyError:
            func_append_r2(r2)
    # concaténation des titres + info sur le merge
    tit = "\n".join([tb1.titr, tb2.titr, "MERGE avec labels=%s" % repr(labels)])
    return Table(rows1, n_para, n_type, tit)


def remove_twins(tab, labels=[], format_r=None):
    """Remove lines if the values of `labels` have already been seen"""
    kept = _unique_values(tab, labels, format_r, keep_first=True)
    removed = set(range(len(tab))).difference(list(kept.values()))
    todel = reversed(sorted(list(removed)))
    for i in todel:
        del tab.rows[i]


def typaster(obj, prev=None, strict=False):
    """Retourne le type Aster ('R', 'I', Kdef) correspondant à l'objet obj.
    Si prev est fourni, on vérifie que obj est du type prev.
    Si strict=False, on autorise que obj ne soit pas du type prev s'ils sont
    tous les deux numériques ; dans ce cas, on retourne le "type enveloppe" 'R'.
    """
    dtyp = {int: "I", float: "R", str: Kdef, type(None): "I"}
    if is_float(obj):
        obj = float(obj)
    if type(obj) in list(dtyp.keys()):
        typobj = dtyp[type(obj)]
        if prev in [None, typobj]:
            return typobj
        elif prev[0] == typobj[0] == "K":
            if len(obj) <= int(prev[1:]):
                return prev
            else:
                raise TypeError(
                    "La longueur de la chaine %s est incompatible avec le type %s"
                    % (repr(obj), repr(prev))
                )
        elif strict:  # prev is not None et typobj != prev et strict
            raise TypeError("La valeur %s n'est pas de type %s" % (repr(obj), repr(prev)))
        elif prev in ("I", "R") and typobj in ("I", "R"):
            return "R"
        else:
            raise TypeError(
                "La valeur %s n'est pas compatible avec le type %s" % (repr(obj), repr(prev))
            )
    else:
        raise TypeError(
            "Une table ne peut contenir que des entiers, réels " "ou chaines de caractères."
        )


def _unique_values(tab, params, format_r, keep_first=False):
    """Return a dict where:
    - the keys are tuples of values restricted to some columns
    - the values are the line_index
    index is set to None if several lines have the same values
    or, if keep_first is True, index is the first line index found.
    """
    _reformat = _reformat_func(format_r)
    rows = tab.rows
    dlab = {}
    for i in range(len(rows)):
        tup = _reformat(list(map(rows[i].__getitem__, params)))
        if dlab.get(tup, "") == "":
            dlab[tup] = i
        elif not keep_first:
            dlab[tup] = None
    return dlab


# fonctions utilitaires
def _reformat_func(format_r):
    """Return a function to format a list of values"""
    if format_r:
        fmtr = fmtF2PY(format_r)

        def _func(values):
            """Convertit les réels en utilisant format_r"""
            conv = []
            for i in values:
                try:
                    vali = fmtr % i
                except TypeError:
                    vali = i
                conv.append(vali)
            return tuple(conv)

    else:

        def _func(values):
            """no change"""
            return tuple(values)

    return _func


def _func_test_abs(v, VALE, PRECISION):
    """Retourne True si v est parmi VALE à PRECISION près en absolu"""
    for x in VALE:
        if v is not None and (x - PRECISION <= v <= x + PRECISION):
            return True
    return False


def _func_test_rela(v, VALE, PRECISION):
    """Retourne True si v est parmi VALE à PRECISION près en relatif"""
    for x in VALE:
        sign = float(x > 0.0) or -1.0
        if v is not None and (
            sign * x * (1.0 - PRECISION) <= sign * v <= sign * x * (1.0 + PRECISION)
        ):
            return True
    return False


def check_nan(values):
    """Raise ValueError exception if nan is found in values."""
    bool_val = numpy.isnan([p if p is not None else 0.0 for p in values])
    if bool_val.any():
        indexes = [str(idx[0]) for idx in numpy.argwhere(bool_val)]
        raise ValueError("NaN present at indexes %s" % ",".join(indexes))


def check_inf(values):
    """Raise ValueError exception if inf is found in values."""
    bool_val = numpy.isinf([p if p is not None else 0.0 for p in values])
    if bool_val.any():
        indexes = [str(idx[0]) for idx in numpy.argwhere(bool_val)]
        raise ValueError("Inf present at indexes %s" % ",".join(indexes))
