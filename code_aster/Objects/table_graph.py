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
__all__ = ["Graph", "AjoutParaCourbe"]

# TODO solve dependency TablePy/Graph
# aslint: disable=C4008

import os
import os.path
import re
import shutil
import sys
import time

import numpy as np

from ..Messages import UTMESS
from ..Utilities import ExecutionParameter, disable_fpe, value_is_sequence

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class Graph:
    """Cette classe définit l'objet Graph pour Code_Aster.

    Important :  Utiliser les méthodes dédiées à la manipulation des données
       (AjoutCourbe, ...) car elles tiennent à jour les attributs "privés"
       relatifs aux données : NbCourbe, les extrema...

    Attributs :
    - Données de chaque courbe :
       .Valeurs   : liste des valeurs de chaque courbe, pour chaque courbe :
          (paramètres, parties réelles [, parties imaginaires])
       .Legendes  : liste des noms de chaque courbe
       .Labels    : liste des noms des colonnes de chaque courbe
       .Styles    : liste des infices de styles de ligne
       .Couleurs  : liste des indices de couleurs
       .Marqueurs : liste des indices de symboles/marqueurs
       .FreqMarq  : liste des fréquences des marqueurs
       .Tri       : liste du tri à effectuer sur les données ('N', 'X', 'Y',
          'XY' ou 'YX')
      Pour Lignes, Couleurs, Marqueurs, FreqMarq, -1 signifie valeur par défaut
      du traceur.

    - Propriétés :
       .Titre : titre du graphique
       .SousTitre : sous-titre (appelé commentaire dans agraf)
    - Axes :
       .Min_X, .Max_X, .Min_Y, .Max_Y : bornes du tracé (une méthode permet de
          renseigner automatiquement ces valeurs avec les extréma globaux)
       .Legende_X, .Legende_Y : légende des axes
       .Echelle_X, .Echelle_Y : type d'échelle (LIN, LOG)
       .Grille_X, .Grille_Y : paramètre de la grille (pas ou fréquence au choix
          de l'utilisateur en fonction du traceur qu'il veut utiliser)

    Attributs privés (modifiés uniquement par les méthodes de la classe) :
       .NbCourbe : nombre de courbes
       .BBXmin, BBXmax, BBYmin, BBYmax : extrema globaux (bounding box)
       .LastTraceArgs, LastTraceFormat : données utilisées lors du dernier tracé
    """

    # ------------------------------------------------------------------------

    def __init__(self):
        """Construction + valeurs par défaut des attributs"""
        self.Valeurs = []
        self.Legendes = []
        self.Labels = []
        self.Styles = []
        self.Couleurs = []
        self.Marqueurs = []
        self.FreqMarq = []
        self.Tri = []
        self.Titre = ""
        self.SousTitre = ""
        self.Min_X = None
        self.Max_X = None
        self.Min_Y = None
        self.Max_Y = None
        self.MinP_X = 1.0e99  # minimum > 0 pour les échelles LOG
        self.MinP_Y = 1.0e99
        self.Legende_X = ""
        self.Legende_Y = ""
        self.Echelle_X = "LIN"
        self.Echelle_Y = "LIN"
        self.Grille_X = -1
        self.Grille_Y = -1
        # attributs que l'utilisateur ne doit pas modifier
        self.NbCourbe = len(self.Valeurs)
        self.BBXmin = 1.0e99
        self.BBXmax = -1.0e99
        self.BBYmin = 1.0e99
        self.BBYmax = -1.0e99
        # pour conserver les paramètres du dernier tracé
        self.LastTraceArgs = {}
        self.LastTraceFormat = ""

    def __get_titre(self):
        """private get method"""
        return self._titre

    def __set_titre(self, value):
        """private set method"""
        if type(value) not in (list, tuple):
            value = [value]
        self._titre = value

    Titre = property(__get_titre, __set_titre)

    # ------------------------------------------------------------------------
    def SetExtremaX(self, marge=0.0, x0=None, x1=None, force=True):
        """Remplit les limites du tracé (Min/Max_X) avec les valeurs de la
        bounding box +/- avec une 'marge'*(Max-Min)/2.
        x0,x1 permettent de modifier la bb.
        """
        if x0 is not None:
            self.BBXmin = min([self.BBXmin, x0])
        if x1 is not None:
            self.BBXmax = max([self.BBXmax, x1])

        dx = max(self.BBXmax - self.BBXmin, 0.01 * self.BBXmax)
        if dx == 0.0:
            dx = 1.0e-6
        if force or self.Min_X is None:
            self.Min_X = self.BBXmin - marge * dx / 2.0
        if force or self.Max_X is None:
            self.Max_X = self.BBXmax + marge * dx / 2.0
        return

    def SetExtremaY(self, marge=0.0, y0=None, y1=None, force=True):
        """Remplit les limites du tracé (Min/Max_Y) avec les valeurs de la
        bounding box +/- avec une 'marge'*(Max-Min)/2.
        y0,y1 permettent de modifier la bb.
        """
        if y0 is not None:
            self.BBYmin = min([self.BBYmin, y0])
        if y1 is not None:
            self.BBYmax = max([self.BBYmax, y1])

        dy = max(self.BBYmax - self.BBYmin, 0.01 * self.BBYmax)
        if dy == 0.0:
            dy = 1.0e-6
        if force or self.Min_Y is None:
            self.Min_Y = self.BBYmin - marge * dy / 2.0
        if force or self.Max_Y is None:
            self.Max_Y = self.BBYmax + marge * dy / 2.0
        return

    def SetExtrema(self, marge=0.0, x0=None, x1=None, y0=None, y1=None, force=True):
        """Remplit les limites du tracé (Min/Max_X/Y) avec les valeurs de la
        bounding box +/- avec une 'marge'*(Max-Min)/2.
        x0,x1,y0,y1 permettent de modifier la bb.
        """
        self.SetExtremaX(marge, x0, x1, force=force)
        self.SetExtremaY(marge, y0, y1, force=force)
        return

    # ------------------------------------------------------------------------

    def AutoBB(self, debut=-1):
        """Met à jour automatiquement la "bounding box"
        (extrema toutes courbes confondues)
        Appelé par les méthodes de manipulation des données
        """
        if debut == -1:
            debut = self.NbCourbe - 1
        if debut == 0:
            X0 = 1.0e99
            X1 = -1.0e99
            Y0 = 1.0e99
            Y1 = -1.0e99
        else:
            X0 = self.BBXmin
            X1 = self.BBXmax
            Y0 = self.BBYmin
            Y1 = self.BBYmax

        for i in range(debut, self.NbCourbe):
            X0 = min([X0] + list(self.Valeurs[i][0]))
            X1 = max([X1] + list(self.Valeurs[i][0]))
            self.MinP_X = min([self.MinP_X] + [x for x in list(self.Valeurs[i][0]) if x > 0])
            for ny in range(1, len(self.Valeurs[i])):
                Y0 = min([Y0] + list(self.Valeurs[i][ny]))
                Y1 = max([Y1] + list(self.Valeurs[i][ny]))
                self.MinP_Y = min([self.MinP_Y] + [y for y in list(self.Valeurs[i][ny]) if y > 0])
        self.BBXmin = X0
        self.BBXmax = X1
        self.BBYmin = Y0
        self.BBYmax = Y1
        return

    # ------------------------------------------------------------------------

    def AjoutCourbe(self, Val, Lab, Leg="", Sty=-1, Coul=-1, Marq=-1, FreqM=-1, Tri="N"):
        """Ajoute une courbe dans les données
           Val   : liste de 2 listes (ou 3 si complexe) : abs, ord[, imag]
           Leg   : une chaine
           Lab   : liste de 2 chaines (ou 3 si complexe)
           Sty   : un entier
           Coul  : un entier
           Marq  : un entier
           FreqM : un entier
           Tri   : chaine de caractères : N, X, Y, XY ou YX
        Met à jour les attributs : NbCourbe, BBXmin/Xmax/Ymin/Ymax
        """
        nbc = len(Val)  # nombre de colonnes : 2 ou 3

        # verifications : "if not (conditions requises)"
        if not (
            2 <= nbc <= 3
            and value_is_sequence(Val[0])
            and value_is_sequence(Val[1])
            and (nbc == 2 or value_is_sequence(Val[2]))
            and len(Val[0]) == len(Val[1])
            and (nbc == 2 or len(Val[0]) == len(Val[2]))
        ):
            UTMESS("F", "GRAPH0_1", valk="Val")

        if len(Lab) != nbc:
            UTMESS("S", "GRAPH0_2", valk="Lab")

        # ajout dans les données
        self.Legendes.append(str(Leg))
        self.Labels.append([str(L) for L in Lab])
        self.Valeurs.append(Val)
        self.Styles.append(Sty)
        self.Couleurs.append(Coul)
        self.Marqueurs.append(Marq)
        self.FreqMarq.append(FreqM)
        self.Tri.append(Tri)

        self.NbCourbe = self.NbCourbe + 1
        self.AutoBB()
        return

    # ------------------------------------------------------------------------

    def Courbe(self, n):
        """Permet de récupérer les données de la courbe d'indice n sous forme
        d'un dictionnaire.
        """
        dico = {
            "Leg": self.Legendes[n],  # légende de la courbe
            "LabAbs": self.Labels[n][0],  # labels des abscisses
            "LabOrd": [self.Labels[n][1]],  # labels des ordonnées
            "NbCol": len(self.Valeurs[n]),  # nombre de colonnes
            "NbPts": len(self.Valeurs[n][0]),  # nombre de points
            "Abs": self.Valeurs[n][0],  # liste des abscisses
            "Ord": [self.Valeurs[n][1]],  # liste des ordonnées
            "Sty": self.Styles[n],  # style de la ligne
            "Coul": self.Couleurs[n],  # couleur
            "Marq": self.Marqueurs[n],  # marqueur
            "FreqM": self.FreqMarq[n],  # fréquence du marqueur
            "Tri": self.Tri[n],  # ordre de tri des données
        }
        if dico["NbCol"] == 3:
            dico["LabOrd"].append(self.Labels[n][2])  # labels de la partie imaginaire
            dico["Ord"].append(self.Valeurs[n][2])
            # liste des ordonnées partie imaginaire
        return dico

    # ------------------------------------------------------------------------

    def Trace(self, FICHIER=None, FORMAT=None, dform=None, **opts):
        """Tracé du Graph selon le format spécifié.
        FICHIER : nom du(des) fichier(s). Si None, on dirige vers stdout
        dform : dictionnaire de formats d'impression (format des réels,
           commentaires, saut de ligne...)
        opts  : voir TraceGraph.
        """
        para = {
            "TABLEAU": {"mode": "a", "driver": TraceTableau},
            "XMGRACE": {"mode": "a", "driver": TraceXmgrace},
            "AGRAF": {"mode": "a", "driver": TraceAgraf},
            "LISS_ENVELOP": {"mode": "a", "driver": TraceMatplotlib},
        }
        kargs = {}
        if self.LastTraceArgs == {}:
            kargs["FICHIER"] = FICHIER
            kargs["dform"] = dform
            kargs["opts"] = opts
        else:
            kargs = self.LastTraceArgs.copy()
            if FORMAT is None:
                FORMAT = self.LastTraceFormat
            if FICHIER is not None:
                kargs["FICHIER"] = FICHIER
            if dform is not None:
                kargs["dform"] = dform
            if opts != {}:
                kargs["opts"] = opts
        if not FORMAT in list(para.keys()):
            UTMESS("A", "GRAPH0_3", valk=FORMAT)
        else:
            kargs["fmod"] = para[FORMAT]["mode"]
            self.LastTraceArgs = kargs.copy()
            self.LastTraceFormat = FORMAT
            # call the associated driver
            para[FORMAT]["driver"](self, **kargs)

    # ------------------------------------------------------------------------

    def __repr__(self):
        """Affichage du contenu d'un Graph"""
        srep = ""
        for attr in [
            "NbCourbe",
            "Legendes",
            "Labels",
            "Valeurs",
            "Min_X",
            "Max_X",
            "Min_Y",
            "Max_Y",
            "BBXmax",
            "BBXmin",
            "BBYmax",
            "BBYmin",
            "Legende_X",
            "Legende_Y",
            "Echelle_X",
            "Echelle_Y",
            "Grille_X",
            "Grille_Y",
            "Tri",
        ]:
            srep = srep + "%-10s : %s\n" % (attr, str(getattr(self, attr)))
        return srep


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------


class TraceGraph:
    """
    Cette classe définit le tracé d'un objet Graph dans un fichier.

    Attributs :
       .NomFich : liste de noms de fichier de sortie

    Attributs privés (modifiés uniquement par les méthodes de la classe) :
       .Fich    : liste des objets 'fichier'
       .Graph   : objet Graph que l'on veut tracer
       .DicForm : dictionnaire des formats de base (séparateur, format des réels...)

    Les méthodes Entete, DescrCourbe, Trace (définition de l'entete, partie descriptive
    d'une courbe, méthode de tracé/impression) sont définies dans une classe dérivée.
    """

    # ------------------------------------------------------------------------

    def __init__(self, graph, FICHIER, fmod="w", dform=None, opts={}):
        """Construction, ouverture du fichier, surcharge éventuelle du formatage
        (dform), mode d'ouverture du fichier (fmod).
        opts  : dictionnaire dont les valeurs seront affectées comme attributs
           de l'objet (A utiliser pour les propriétés spécifiques
           à un format, exemple 'PILOTE' pour Xmgrace).
        """
        # attributs optionnels (au début pour éviter un écrasement maladroit !)
        for k, v in list(opts.items()):
            setattr(self, k, v)

        # Ouverture du(des) fichier(s)
        self.NomFich = []
        if type(FICHIER) is bytes or type(FICHIER) is str:
            self.NomFich.append(FICHIER)
        elif type(FICHIER) in (list, tuple):
            self.NomFich = FICHIER[:]
        else:
            # dans ce cas, on écrira sur stdout (augmenter le 2 éventuellement)
            self.NomFich = [None] * 2
        self.Fich = []
        for ff in self.NomFich:
            if ff is not None:
                self.Fich.append(open(ff, fmod))
            else:
                self.Fich.append(sys.stdout)

        # objet Graph sous-jacent
        self.Graph = graph
        # si Min/Max incohérents
        if graph.Min_X is None or graph.Max_X is None or graph.Min_X > graph.Max_X:
            graph.SetExtremaX(marge=0.05, force=True)
        if graph.Min_Y is None or graph.Max_Y is None or graph.Min_Y > graph.Max_Y:
            graph.SetExtremaY(marge=0.05, force=True)

        if graph.Echelle_X == "LOG":
            graph.Grille_X = 10
            # verif si Min<0 à cause de la marge
            if graph.Min_X < 0.0:
                if graph.BBXmin < 0.0:
                    UTMESS("A", "GRAPH0_4")
                graph.Min_X = graph.MinP_X
        if graph.Echelle_Y == "LOG":
            graph.Grille_Y = 10
            if graph.Min_Y < 0.0:
                if graph.BBYmin < 0.0:
                    UTMESS("A", "GRAPH0_5")
                graph.Min_Y = graph.MinP_Y

        # formats de base (identiques à ceux du module Table)
        self.DicForm = {
            "csep": " ",  # séparateur
            "ccom": "#",  # commentaire
            "ccpara": "",  # commentaire des labels
            "cdeb": "",  # début de ligne
            "cfin": "\n",  # fin de ligne
            "sepch": ";",  # remplace les sauts de ligne à l'intérieur d'une cellule
            "formK": "%-12s",  # chaines
            "formR": "%12.5E",  # réels
            "formI": "%12d",  # entiers
        }
        if dform is not None and type(dform) == dict:
            self.DicForm.update(dform)

        # let's go
        self.Trace()

    # ------------------------------------------------------------------------
    def __del__(self):
        """Fermeture du(des) fichier(s) à la destruction"""
        if hasattr(self, "Fich"):
            self._FermFich()

    # ------------------------------------------------------------------------

    def _FermFich(self):
        """Fermeture du(des) fichier(s)"""
        for fp in self.Fich:
            if fp != sys.stdout:
                fp.close()

    # ------------------------------------------------------------------------

    def _OuvrFich(self):
        """Les fichiers sont ouverts par le constructeur. S'ils ont été fermés,
        par un appel au Tracé, _OuvrFich ouvre de nouveau les fichiers dans le
        même mode"""
        n = len(self.NomFich)
        for i in range(n):
            if self.Fich[i].closed:
                self.Fich[i] = open(self.NomFich[i], self.Fich[i].mode)

    # ------------------------------------------------------------------------
    def Entete(self):
        """Retourne l'entete"""
        raise NotImplementedError("Cette méthode doit être définie par la classe fille.")

    # ------------------------------------------------------------------------

    def DescrCourbe(self, **args):
        """Retourne la chaine de caractères décrivant les paramètres de la courbe."""
        raise NotImplementedError("Cette méthode doit être définie par la classe fille.")

    # ------------------------------------------------------------------------

    def Trace(self):
        """Méthode pour 'tracer' l'objet Graph dans un fichier.
        Met en page l'entete, la description des courbes et les valeurs selon
        le format et ferme le fichier.
        """
        raise NotImplementedError("Cette méthode doit être définie par la classe fille.")


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------
class TraceTableau(TraceGraph):
    """
    Impression d'un objet Graph sous forme d'un tableau de colonnes,
    on suppose que les courbes partagent la même liste d'abscisse à 'EPSILON'
    près, sinon on alarme.
    """

    EPSILON = 1.0e-4

    def Trace(self):
        """Méthode pour 'tracer' l'objet Graph dans un fichier.
        Met en page l'entete, la description des courbes et les valeurs selon
        le format et ferme le fichier.
        L'ouverture et la fermeture du fichier sont gérées par l'objet TablePy.
        """
        from .table_py import Table as TablePy

        g = self.Graph
        msg = []
        if g.NbCourbe > 0:
            # validité des données (abscisses identiques)
            t0 = np.array(g.Courbe(0)["Abs"])
            max0 = max(abs(t0))
            for i in range(1, g.NbCourbe):
                if g.Courbe(i)["NbPts"] != g.Courbe(0)["NbPts"]:
                    msg.append(
                        "La courbe %d n'a pas le même " "nombre de points que la 1ère." % (i + 1)
                    )
                else:
                    ti = np.array(g.Courbe(i)["Abs"])
                    if max(abs((ti - t0).ravel())) > self.EPSILON * max0:
                        msg.append(
                            "Courbe %d : écart entre les "
                            "abscisses supérieur à %9.2E" % (i + 1, self.EPSILON)
                        )
                        msg.append(
                            "     Utilisez IMPR_FONCTION pour interpoler "
                            "les valeurs sur la première liste d'abscisses."
                        )
            # objet TablePy
            Tab = TablePy()
            # titre / sous-titre
            tit = []
            tit.extend([self.DicForm["ccom"] + " " + line for line in g.Titre])
            tit.append(self.DicForm["ccom"] + " " + g.SousTitre)
            # legendes
            for i in range(g.NbCourbe):
                tit.append(self.DicForm["ccom"] + " Courbe " + str(i))
                tit.extend(
                    [self.DicForm["ccom"] + " " + leg for leg in g.Legendes[i].split(os.linesep)]
                )
            Tab.titr = self.DicForm["cfin"].join(tit)
            # noms des paramètres/colonnes
            Tab.para.append(g.Labels[0][0])
            for i in range(g.NbCourbe):
                for lab in g.Labels[i][1:]:
                    Tab.para.append(lab)
            # types
            Tab.type = ["R"] * len(Tab.para)
            # lignes de la Table
            dC0 = g.Courbe(0)
            for j in range(dC0["NbPts"]):
                row = {}
                row[dC0["LabAbs"]] = dC0["Abs"][j]
                for i in range(g.NbCourbe):
                    dCi = g.Courbe(i)
                    for k in range(dCi["NbCol"] - 1):
                        try:
                            row[dCi["LabOrd"][k]] = dCi["Ord"][k][j]
                        except IndexError:
                            row[dCi["LabOrd"][k]] = None
                Tab.append(row)
            Tab.Impr(FICHIER=self.NomFich[0], FORMAT="TABLEAU", dform=self.DicForm)
            # erreurs ?
            if msg:
                UTMESS("A", "GRAPH0_6", valk="\n".join(msg))
        return


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------


class TraceXmgrace(TraceGraph):
    """
    Impression d'un objet Graph au format XMGRACE.
    Attribut supplémentaire : .PILOTE
    """

    PILOTE = ""
    # ------------------------------------------------------------------------

    def Entete(self):
        """Retourne l'entete du fichier .agr correspondant à la mise en forme"""
        dic_ech = {"LIN": "Normal", "LOG": "Logarithmic"}
        g = self.Graph
        entete = []
        entete.append(
            """
# Grace project file
#
@version 50100
@page size 842, 595
@page scroll 5%
@page inout 5%
@link page off
@map font 0 to "Times-Roman", "Times-Roman"
@map font 1 to "Times-Italic", "Times-Italic"
@map font 2 to "Times-Bold", "Times-Bold"
@map font 3 to "Times-BoldItalic", "Times-BoldItalic"
@map font 4 to "Helvetica", "Helvetica"
@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"
@map font 6 to "Helvetica-Bold", "Helvetica-Bold"
@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"
@map font 8 to "Courier", "Courier"
@map font 9 to "Courier-Oblique", "Courier-Oblique"
@map font 10 to "Courier-Bold", "Courier-Bold"
@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"
@map font 12 to "Symbol", "Symbol"
@map font 13 to "ZapfDingbats", "ZapfDingbats"
@map color 0 to (255, 255, 255), "white"
@map color 1 to (0, 0, 0), "black"
@map color 2 to (255, 0, 0), "red"
@map color 3 to (0, 255, 0), "green"
@map color 4 to (0, 0, 255), "blue"
@map color 5 to (255, 255, 0), "yellow"
@map color 6 to (188, 143, 143), "brown"
@map color 7 to (220, 220, 220), "grey"
@map color 8 to (148, 0, 211), "violet"
@map color 9 to (0, 255, 255), "cyan"
@map color 10 to (255, 0, 255), "magenta"
@map color 11 to (255, 165, 0), "orange"
@map color 12 to (114, 33, 188), "indigo"
@map color 13 to (103, 7, 72), "maroon"
@map color 14 to (64, 224, 208), "turquoise"
@map color 15 to (0, 139, 0), "green4"
@reference date 0
@date wrap off
@date wrap year 1950
@timestamp off
@default linewidth 1.0
@default linestyle 1
@default color 1
@default pattern 1
@default font 0
@default char size 1.000000
@default symbol size 1.000000
@default sformat "%.8g"
@background color 0
@page background fill on
@r0 off
@link r0 to g0
@r0 type above
@r0 linestyle 1
@r0 linewidth 1.0
@r0 color 1
@r0 line 0, 0, 0, 0
@r1 off
@link r1 to g0
@r1 type above
@r1 linestyle 1
@r1 linewidth 1.0
@r1 color 1
@r1 line 0, 0, 0, 0
@r2 off
@link r2 to g0
@r2 type above
@r2 linestyle 1
@r2 linewidth 1.0
@r2 color 1
@r2 line 0, 0, 0, 0
@r3 off
@link r3 to g0
@r3 type above
@r3 linestyle 1
@r3 linewidth 1.0
@r3 color 1
@r3 line 0, 0, 0, 0
@r4 off
@link r4 to g0
@r4 type above
@r4 linestyle 1
@r4 linewidth 1.0
@r4 color 1
@r4 line 0, 0, 0, 0
@g0 on
@g0 hidden false
@g0 type XY
@g0 stacked false
@g0 bar hgap 0.000000
@with g0
@    stack world 0, 0, 0, 0
@    znorm 1
@    view xmin 0.150000
@    view xmax 1.150000
@    view ymin 0.150000
@    view ymax 0.850000
@    title font 0
@    title size 1.500000
@    title color 1
@    subtitle font 0
@    subtitle size 1.000000
@    subtitle color 1
@    xaxes invert off
@    yaxes invert off
@    xaxis  on
@    xaxis  type zero false
@    xaxis  offset 0.000000 , 0.000000
@    xaxis  bar on
@    xaxis  bar color 1
@    xaxis  bar linestyle 1
@    xaxis  bar linewidth 1.0
@    xaxis  label layout para
@    xaxis  label place auto
@    xaxis  label char size 1.000000
@    xaxis  label font 0
@    xaxis  label color 1
@    xaxis  label place normal
@    xaxis  tick on
@    xaxis  tick minor ticks 1
@    xaxis  tick default 6
@    xaxis  tick place rounded true
@    xaxis  tick in
@    xaxis  tick major size 1.000000
@    xaxis  tick major color 1
@    xaxis  tick major linewidth 1.0
@    xaxis  tick major linestyle 2
@    xaxis  tick major grid on
@    xaxis  tick minor color 1
@    xaxis  tick minor linewidth 1.0
@    xaxis  tick minor linestyle 2
@    xaxis  tick minor grid off
@    xaxis  tick minor size 0.500000
@    xaxis  ticklabel on
@    xaxis  ticklabel format general
@    xaxis  ticklabel prec 5
@    xaxis  ticklabel angle 0
@    xaxis  ticklabel skip 0
@    xaxis  ticklabel stagger 0
@    xaxis  ticklabel place normal
@    xaxis  ticklabel offset auto
@    xaxis  ticklabel offset 0.000000 , 0.010000
@    xaxis  ticklabel start type auto
@    xaxis  ticklabel start 0.000000
@    xaxis  ticklabel stop type auto
@    xaxis  ticklabel stop 0.000000
@    xaxis  ticklabel char size 0.800000
@    xaxis  ticklabel font 0
@    xaxis  ticklabel color 1
@    xaxis  ticklabel formula ""
@    xaxis  ticklabel append ""
@    xaxis  ticklabel prepend ""
@    xaxis  tick place both
@    xaxis  tick spec type none
@    yaxis  on
@    yaxis  type zero false
@    yaxis  offset 0.000000 , 0.000000
@    yaxis  bar on
@    yaxis  bar color 1
@    yaxis  bar linestyle 1
@    yaxis  bar linewidth 1.0
@    yaxis  label layout para
@    yaxis  label place auto
@    yaxis  label char size 1.000000
@    yaxis  label font 0
@    yaxis  label color 1
@    yaxis  label place normal
@    yaxis  tick on
@    yaxis  tick minor ticks 1
@    yaxis  tick default 6
@    yaxis  tick place rounded true
@    yaxis  tick in
@    yaxis  tick major size 1.000000
@    yaxis  tick major color 1
@    yaxis  tick major linewidth 1.0
@    yaxis  tick major linestyle 2
@    yaxis  tick major grid on
@    yaxis  tick minor color 1
@    yaxis  tick minor linewidth 1.0
@    yaxis  tick minor linestyle 1
@    yaxis  tick minor grid off
@    yaxis  tick minor size 0.500000
@    yaxis  ticklabel on
@    yaxis  ticklabel format general
@    yaxis  ticklabel prec 5
@    yaxis  ticklabel angle 0
@    yaxis  ticklabel skip 0
@    yaxis  ticklabel stagger 0
@    yaxis  ticklabel place normal
@    yaxis  ticklabel offset auto
@    yaxis  ticklabel offset 0.000000 , 0.010000
@    yaxis  ticklabel start type auto
@    yaxis  ticklabel start 0.000000
@    yaxis  ticklabel stop type auto
@    yaxis  ticklabel stop 0.000000
@    yaxis  ticklabel char size 0.800000
@    yaxis  ticklabel font 0
@    yaxis  ticklabel color 1
@    yaxis  ticklabel formula ""
@    yaxis  ticklabel append ""
@    yaxis  ticklabel prepend ""
@    yaxis  tick place both
@    yaxis  tick spec type none
@    altxaxis  off
@    altyaxis  off
@    legend on
@    legend loctype view
@    legend 0.85, 0.8
@    legend box color 1
@    legend box pattern 1
@    legend box linewidth 1.0
@    legend box linestyle 1
@    legend box fill color 0
@    legend box fill pattern 1
@    legend font 0
@    legend char size 0.750000
@    legend color 1
@    legend length 4
@    legend vgap 1
@    legend hgap 1
@    legend invert false
@    frame type 0
@    frame linestyle 1
@    frame linewidth 1.0
@    frame color 1
@    frame pattern 1
@    frame background color 0
@    frame background pattern 0
"""
        )
        entete.append('@    title "' + " ".join(g.Titre) + '"')
        entete.append('@    subtitle "' + g.SousTitre + '"')
        entete.append('@    xaxis  label "' + g.Legende_X + '"')
        entete.append('@    yaxis  label "' + g.Legende_Y + '"')
        entete.append("@    xaxes scale " + dic_ech[g.Echelle_X])
        entete.append("@    yaxes scale " + dic_ech[g.Echelle_Y])
        entete.append("@    xaxis  tick major " + str(g.Grille_X))
        entete.append("@    yaxis  tick major " + str(g.Grille_Y))
        entete.append("@    world xmin " + str(g.Min_X))
        entete.append("@    world xmax " + str(g.Max_X))
        entete.append("@    world ymin " + str(g.Min_Y))
        entete.append("@    world ymax " + str(g.Max_Y))
        return entete

    # ------------------------------------------------------------------------

    def DescrCourbe(self, **args):
        """Retourne la chaine de caractères décrivant les paramètres de la courbe."""
        # valeurs par défaut
        sty = str(ValCycl(args["Sty"], 0, 8, 1))
        color = str(ValCycl(args["Coul"], 1, 15, args["NumSet"] + 1))
        symbol = str(ValCycl(args["Marq"], 0, 10, args["NumSet"]))
        freqm = str(ValCycl(args["FreqM"], 0, -1, 0))

        sn = str(args["NumSet"])
        descr = []
        descr.append(
            """
@    s0 hidden false
@    s0 type xy
@    s0 symbol size 1.000000
@    s0 symbol pattern 1
@    s0 symbol linestyle 1
@    s0 symbol fill pattern 0
@    s0 symbol linewidth 1.0
@    s0 symbol char 65
@    s0 symbol char font 0
@    s0 line type 1
@    s0 line linewidth 1.0
@    s0 line pattern 1
@    s0 baseline type 0
@    s0 baseline off
@    s0 dropline off
@    s0 fill type 0
@    s0 fill rule 0
@    s0 fill pattern 1
@    s0 avalue off
@    s0 avalue type 2
@    s0 avalue char size 1.000000
@    s0 avalue font 0
@    s0 avalue rot 0
@    s0 avalue format general
@    s0 avalue prec 3
@    s0 avalue prepend ""
@    s0 avalue append ""
@    s0 avalue offset 0.000000 , 0.000000
@    s0 errorbar on
@    s0 errorbar place both
@    s0 errorbar pattern 1
@    s0 errorbar size 1.000000
@    s0 errorbar linewidth 1.0
@    s0 errorbar linestyle 1
@    s0 errorbar riser linewidth 1.0
@    s0 errorbar riser linestyle 1
@    s0 errorbar riser clip off
@    s0 errorbar riser clip length 0.100000

@    s0 comment ""
""".replace(
                " s0 ", " s" + sn + " "
            )
        )
        descr.append("@    s" + sn + " symbol " + symbol)
        descr.append("@    s" + sn + " symbol color " + color)
        descr.append("@    s" + sn + " symbol skip " + freqm)
        descr.append("@    s" + sn + " symbol fill color " + color)
        descr.append("@    s" + sn + " line linestyle " + sty)
        descr.append("@    s" + sn + " line color " + color)
        descr.append("@    s" + sn + " fill color " + color)
        descr.append("@    s" + sn + " avalue color " + color)
        descr.append("@    s" + sn + " errorbar color " + color)
        descr.append("@    s" + sn + ' legend "' + args["Leg"].replace(os.linesep, " ; ") + '"')
        return descr

    # ------------------------------------------------------------------------

    def Trace(self):
        """Méthode pour 'tracer' l'objet Graph dans un fichier.
        Met en page l'entete, la description des courbes et les valeurs selon
        le format et ferme le fichier.
        """
        g = self.Graph
        if self.PILOTE.startswith("INTERACTIF"):
            self.NomFich[0] = "Trace_%s.dat" % time.strftime("%y%m%d%H%M%S", time.localtime())
            self.Fich[0] = open(self.NomFich[0], "w")
        # initialise le graph
        self._FermFich()
        nbsets, x0, x1, y0, y1 = IniGrace(self.NomFich[0])
        NumSetIni = nbsets + 1
        g.SetExtrema(0.05, x0, x1, y0, y1, force=False)
        # si Min/Max incohérents
        if g.Echelle_X == "LOG":
            g.Grille_X = 10
            if g.Min_X < 0.0:
                if g.BBXmin < 0.0:
                    UTMESS("A", "GRAPH0_4")
                g.Min_X = g.MinP_X
        if g.Echelle_Y == "LOG":
            g.Grille_Y = 10
            if g.Min_Y < 0.0:
                if g.BBYmin < 0.0:
                    UTMESS("A", "GRAPH0_5")
                g.Min_Y = g.MinP_Y

        if g.NbCourbe < 1:
            self._FermFich()
            return
        # cohérence des valeurs par défaut
        if g.Grille_X < 0 or g.Grille_Y < 0:
            deltaX = g.Max_X - g.Min_X
            deltaY = g.Max_Y - g.Min_Y
            g.Grille_X = deltaX / 5.0
            g.Grille_Y = deltaY / 5.0
            if deltaX > 4:
                g.Grille_X = int(round(g.Grille_X))
            if deltaY > 4:
                g.Grille_Y = int(round(g.Grille_Y))
            if g.Grille_X == 0.0:
                g.Grille_X = 1.0e-6
            if g.Grille_Y == 0.0:
                g.Grille_Y = 1.0e-6
        # entete
        content = self.Entete()
        content.append("")
        # valeurs
        it = -1
        for i in range(g.NbCourbe):
            dCi = g.Courbe(i)
            for k in range(dCi["NbCol"] - 1):
                it = it + 1
                dCi["NumSet"] = NumSetIni + it
                content.extend(self.DescrCourbe(**dCi))
                content.append("")
        # partie données (.dat)
        it = -1
        for i in range(g.NbCourbe):
            dCi = g.Courbe(i)
            for k in range(dCi["NbCol"] - 1):
                it = it + 1
                content.append("@target g0.s%d" % (NumSetIni + it))
                content.append("@type xy")
                listX, listY = Tri(g.Tri, lx=dCi["Abs"], ly=dCi["Ord"][k])
                for j in range(dCi["NbPts"]):
                    svX = self.DicForm["formR"] % listX[j]
                    svY = self.DicForm["formR"] % listY[j]
                    content.append(
                        self.DicForm["formR"] % listX[j] + " " + self.DicForm["formR"] % listY[j]
                    )
                content.append("&")
        content.append("")

        # Production du fichier postscript, jpeg ou lancement interactif
        pilo = self.PILOTE
        if pilo == "":
            self._OuvrFich()
            self.Fich[0].write("\n".join(content))
            self._FermFich()
        else:
            xmgr = ExecutionParameter().get_option("prog:xmgrace")
            nfwrk = self.NomFich[0] + ".wrk"
            with open(nfwrk, "w") as f:
                f.write("\n".join(content))
            nfhard = self.NomFich[0] + ".hardcopy"
            # nom exact du pilote
            bg = pilo == "INTERACTIF_BG"
            if pilo == "POSTSCRIPT":
                pilo = "PostScript"
            elif pilo.startswith("INTERACTIF"):
                pilo = "X11"
            # ligne de commande
            if pilo == "X11":
                lcmde = "%s %s" % (xmgr, nfwrk)
                if bg:
                    lcmde += " &"
                if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
                    os.environ["DISPLAY"] = ":0.0"
                    UTMESS("I", "GRAPH0_7")
                UTMESS("I", "GRAPH0_8", valk=os.environ["DISPLAY"])
            else:
                if ExecutionParameter().get_option("prog:gracebat") != "gracebat":
                    xmgr = ExecutionParameter().get_option("prog:gracebat")
                lcmde = "%s -hdevice %s -hardcopy -printfile %s %s" % (xmgr, pilo, nfhard, nfwrk)
            # appel xmgrace
            UTMESS("I", "EXECLOGICIEL0_8", valk=lcmde)
            if not os.path.exists(xmgr):
                UTMESS("S", "EXECLOGICIEL0_6", valk=xmgr)
            iret = os.system(lcmde)
            if iret == 0 or os.path.exists(nfhard):
                if pilo not in ("", "X11"):
                    with open(nfhard, "rb") as f:
                        new = f.read()
                    with open(self.NomFich[0], "ab") as f:
                        f.write(new)
            else:
                UTMESS("A", "GRAPH0_9", valk=pilo)
        # menage
        if self.PILOTE.startswith("INTERACTIF"):
            os.remove(self.NomFich[0])
        return


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------


class TraceAgraf(TraceGraph):
    """
    Impression d'un objet Graph au format AGRAF.
    """

    # ------------------------------------------------------------------------

    def Entete(self):
        """Retourne l'entete des directives Agraf"""
        dic_ech = {"LIN": "0", "LOG": "1"}
        g = self.Graph
        entete = []
        entete.append(
            """
ASPECT_GRAPHIQUE:
  En-tete :Departement Analyses Mecaniques et Acoustique
  Aspect :0
  Nombre de vues :1
  Cesure commentaire :40
  MinMax :0
  Fonte Titre :%helvetica-14
  Fonte Axes :%courier-12
  Fonte Autre :%times-12

  DEFAUT_COURBE:
    Couleur (rvb) :     0     0     0

  DEFAUT_COURBE:
    Couleur (rvb) : 65535     0     0

  DEFAUT_COURBE:
    Couleur (rvb) : 11822 35723 22359

  DEFAUT_COURBE:
    Couleur (rvb) :     0     0 65535

  DEFAUT_COURBE:
    Couleur (rvb) : 65535     0 65535

  DEFAUT_COURBE:
    Couleur (rvb) :     0 65535 65535

  DEFAUT_COURBE:
    Couleur (rvb) :     0 65535     0

  DEFAUT_COURBE:
    Couleur (rvb) : 41120 21074 11565

  DEFAUT_COURBE:
    Couleur (rvb) : 65535 42405     0

  DEFAUT_COURBE:
    Couleur (rvb) : 41120  8224 61680

  DEFAUT_COURBE:
    Couleur (rvb) : 65535 65535     0

  DEFAUT_COURBE:
    Couleur (rvb) : 53970 46260 35980

GRAPHIQUE:
"""
        )
        if g.Titre == "":
            g.Titre = "GRAPHIQUE CODE_ASTER"
        entete.append("Titre :" + " ".join(g.Titre) + "\n")
        if g.SousTitre != "":
            entete.append("Commentaire :" + g.SousTitre + "\n")
        entete.append("Frequence Grille X :" + str(int(g.Grille_X)) + "\n")
        entete.append("Frequence Grille Y :" + str(int(g.Grille_Y)) + "\n")
        entete.append("Echelle X :" + dic_ech[g.Echelle_X] + "\n")
        entete.append("Echelle Y :" + dic_ech[g.Echelle_Y] + "\n")
        if g.Legende_X != "":
            entete.append("Legende X :" + g.Legende_X + "\n")
        if g.Legende_Y != "":
            entete.append("Legende Y :" + g.Legende_Y + "\n")
        entete.append("Min X : " + str(g.Min_X) + "\n")
        entete.append("Max X : " + str(g.Max_X) + "\n")
        entete.append("Min Y : " + str(g.Min_Y) + "\n")
        entete.append("Max Y : " + str(g.Max_Y) + "\n")

        return entete

    # ------------------------------------------------------------------------

    def DescrCourbe(self, **args):
        """Retourne la chaine de caractères décrivant les paramètres de la courbe."""
        # valeurs par défaut
        sty = str(ValCycl(args["Sty"], 0, 2, 0))
        color = str(ValCycl(args["Coul"], 0, 12, args["NumSet"]))
        symbol = str(ValCycl(args["Marq"], 0, 12, args["NumSet"]))
        freqm = str(ValCycl(args["FreqM"], 0, -1, 0))

        descr = []
        descr.append("  COURBE:\n")
        descr.append("     Trait :" + sty + "\n")
        descr.append("     Couleur :" + color + "\n")
        descr.append("     Marqueur :" + symbol + "\n")
        descr.append("     Frequence Marqueur :" + freqm + "\n")
        if args["Leg"] != "":
            descr.append("     Legende :" + args["Leg"] + "\n")
        descr.append("     Tri :" + args["Tri"] + "\n")
        descr.append("     Abscisses : [ " + str(args["Bloc"]) + ", " + str(args["ColX"]) + "]\n")
        descr.append("     Ordonnees : [ " + str(args["Bloc"]) + ", " + str(args["ColY"]) + "]\n")
        return descr

    # ------------------------------------------------------------------------

    def Trace(self):
        """Méthode pour 'tracer' l'objet Graph dans un fichier.
        Met en page l'entete, la description des courbes et les valeurs selon
        le format et ferme le fichier.
        """
        self._OuvrFich()
        fdogr = self.Fich[0]
        fdigr = self.Fich[1]
        g = self.Graph
        if g.NbCourbe > 0:
            # cohérence des valeurs par défaut
            if g.Grille_X < 0 or g.Grille_Y < 0:
                g.Grille_X = 0
                g.Grille_Y = 0
            # entete
            for lig in self.Entete():
                fdigr.write(lig)
            # valeurs
            for i in range(g.NbCourbe):
                dCi = g.Courbe(i)
                dCi["NumSet"] = i
                # partie directives (.digr)
                for k in range(dCi["NbCol"] - 1):
                    dCi["Bloc"] = i + 1
                    dCi["ColX"] = 1
                    dCi["ColY"] = k + 2
                    for lig in self.DescrCourbe(**dCi):
                        fdigr.write(lig)
                # partie données (.dogr)
                if dCi["Leg"] != "":
                    leg = dCi["Leg"]
                else:
                    leg = "COURBE_" + str(i)
                fdogr.write("#NOM DE LA FONCTION: " + leg + "\n")
                for j in range(dCi["NbPts"]):
                    for k in range(dCi["NbCol"]):
                        sv = self.DicForm["formR"] % g.Valeurs[i][k][j]
                        fdogr.write(" " + sv)
                    fdogr.write("\n")
                fdogr.write("\n")
        self._FermFich()


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------


class TraceMatplotlib(TraceGraph):
    def Entete(self):
        entete = " "
        return entete

    def DescrCourbe(self, **args):
        desc = " "
        return desc

    def Trace(self):
        fichier = self.NomFich
        g = self.Graph

        figsize, dpi = [160.0 / 2.54 / 2 * 1.35, 100.0 / 2.54 * 1.35], 160
        fig = plt.figure(1, figsize=figsize, dpi=dpi)
        fig.clear()
        ax1 = fig.add_subplot(1, 1, 1)
        liss_nappe = {}
        for i in range(g.NbCourbe):
            dCi = g.Courbe(i)
            Leg = dCi["Leg"]
            if Leg.startswith("NAPPE_LISSEE"):
                stycol = "-ok"
                AMOR = Leg.split("=")[-1]
                liss_nappe[float(AMOR)] = [dCi["Abs"], dCi["Ord"][0]]
            else:
                stycol = "--b"
            ax1.plot(dCi["Abs"], dCi["Ord"][0], stycol, linewidth=2)

        if g.Echelle_X == "LOG":
            ax1.set_xscale("log")
        if g.Echelle_Y == "LOG":
            ax1.set_yscale("log")

        plt.title(g.Titre[0], fontsize=32)
        if g.Min_X is not None and g.Max_X is not None:
            plt.xlim([g.Min_X, g.Max_X])
        if g.Min_Y is not None and g.Max_Y is not None:
            plt.ylim([g.Min_Y, g.Max_Y])
        plt.xlabel(g.Legende_X, fontsize=32)
        plt.ylabel(g.Legende_Y, fontsize=32)
        ax = plt.gca()
        ax.xaxis.set_label_coords(0.5, 0.01)
        for lab in ax.xaxis.get_ticklabels():
            lab.set_fontsize(32)
        for lab in ax.yaxis.get_ticklabels():
            lab.set_fontsize(32)
        plt.grid(True, which="both")

        if len(liss_nappe) > 0:
            tablefig = []
            listAmor = list(liss_nappe.keys())
            listAmor.sort()
            listFreq = liss_nappe[listAmor[0]][0]
            for i in range(len(listFreq)):
                row = []
                row.append("%.2f" % listFreq[i])
                for amor in listAmor:
                    row.append("%.3f" % liss_nappe[amor][1][i])
                tablefig.append(row)
            labelc = ["Freq \n[Hz]"]
            for amor in listAmor:
                labelc.append("Damp \n %.1f" % (amor * 100) + "%")

            agg = plt.table(
                cellText=tablefig,
                colWidths=[0.06] * 9,
                colColours=["white"] * 9,
                cellLoc="center",
                colLabels=labelc,
                colLoc="center",
                loc="top",
                alpha=1.0,
                zorder=10,
            )
            table_cells = agg.get_children()
            for cell in table_cells:
                txt = cell.get_text().get_text()
                if "Freq" in txt or "Damp" in txt:
                    cell.set_height(0.03)
                else:
                    cell.set_height(0.018)

            agg.set_fontsize(28)
            agg.scale(1.2, 1.2)

        titre = g.Titre[0]
        sous_titre = g.SousTitre
        if "," in sous_titre:
            sssplit = sous_titre.split(",")
            sous_titre_1 = sssplit[0]
            sous_titre_2 = sssplit[1]
        else:
            sous_titre_1 = sous_titre
            sous_titre_2 = ""
        legend = plt.table(
            cellText=[[titre], [sous_titre_1], [sous_titre_2]],
            colWidths=[0.3] * 1,
            cellLoc="center",
            rowLoc="center",
            loc="upper left",
            zorder=10,
            fontsize=32,
        )
        table_cells = legend.get_children()
        table_cells[2].set_height(0.05)
        table_cells[1].set_height(0.05)
        table_cells[0].set_height(0.1)
        legend.set_fontsize(40)
        legend.scale(1.3, 1.3)

        with disable_fpe():
            plt.savefig(fichier[0], format="png", bbox_inches="tight")


def ValCycl(val, vmin, vmax, vdef):
    """
    Retourne une valeur entre vmin et vmax (bornes incluses) :
       - si val<vmin, on utilise val=vdef,
       - si val>vmax, on cycle tel que val=vmax+1 retourne vmin, etc.
       - si vmax<vmin, il n'y a pas de max
    """
    if val is None:
        val = vmin - 1
    if val < vmin:
        val = vdef
    if vmax < vmin:
        return val
    else:
        return ((val - vmin) % (vmax + 1 - vmin)) + vmin


# ------------------------------------------------------------------------


def Tri(tri, lx, ly):
    """Retourne les listes triées selon la valeur de tri ('X','Y','XY','YX')."""
    dNumCol = {"X": 0, "Y": 1}
    tab = np.array((lx, ly))
    tab = np.transpose(tab)
    li = list(range(len(tri)))
    li.reverse()
    for i in li:
        if tri[-i] in list(dNumCol.keys()):
            icol = dNumCol[tri[-i]]
            tab = np.take(tab, np.argsort(tab[:, icol]))
    return [tab[:, 0].tolist(), tab[:, 1].tolist()]


# ------------------------------------------------------------------------


def AjoutParaCourbe(dCourbe, args):
    """Ajoute les arguments fournis par l'utilisateur (args) dans le dictionnaire
    décrivant la courbe (dCourbe).
    """
    # correspondance : mot-clé Aster / clé du dico de l'objet Graph
    keys = {
        "LEGENDE": "Leg",
        "STYLE": "Sty",
        "COULEUR": "Coul",
        "MARQUEUR": "Marq",
        "FREQ_MARQUEUR": "FreqM",
        "TRI": "Tri",
    }
    for mc, key in list(keys.items()):
        if mc in args:
            dCourbe[key] = args[mc]


# ------------------------------------------------------------------------


def IniGrace(fich):
    """Retourne le numéro de la dernière courbe d'un fichier xmgrace (sinon 0)."""
    ns = 0
    x0 = None
    x1 = None
    y0 = None
    y1 = None
    if os.path.exists(fich) and os.stat(fich).st_size != 0:
        assert not is_binary(fich), "Can not append text to a binary file"
        shutil.copy(fich, fich + ".prev")
        fpre = open(fich + ".prev", "r")
        fnew = open(fich, "w")
        for line in fpre:
            ikeep = True
            mat = re.search(r"@target g[0-9]+\.s([0-9]+)", line)
            if mat is not None and int(mat.group(1)) > ns:
                ns = int(mat.group(1))
            mat = re.search(r"@[ ]+world[ ]+xmin[ ]+([\-\+\.0-9eEdD]+)", line)
            if mat is not None:
                try:
                    x0 = float(mat.group(1))
                    ikeep = False
                except ValueError:
                    pass
            mat = re.search(r"@[ ]+world[ ]+xmax[ ]+([\-\+\.0-9eEdD]+)", line)
            if mat is not None:
                try:
                    x1 = float(mat.group(1))
                    ikeep = False
                except ValueError:
                    pass
            mat = re.search(r"@[ ]+world[ ]+ymin[ ]+([\-\+\.0-9eEdD]+)", line)
            if mat is not None:
                try:
                    y0 = float(mat.group(1))
                    ikeep = False
                except ValueError:
                    pass
            mat = re.search(r"@[ ]+world[ ]+ymax[ ]+([\-\+\.0-9eEdD]+)", line)
            if mat is not None:
                try:
                    y1 = float(mat.group(1))
                    ikeep = False
                except ValueError:
                    pass
            if ikeep:
                fnew.write(line)
        fpre.close()
        fnew.close()
        if None not in (x0, x1, y0, y1):
            UTMESS("I", "GRAPH0_10", valk=fich, vali=ns, valr=(x0, x1, y0, y1))
        else:
            UTMESS("A", "GRAPH0_11", valk=fich)
    return ns, x0, x1, y0, y1


def is_binary(fname):
    """Tell if a file appears to be binary.
    All non UTF-8 files are considered binary.
    """
    with open(fname, "rb") as fobj:
        try:
            fobj.read().decode("utf-8")
        except UnicodeDecodeError:
            return True
    return False
