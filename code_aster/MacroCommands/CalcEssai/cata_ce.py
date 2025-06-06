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

# person_in_charge: albert.alarcon at edf.fr

# On definit dans ce module python des classes associees a des objets aster,
# permettant leur manipulation en python (pour le calcul numerique).
# - Classe Resultat : regroupe tous les objets sd_resultat. Elle est devouee
# mode_meca, grace a des methodes telles get_modes() qui recupere les
# caracteristiques modales
# - CaraElem
# - ChampMateriau
# - InterSpectre : regroupe les matrices python et l'objet Aster. Les calculs inverses
# d'effort sont faits sur des matrices python. Les concepts Aster sont
# crees par defaut.


import aster
import numpy

from ...Cata.DataStructure import (
    cara_elem,
    cham_mater,
    dyna_harmo,
    interspectre,
    maillage_sdaster,
    matr_asse_depl_r,
    matr_asse_gene_r,
    mode_meca,
    modele_sdaster,
    nume_ddl_sdaster,
    table_fonction,
    table_sdaster,
)
from ...Cata.Syntax import _F
from ...CodeCommands import CREA_CHAMP, DEFI_FONCTION, DEFI_INTE_SPEC, RECU_FONCTION
from ...Messages import UTMESS
from ...Objects import DataStructure, FullResult, HarmoGeneralizedResult, TransientGeneralizedResult


class Resultat:

    """Classe chapeau a toutes les sd_resultat traitees dans CALC_ESSAI"""

    # TODO : mettre ici les attributs et procedures communs a tous les resus

    def __init__(self, objects, nom, obj_ast, mess):
        self.objects = objects
        self.nom = nom.strip()
        self.obj = obj_ast
        self.mess = mess
        self.cara_mod = None

    def get_modele(self):
        """Recherche le modele associe au resultat. Routine generique pour les dyna_harmo et mode_meca"""
        if not self.modele:
            model = self.obj.getModel()
            if not model:
                model = self.obj.getDOFNumbering().getModel()
            if model:
                modele_name = model.getName()
                if modele_name[0:1] != "#":
                    self.modele_name = modele_name.rstrip()
                    self.modele = self.objects.modeles[self.modele_name]
                    return
        # Si cela ne marche pas, on passe par le maillage
        if not self.modele:
            for _mod in list(self.objects.modeles.values()):
                if not _mod.maya_name:
                    _mod.get_maillage()

        if not self.modele:
            self.mess.disp_mess("On ne trouve pas le modele associe au resultat " + self.nom)
            UTMESS("A", "CALCESSAI0_8", valk=(self.nom))


# raise IndexError


class ModeMeca(Resultat):

    """!Gestion des sd_resultat d'aster

    Permet de conserver une reference sur un objet aster sd_resultat
    et permet aussi de recuperer facilement les concepts aster associes
    tels le modele, le maillage, la numerotation, matrices de masse et raideur
    """

    def __init__(
        self,
        objects=None,  # macro CalcEssaiObjects parente
        nom=None,  # nom du concept aster
        obj_ast=None,  # concept Aster
        mess=None,  # fenetre de messages
    ):
        """Constructeur"""
        Resultat.__init__(self, objects, nom, obj_ast, mess)
        self.modele_name = ""
        self.modele = None
        self.maya_name = ""
        self.maya = None
        self.nume_name = ""
        self.nume = None
        self.mass_name = ""
        self.mass = None
        self.kass_name = ""
        self.kass = None
        self.cass_name = ""
        self.cass = None
        self.nom_cham = None

        try:
            self.get_nom_cham()
        except:
            pass

        # self.show_linked_concepts()

    def get_nume(self):
        """Recuperation de la numerotation et du nume_ddl"""
        if isinstance(self.obj, (FullResult, HarmoGeneralizedResult, TransientGeneralizedResult)):
            self.nume_name = self.obj.getDOFNumbering().getName()
            if self.nume_name:
                self.nume = self.objects.nume_ddl[self.nume_name]

    def get_maillage(self):
        """Recherche du maillage associe au resultat"""
        if not self.nume_name:
            self.get_nume()
        if self.nume_name:  # methode 1
            self.maya_name = self.nume.getMesh().getName()
            self.maya = self.objects.maillages[self.maya_name]
        else:  # methode 2
            self.maya_name = self.obj.getMesh().getName()
            self.maya = self.objects.maillages[self.maya_name]

    def show_linked_concepts(self):
        """!Affichage du concept resultats et des concepts lies"""
        self.mess.disp_mess((self.nom + " : "))
        self.mess.disp_mess((".modele" + self.modele_name))
        self.mess.disp_mess((".maillage" + self.maya_name))
        self.mess.disp_mess((".nume" + self.nume_name))
        self.mess.disp_mess((".mass" + self.mass_name))
        self.mess.disp_mess((".kass" + self.kass_name))
        self.mess.disp_mess((" "))

    def get_modes_data(self):
        """retourne pour tout type de mode_meca les caras
        dynamiques (frequences, amor...) et les nom_cmp pour
        les modes statiques quand ils le sont"""
        cara_mod = {
            "NUME_ORDRE": [],
            "FREQ": [],
            "AMOR_REDUIT": [],
            "AMOR_GENE": [],
            "RIGI_GENE": [],
            "MASS_GENE": [],
            "NUME_MODE": [],  # VOIR SI ON A VRAIMENT BESOIN DE CELUI-LA. SIMPLIFIER AU MAX
            "NOEUD_CMP": [],
        }
        resu_stat = self.obj.LIST_PARA()["NOEUD_CMP"]
        nume_ordr = self.obj.LIST_PARA()["NUME_ORDRE"]
        cara_mod["NUME_MODE"] = self.obj.LIST_PARA()["NUME_MODE"]
        cara_mod["NUME_ORDRE"] = self.obj.LIST_PARA()["NUME_ORDRE"]

        # Rangement des donnees modales selon le type : statique ou dynamique
        # Si une donnee est incompatible avec le type (exemple : la frequence d'un mode statique), on remplit par None
        # TODO : depuis issue15113, le "if/else" n'est peut-etre plus
        # necessaire.
        for ind_ordr in range(len(nume_ordr)):
            liste = ["FREQ", "AMOR_REDUIT", "AMOR_GENE", "RIGI_GENE", "MASS_GENE"]
            if resu_stat[ind_ordr] is not None:  # mode statique
                for ind_list in liste:
                    cara_mod[ind_list].append(None)
                cara_mod["NOEUD_CMP"].append(resu_stat[ind_ordr])
            else:  # mode dynamique
                for ind_list in liste:
                    cara_mod[ind_list].append(self.obj.LIST_PARA()[ind_list][ind_ordr])
                cara_mod["NOEUD_CMP"].append(None)

        self.cara_mod = cara_mod
        return cara_mod

    def extr_matr(self):
        """Extrait les champs de deformees contenus dans le resultat"""

        self.nume_phy = nume_ddl_phy(self)

        if not self.cara_mod:
            self.cara_mod = self.get_modes_data()
        nb_mod = len(self.cara_mod["NUME_MODE"])  # nb de modes

        matrice = []
        for ind_mod in range(1, nb_mod + 1):
            champ = crea_champ(self.obj, ind_mod)
            matrice.append(champ)

        matrice = numpy.transpose(numpy.array(matrice))

        return matrice

    def get_nom_cham(self):
        """Recherche le type de champ rempli dans la sd ('ACCE', 'DEPL'...)"""
        self.nom_cham = self.obj.getFieldsOnNodesRealNames()[0]
        # cette methode sert a PROJ_MESUR_MODAL dans MACRO_EXPANS : on ne garde
        # donc qu'un seul nom symbolique. Par ordre decroissant de priorite :
        # 'DEPL', 'VITE', 'ACCE', etc...


# -------------------------------------------------------------------------


class DynaHarmo(Resultat):

    """pour les resultats de type dyna_harmo"""

    def __init__(
        self,
        objects=None,  # macro CalcEssaiObjects parente
        nom=None,  # nom du concept aster
        obj_ast=None,  # concept Aster
        mess=None,  # fenetre de messages
    ):
        Resultat.__init__(self, objects, nom, obj_ast, mess)
        """Constructeur"""
        self.cara_mod = None
        self.modele_name = ""
        self.modele = None
        self.maya_name = ""
        self.maya = None
        self.nume_name = ""
        self.nume = None
        self.mass_name = ""
        self.mass = None
        self.kass_name = ""
        self.kass = None
        self.mess = mess
        self.nom_cham = None
        self.nume_ddl = None

        try:
            self.get_nom_cham()
        except AttributeError:
            pass

    def extr_freq(self):
        vari_acces = list(
            zip(self.obj.LIST_VARI_ACCES()["NUME_ORDRE"], self.obj.LIST_VARI_ACCES()["FREQ"])
        )
        return vari_acces

    def get_maillage(self):
        self.maya_name = self.obj.getMesh().getName()
        self.maya = self.objects.maillages[self.maya_name]

    def get_nume(self):
        """Recherche d'un nume_ddl. Le pb est qu'il n'est pas necessaire
        d'avoir un nume_ddl pour faire un dyna_harmo. On cherche donc
        le maillage associe au resultat et on regarde quel nume
        possede ce meme maillage"""

        self.get_maillage()
        self.nume_ddl_name = self.obj.getDOFNumbering().getName()
        self.nume_ddl = self.objects.nume_ddl[nume]
        return self.nume_ddl

    def get_nom_cham(self):
        """Recherche le type de champ rempli dans la sd ('ACCE', 'DEPL'...)"""
        self.nom_cham = self.obj.getFieldsOnNodesRealNames()[0]
        # cette methode sert a PROJ_MESU_MODAL dans MACRO_EXPANS : on ne garde
        # donc qu'un seul nom symbolique. Par ordre decroissant de priorite :
        # 'DEPL', 'VITE', 'ACCE', etc...

    def get_modes_data(self):
        """retourne pour tout type de mode_meca les caras
        dynamiques (frequences, amor...) et les nom_cmp pour
        les modes statiques quand ils le sont"""
        cara_mod = {"NUME_ORDRE": [], "FREQ": [], "NOEUD_CMP": []}

        cara_mod["NUME_ORDRE"] = self.obj.LIST_PARA()["NUME_ORDRE"]
        cara_mod["FREQ"] = self.obj.LIST_PARA()["FREQ"]
        cara_mod["NOEUD_CMP"] = [None] * len(cara_mod["FREQ"])

        self.cara_mod = cara_mod
        return cara_mod


# -------------------------------------------------------------------------


class InterSpectre:

    """!Gestion des concepts de type interspectre

    Regroupe les concepts aster de type interspectre :
    - Extrait les sd_interspectre sous forme d'une matrice python interspectrale
    - Cree un lien entre les numerotations des ddl de cette matrice
      avec les numerotations d'un modele EF. ex : la ligne/colonne 3 de la
      matrice interspectrale correspond au noeud 1DZ du modele
    - Cree une sd_interspectre a partir d'une matrice python
        On peut creer une table avec le nom, la table format aster (obj_ast) ou la matrice
        format python (mat)
    """

    def __init__(self, nom=None, obj_ast=None, mat=None, frequences=[], mess=None, var_opt=None):
        self.nom = nom.strip()  # nom aster de la sd
        self.obj = obj_ast  # objet aster (a remplir ou a fabriquer)
        self.matr_inte_spec = mat  # matrice inter-spectrale format python
        self.intsp = 0  # vaut 1 si la table est un intsp, 0 sinon
        self.var_opt(var_opt)
        # definit self.opt : vaut 1 si les modeles sont definis
        self.f = frequences
        self.resu_name = ""
        self.resu = None
        self.maya_name = ""
        self.maya = None
        self.nume_name = ""
        self.nume_phy = None
        self.nume_gene = None
        self.mass_name = ""
        self.mass = None
        self.kass_name = ""
        self.kass = None
        self.options_meth = ["Efforts discrets localises", "Efforts et moments discrets"]
        self.mess = mess

        # -- genere une erreur si pas OK -> on passe par Tempo

        if len(self.f) == 0:
            self.extr_freq()
        self.intsp = 1

    def make_inte_spec(self, titre, paras_out):
        """
        Fabrique un objet aster de type sd_interspectre :
        """
        dim = self.matr_inte_spec.shape[1]
        nb_freq = len(self.f)
        l_fonc = []
        nume_ordr = []
        for i in range(dim):
            for j in range(i, dim):
                fonc = []
                for ind_freq in range(nb_freq):
                    fonc.append(self.f[ind_freq])
                    fonc.append(self.matr_inte_spec[ind_freq, i, j].real)
                    fonc.append(self.matr_inte_spec[ind_freq, i, j].imag)
                nume_ordr.append([int(i), int(j)])

                _fonc = DEFI_FONCTION(
                    NOM_PARA="FREQ", NOM_RESU="DSP", INTERPOL="LIN", INFO=1, VALE_C=fonc
                )
                l_fonc.append(_fonc)

        nume_ordr = numpy.array(nume_ordr)
        nume_i = nume_ordr[:, 0]
        nume_j = nume_ordr[:, 1]
        mcfact = []
        if isinstance(self.resu, Resultat):
            # Si on associe l'inter-spectre a un resultat,
            # on range les fonctions par rapport aux noeuds et composantes
            if not self.nume_phy:
                self.nume_phy = nume_ddl_phy(self.resu)
            ddl = self.nume_phy
            noeu_i = []
            noeu_j = []
            cmp_i = []
            cmp_j = []
            for ind in range(len(nume_i)):
                li = ddl[nume_i[ind]].split("_")
                lj = ddl[nume_j[ind]].split("_")
                noeu_i.append(li[0])
                noeu_j.append(lj[0])
                cmp_i.append(li[1])
                cmp_j.append(lj[1])
            for k in range(len(nume_i)):
                mcfact.append(
                    _F(
                        NOEUD_I=noeu_i[k],
                        NOEUD_J=noeu_j[k],
                        NOM_CMP_I=cmp_i[k],
                        NOM_CMP_J=cmp_j[k],
                        FONCTION=l_fonc[k],
                    )
                )
        else:
            for k in range(len(nume_i)):
                mcfact.append(
                    _F(NUME_ORDRE_I=1 + nume_i[k], NUME_ORDRE_J=1 + nume_j[k], FONCTION=l_fonc[k])
                )

        __inte_out = CreaTable(mcfact, titre, paras_out, self.mess)

    def var_opt(self, opt):
        if opt == "Efforts discrets localises":
            self.opt = 0
        elif opt == "Efforts et moments discrets":
            self.opt = 1
        else:
            self.opt = 0

    def set_model(self, resu):
        """Lie l'inter-spectre au concept mode_meca OBS. Permet de lier les
        lignes et colonnes de l'inter-spectre aux DDL des deformees modales
        et de tout ranger dans le bon ordre. Si l'inter-spectre est defini
        avec des numeros d'ordre, alors on suppose qu'ils sont ranges dans le bon
        ordre.
        """
        self.resu = resu

    def extr_nume_ordr(self):
        """! Extraction des numeros d'ordre de l'inter-spectre"""
        coupl_ddl = []

        # Cas ou l'inter-spectre est defini par ses noeuds et composantes
        noeudi = self.obj.getLineNodes()
        noeudj = self.obj.getColumnNodes()
        cmpi = self.obj.getLineComponents()
        cmpj = self.obj.getColumnComponents()

        # l'inter-spectre n'est defini qu'avec des numeros d'ordre independants
        # du modele
        numi = self.obj.getLineIndexes()
        numj = self.obj.getColumnIndexes()

        if noeudi:
            self.isnume = 1
            for ind in range(len(noeudi)):
                coupl_ddl.append(
                    (
                        noeudi[ind].split()[0] + "_" + cmpi[ind].split()[0],
                        noeudj[ind].split()[0] + "_" + cmpj[ind].split()[0],
                    )
                )
        elif numi:
            self.isnume = 0
            for ind in range(len(numi)):
                coupl_ddl.append((numi[ind], numj[ind]))

        return coupl_ddl

    def extr_inte_spec(self, resu=None, intersp=1):
        """!Extraction d'une matrice inter-spectrale a partir d'une sd_interspectre
        Verification de la coherence entre les ddl de l'inter-spectre et du concept resultat
         - Si resi is not None, on peut associer les numeros d'ordre de l'IS a des DDL du resu,
         - Si intersp = 1, on a un vrai inter-spectre. Si = 0, alors c'est une FRF ou une coherence
           (exemple, CALC_SEPC calcule des FRF mais rend un concept de type inter-spectre).
           Dans ce cas, on ne calcule pas la partie symetrique de la matrice."""
        self.mess.disp_mess("Extraction de l'inter-spectre " + self.nom)
        self.mess.disp_mess(" ")

        coupl_ddl = self.extr_nume_ordr()
        nb_fonc = len(coupl_ddl)
        nb_freq = len(self.f)

        nume_ordr_l = list(set([kk[0] for kk in coupl_ddl]))
        nb_l = len(nume_ordr_l)
        nume_ordr_c = list(set([kk[1] for kk in coupl_ddl]))
        nb_c = len(nume_ordr_c)

        if intersp:
            if nb_l != nb_c:
                self.mess.disp_mess("Inter-spectre non valide")
                return

        # Initialisation
        self.matr_inte_spec = numpy.zeros((nb_freq, nb_l, nb_c), complex)

        # Y a-t-il une numerotation associee a l'inter-spectre. Si oui
        if resu:
            self.set_model(resu)
            self.nume_phy = nume_ddl_phy(resu)
            nb_mes = len(self.nume_phy)
            # verification de la coherence entre la taille de l'inter-spectre et du DDL du resu
            # TODO : retirer la verif ici et la mettre ailleurs
            if nb_mes * (nb_mes + 1) // 2 != nb_fonc:
                nb_mes_intsp = 0.5 * (-1 + numpy.sqrt(1 + 8 * nb_fonc))
                self.mess.disp_mess(" Nombre de mesures de CPhi : " + str(int(nb_mes)))
                self.mess.disp_mess(
                    " Nombre de mesures de l'inter-spectre : " + str(int(nb_mes_intsp))
                )
                self.mess.disp_mess(" ")
                raise TypeError

        for nume_i, nume_j in coupl_ddl:
            nume_i2 = nume_i
            if nume_i[0] == "N":
                nume_i2 = nume_i[1:]
            nume_j2 = nume_j
            if nume_j[0] == "N":
                nume_j2 = nume_j[1:]
            if not self.nume_phy or not self.isnume:
                # rangement alpha-numerique des donnees de l'inter-spectre
                ind_l = nume_ordr_l.index(nume_i2)
                ind_c = nume_ordr_c.index(nume_j2)
                __fonc = RECU_FONCTION(INTE_SPEC=self.obj, NUME_ORDRE_I=nume_i, NUME_ORDRE_J=nume_j)
            elif self.nume_phy:
                # rangement selon l'ordre des DDL donne par self.nume
                ind_l = self.nume_phy.index(nume_i2)
                ind_c = self.nume_phy.index(nume_j2)
                if nume_i == nume_j:
                    __fonc = RECU_FONCTION(
                        INTE_SPEC=self.obj,
                        NOEUD_I=nume_i.split("_")[0],
                        NOM_CMP_I=nume_i.split("_")[1],
                    )
                else:
                    __fonc = RECU_FONCTION(
                        INTE_SPEC=self.obj,
                        NOEUD_I=nume_i.split("_")[0],
                        NOM_CMP_I=nume_i.split("_")[1],
                        NOEUD_J=nume_j.split("_")[0],
                        NOM_CMP_J=nume_j.split("_")[1],
                    )
            else:
                self.mess.disp_mess("Erreur dans l'extraction de l'inter-spectre : cas non-traite")
            fonc_py = __fonc.convert("complex")
            ordo = numpy.array(fonc_py.vale_y)
            if ind_l != ind_c:
                self.matr_inte_spec[:, ind_l, ind_c] = ordo
                self.matr_inte_spec[:, ind_c, ind_l] = numpy.conjugate(ordo)
            else:
                self.matr_inte_spec[:, ind_l, ind_c] = ordo

    def extr_freq(self):
        """Extraction des frequences d'etude dans la tabl_intsp qui contient
        les inter-spectres mesures"""
        self.f = self.obj.getNumberOfFrequencies()
        self.intsp = 1


# -------------------------------------------------------------------------


class Tempo:

    """!Gestion des concepts de type table_fonction contenant des temporels

    Regroupe les concepts aster de type table_fonction :
    - Differencie les tabl_fonc des autres tables tabl_sdaster,
    - Extrait les tabl_fonc sous forme d'une matrice python de temporels
    - Cree un lien entre les numerotations des ddl de cette matrice
      avec les numerotations d'un modele EF. ex : la ligne/colonne 3 de la
      matrice interspectrale correspond au noeud 1DZ du modele
    """

    def __init__(
        self,
        nom=None,
        obj_ast=None,
        tempo=None,
        nume_ordr=None,
        nume_mes=None,
        temps=[],
        mess=None,
        var_opt=None,
    ):
        self.nom = nom.strip()  # nom aster de la sd
        self.obj = obj_ast  # objet aster (a remplir ou a fabriquer)
        self.list_tempo = tempo  # liste de liste de liste format python contenant les temporels
        self.nume_ordr = (
            nume_ordr  # liste contenant les numeros d'ordre auxquels sont rataches les tempo
        )
        self.nume_mes = (
            nume_mes  # liste contenant les numeros de mesures pour chaque numero d'ordre
        )
        self.tempo = 0  # vaut 1 si la table est un tempo, 0 sinon
        self.var_opt(var_opt)
        # definit self.opt : vaut 1 si les modeles sont definis
        self.t = temps
        self.resu_name = ""
        self.resu = None
        self.maya_name = ""
        self.maya = None
        self.nume_name = ""
        self.nume_phy = None
        self.nume_gene = None
        self.mass_name = ""
        self.mass = None
        self.kass_name = ""
        self.kass = None
        self.options_meth = ["Efforts discrets localises", "Efforts et moments discrets"]
        self.mess = mess

        try:
            if len(self.t) == 0:
                self.extr_temps()
            self.ech_t = 1  # self.intsp = 1

        except KeyError:
            # Cas ou la table_sdaster n'est pas un Tempo
            pass  # TODO : faire en sorte que cette table ne soit pas visible

    def var_opt(self, opt):
        if opt == "Efforts discrets localises":
            self.opt = 0
        elif opt == "Efforts et moments discrets":
            self.opt = 1
        else:
            self.opt = 0

    def set_model(self, resu):
        """Lie le temporel au concept mode_meca OBS. Permet de lier les
        catalogues de temporels aux DDL des deformees modales
        et de tout ranger dans le bon ordre. Si les temporels sont definis
        avec des numeros d'ordre, alors on suppose qu'ils sont ranges dans le bon
        ordre.
        """
        self.resu = resu

    def extr_temps(self):
        """Extraction des instants d'etude dans la Tempo qui contient
        les temporels mesures"""
        tabl_py = self.obj.EXTR_TABLE()
        toto = tabl_py["FONCTION"]
        nom_fonc = toto.values()["FONCTION"][0]
        __FONC = RECU_FONCTION(
            TABLE=self.obj,
            NOM_PARA_TABL="FONCTION",
            FILTRE=_F(NOM_PARA="FONCTION", VALE_K=nom_fonc),
        )
        temps = __FONC.Absc()
        self.t = temps
        self.tempo = 1


# -------------------------------------------------------------------------


class Modele:

    """!Gestion des concepts de type modele_sdaster
    Notamment une routine qui permet de fabriquer un nume_ddl
    a partir d'un modele pour les rojtines de type PROJ_CHAMP
    """

    def __init__(self, objects=None, nom=None, obj_ast=None, mess=None, nume_ddl=None):
        self.objects = objects  # les concepts existants dans le jdc
        self.nom = nom.strip()  # nom aster de la sd
        self.obj = obj_ast  # objet aster
        self.nume_ddl = nume_ddl  # Nom d'un nume_ddl associe au modele
        self.mess = mess  # fenetre de messages
        self.maya = None  # concept maillage associe
        self.maya_name = ""  # nom du maillage
        self.kass = None  # matrice de raideur
        self.kass_name = None  # nom de la matrice
        self.mass = None  # matrice de masse
        self.mass_name = None  # nom de ma matrice

    def get_maillage(self):
        self.maya = self.obj.getMesh()
        if self.maya:
            self.maya_name = self.maya.getName()

    def get_nume(self):
        """Recherche des nume_ddl qui dependent de ce modele
        NB : attention, si plusieurs nume_ddl en dependent, seul le premier
        deviendra self.nume_ddl. C'est utile pour les modeles experimentaux
        pour lesquels le ddl est pipo.
        """
        if self.nume_ddl is None:
            for nume_name, nume in list(self.objects.nume_ddl.items()):
                model = nume.getModel()
                if not model:
                    pass
                elif model.getName() == self.nom:
                    self.nume_ddl = nume
                    self.nume_ddl_name = nume_name
                    return
            # TODO : creation automatique d'un nume_ddl pour les resu exp
            # avec des caras bidons.

    def get_matrices(self):
        """Recherche des nume_ddl qui dependent de ce modele
        NB : attention, si plusieurs nume_ddl en dependent, seul le premier
        deviendra self.nume_ddl. C'est utile pour les modeles experimentaux
        pour lesquels le ddl est pipo.
        """
        if self.mass is None or self.kass is None:
            for matr_name, matr in list(self.objects.matrices.items()):
                nom_modele = matr.getModel().getName()
                if nom_modele == self.nom.strip():
                    if matr.getCalculOption() == "RIGI_MECA":
                        self.kass = matr
                        self.kass_name = matr_name
                    if matr.getCalculOption() == "MASS_MECA":
                        self.mass = matr
                        self.mass_name = matr_name

        if self.mass is None or self.kass is None:
            self.mess.disp_mess("On ne trouve pas de matrices associees au modele" + self.nom)
            self.mess.disp_mess("Certains calculs ne seront pas realisables (MAC_MODE)")
        return


# -------------------------------------------------------------------------
class CalcEssaiObjects:

    """!Classe qui recupere les objets pouvant etre utilises par
    CALC_ESSAI dans le catalogue aster"""

    def __init__(self, mess):
        """Constructeur

        Arguments:
            mess (MessageBox): Objet chargé des impressions.
        """
        self.mess = mess
        self.groupno = {}
        self.maillage_modeles = {}
        self.groupno_maillage = {}
        self._init()

    def _init(self):
        self.modeles = {}
        self.mode_meca = {}
        self.dyna_harmo = {}
        self.maillages = {}
        self.matrices = {}
        self.cara_elem = {}
        self.cham_mater = {}
        self.inter_spec = {}
        self.list_tempo = {}
        self.nume_ddl = {}
        self.resultats = {}

    def recup_objects(self, context):
        """Constructeur

        Arguments:
            args (dict): Mots-clés.
        """
        self._init()
        for username, obj in context.items():
            if not isinstance(obj, DataStructure):
                continue
            name = obj.getName()
            self.update(name, obj)

    def update(self, name, obj):
        """Ajout d'un nouvel objet dans self"""
        if obj.getType() == "MODELE_SDASTER":
            self.modeles[name] = Modele(self, name, obj, self.mess)
            self.modeles[name].get_maillage()
            self.modeles[name].get_nume()
        elif obj.getType() == "MAILLAGE_SDASTER":
            self.maillages[name] = obj
        elif obj.getType() in ("MATR_ASSE_DEPL_R", "MATR_ASSE_GENE_R"):
            self.matrices[name] = obj
        elif obj.getType() == "CARA_ELEM":
            self.cara_elem[name] = obj
        elif obj.getType() == "CHAM_MATER":
            self.cham_mater[name] = obj
        elif obj.getType() == "NUME_DDL_SDASTER":
            self.nume_ddl[name] = obj
        elif obj.getType() == "INTERSPECTRE":
            self.inter_spec[name] = InterSpectre(nom=name, obj_ast=obj, mess=self.mess)
        elif obj.getType() == "MODE_MECA":
            self.mode_meca[name] = ModeMeca(self, name, obj, self.mess)
            self.mode_meca[name].get_modele()
            self.mode_meca[name].get_nume()
            self.mode_meca[name].get_maillage()
        elif obj.getType() == "DYNA_HARMO":
            self.dyna_harmo[name] = DynaHarmo(self, name, obj, self.mess)
            self.dyna_harmo[name].get_modele()
            self.dyna_harmo[name].get_nume()
            self.dyna_harmo[name].get_maillage()
        elif obj.getType() == "TABLE_SDASTER":
            # test du type de table :
            typ_table = self.test_table(obj)
            if typ_table == "tempo":
                self.list_tempo[name] = Tempo(nom=name, obj_ast=obj, mess=self.mess)
        else:
            return

        # dict ou on met toutes les sd resu
        self.resultats = {}
        self.resultats.update(self.mode_meca)
        self.resultats.update(self.dyna_harmo)

    def debug(self):
        self.mess.disp_mess(("Modeles" + self.modeles))
        self.mess.disp_mess(("Maillages" + self.maillages))
        self.mess.disp_mess(("Matrices" + self.matrices))
        self.mess.disp_mess(("Resultats"))
        self.mess.disp_mess((" "))
        for i in list(self.resultats.values()):
            i.show_linked_concepts()

    def test_table(self, obj):
        """test si la table est composee de fonctions et si ces
        fonctions sont temporelles ou frequentielles"""
        table_py = obj.EXTR_TABLE()
        paras = table_py.para
        if "FONCTION_C" in paras:
            type_fonc = "FONCTION_C"
        elif "FONCTION" in paras:
            type_fonc = "FONCTION"
        else:
            return
        toto = table_py[type_fonc]
        try:
            nom_fonc = toto.values()[type_fonc][0]
        except IndexError:
            return  # deja vu : des tables avec lacolonne fonction vide...
        __FONC = RECU_FONCTION(
            TABLE=obj, NOM_PARA_TABL=type_fonc, FILTRE=_F(NOM_PARA=type_fonc, VALE_K=nom_fonc)
        )

        nom_para = __FONC.Parametres()["NOM_PARA"]
        if nom_para == "INST":
            return "tempo"
        elif nom_para == "FREQ":
            return "intespec"
        else:
            return

    def get_mode_meca(self, name):
        """!Renvoie un objet resultat identifie par son nom"""
        return self.mode_meca[name]

    def get_model(self, name):
        """!Renvoie un modele"""
        return self.modeles[name]

    def get_inter_spec(self, name):
        return self.inter_spec[name]

    def get_matr(self, name):
        """!Renvoie une matrice de masse ou raideur ou None"""
        return self.matrices.get(name)

    def get_cara_elem(self, name):
        """recup d'une sd resu dans la liste ci-dessus"""
        return self.cara_elem[name]

    def get_cham_mater(self, name):
        """recup d'une sd resu dans la liste ci-dessus"""
        return self.cham_mater[name]


#
#
#                          PETITS UTILITAIRES
#
#
def crea_champ(resu, ind_mod):
    """!Extrait les champs de deplacement d'une sd_resultat aster
    a partir des DDL de mesure pour un mode donne.
    Ces DDL sont identiques a ceux de la macro OBSERVATION
    ayant servi a obtenir le resultat."""
    __CHANO = CREA_CHAMP(
        TYPE_CHAM="NOEU_DEPL_R",
        OPERATION="EXTR",
        RESULTAT=resu,
        NOM_CHAM="DEPL",
        NUME_ORDRE=ind_mod,
    )

    vale, _ = __CHANO.getValuesWithDescription()

    return vale


def nume_ddl_phy(resu):
    """Fabrication d'une numerotation des DDL actifs associes a une sd resultat
    retourne ['N1_DX','N1_DY','N1_DZ','N2_DZ','N3_DZ'...] dans l'ordre alpha-numerique"""

    __CHAMP0 = CREA_CHAMP(
        RESULTAT=resu.obj, OPERATION="EXTR", NUME_ORDRE=1, TYPE_CHAM="NOEU_DEPL_R", NOM_CHAM="DEPL"
    )
    _, description = __CHAMP0.getValuesWithDescription()
    nb_ddl = len(description[0])

    # Recherche des noms des noeuds associes a leur numero
    maya = resu.maya
    nume = []
    for ind in range(nb_ddl):
        nom_no = maya.getNodeName(description[0][ind])
        nume.append(nom_no.strip() + "_" + description[1][ind])

    return nume


def nume_ddl_gene(resu, extract_mode=None):
    """
    Cree le meme vecteur de numerotation avec une numerotation par modes.
    Retourne "MO"+#mode
    """
    modes = []
    nume_mode = resu.get_modes_data()["NUME_MODE"]
    for mod in nume_mode:
        modes.append("MO" + str(int(mod)))
    return modes


def CreaTable(mcfact, titre, paras_out, mess):
    """!Sortie des donnees sous forme de sd_interspectre"""
    tablesOut = paras_out["TablesOut"]
    register = paras_out["Register"]
    compteur = paras_out["ComptTable"]
    paras_out["ComptTable"] = paras_out["ComptTable"] + 1

    if paras_out["ComptTable"] > len(tablesOut):
        mess.disp_mess("!! Il n'y a plus de noms de concepts     !!")
        mess.disp_mess("!! disponibles pour sortir des resultats !!")
        mess.disp_mess(" ")
        return

    tab = DEFI_INTE_SPEC(PAR_FONCTION=mcfact, TITRE=titre)
    register(tab, tablesOut[compteur])

    mess.disp_mess("Les resultats sont sauves dans la table " + tablesOut[compteur].userName)
    mess.disp_mess("Cette table porte pour titre : " + titre)
    mess.disp_mess(" ")

    return paras_out
