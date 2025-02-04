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

from ..Utilities import _

cata_msg = {
    1: _(
        """
Formulation continue.
Vous avez activé le frottement de Coulomb (FROTTEMENT='COULOMB') pourtant toutes les zones de contact
portent un coefficient de frottement nul.

Le frottement est donc désactivé.

Conseil : vérifiez que vous avez correctement défini le coefficient de frottement (mot-clé COULOMB) dans chaque zone.
"""
    ),
    2: _(
        """
La maille < %(k1)s > est de type POI1. C'est impossible dans la méthode continue.
"""
    ),
    3: _(
        """
Le frottement n'est autorisé qu'avec une formulation pénalisée sur la partie contact.
"""
    ),
    4: _(
        """
    Erreur d'utilisation de la méthode LAC.
  Le mot-clef < %(k1)s > doit avoir la même valeur sur toutes les zones
  de contact
"""
    ),
    5: _("""   Il y a %(i1)d points initialement en contact et %(i2)d points exclus."""),
    6: _("""   Il y a %(i1)d éléments initialement en contact."""),
    15: _(
        """
La direction d'appariement fixe données par le vecteur DIRE_APPA est nulle !
"""
    ),
    16: _(
        """
Contact méthode continue.
La méthode d'intégration n'est pas aux noeuds, le champ CONT_NOEU n'est pas créé.
"""
    ),
    18: _(
        """
Contact méthode continue.
La direction d'exclusion du frottement fournie pour la zone de contact numéro %(i1)d (%(r1)f,%(r2)f,%(r3)f) est perpendiculaire au plan de contact sur la maille %(k1)s.

Conseil :
   - Vérifiez le vecteur DIRE_EXCL_FROT. Sa projection sur le plan tangent de contact doit exister
     pour indiquer une direction à exclure.

"""
    ),
    19: _(
        """
La surface maître est fortement facétisée sur %(i1)d noeuds.
A cause des variations brusques de normale d'un noeud à l'autre, vous risquez d'avoir des problèmes de convergence du contact, surtout dans le cas de grands glissements relatifs des deux surfaces,.
Vous pouvez raffiner le maillage, utiliser un maillage quadratique (si c'est possible) ou activer le lissage (LISSAGE='OUI' dans DEFI_CONTACT).

Pour information, voici la liste des noeuds avec l'angle (en degré) mesuré entre la normale nodale et les normales aux éléments:
Noeud      Angle
"""
    ),
    23: _(
        """
Le vecteur normal est nul au niveau du projeté du point de contact de coordonnées
  (%(r1)f,%(r2)f,%(r3)f)
sur la maille %(k1)s,
Erreur de définition de la maille ou projection difficile. Contactez l'assistance dans ce dernier cas.
"""
    ),
    24: _(
        """
Le vecteur normal est nul sur la maille %(k1)s,
Erreur de définition de la maille ou projection difficile. Contactez l'assistance dans ce dernier cas.
"""
    ),
    26: _(
        """
Le vecteur normal est nul au niveau du noeud %(k1)s.
"""
    ),
    27: _(
        """
Le vecteur normal est nul au niveau de la maille %(k1)s.
"""
    ),
    31: _(
        """
Les vecteurs tangents sont nuls au niveau quand on projette le noeud esclave
%(k1)s sur la maille %(k2)s
"""
    ),
    32: _(
        """
Le mot-clef DIST_POUT ne fonctionne qu'avec des sections circulaires définies dans AFFE_CARA_ELEM.
 """
    ),
    35: _(
        """
Les vecteurs tangents sont nuls au niveau quand on projette le noeud esclave
%(k1)s sur la noeud maître %(k2)s
"""
    ),
    37: _(
        """
La section de la poutre n'est pas constante sur l'élément. On prend la moyenne.
 """
    ),
    39: _(
        """
Problème pour récupérer l'épaisseur de la coque pour la maille  %(k1)s
"""
    ),
    40: _(
        """
L'excentricité de la coque pour la maille %(k1)s ne peut pas être traitée
"""
    ),
    41: _(
        """
Problème pour récupérer l'excentricité de la coque pour la maille  %(k1)s
"""
    ),
    43: _(
        """
La normale que vous avez prédéfinie par (VECT_* = 'VECT_Y') sur la maille %(k1)s n'est pas utilisable en 2D.
Utilisez plutôt VECT_* = 'FIXE'
"""
    ),
    45: _(
        """
Pour les formulations discrètes et unilatérales, il faut utiliser la méthode de distribution centralisée "DISTRIBUTION=_F(METHODE='CENTRALISE')" dans AFFE_MODELE.
"""
    ),
    46: _(
        """
Pour la dynamique non linéaire implicite avec contact, il faut utiliser la méthode de distribution centralisée
"DISTRIBUTION=_F(METHODE='CENTRALISE')" dans AFFE_MODELE.
"""
    ),
    50: _(
        """
Avec l'option VECT_MAIT = 'FIXE', seule l'option NORMALE = 'MAIT' est possible.
"""
    ),
    51: _(
        """
Avec l'option VECT_MAIT = 'VECT_Y', seule l'option NORMALE = 'MAIT' est possible.
"""
    ),
    52: _(
        """
Avec l'option VECT_ESCL = 'FIXE', seule l'option NORMALE = 'ESCL' est possible.
"""
    ),
    53: _(
        """
Avec l'option VECT_ESCL = 'VECT_Y', seule l'option NORMALE = 'ESCL' est possible.
"""
    ),
    54: _(
        """
Le LISSAGE n'est possible qu'avec des normales automatiques VECT_ESCL='AUTO' et/ou VECT_MAIT='AUTO'.
"""
    ),
    54: _(
        """
Le LISSAGE n'est possible qu'avec des normales automatiques VECT_ESCL='AUTO' et/ou VECT_MAIT='AUTO'.
"""
    ),
    60: _(
        """
La maille %(k1)s est de type 'POI1', elle nécessite l'utilisation de l'option
NORMALE='FIXE' avec une normale non nulle.
"""
    ),
    61: _(
        """
La maille %(k1)s est de type poutre, elle nécessite la définition d'une base locale.
Utilisez NORMALE='FIXE' ou NORMALE='VECT_Y' dans DEFI_CONTACT.
"""
    ),
    75: _(
        """
La maille %(k1)s est de type 'POI1', elle ne peut pas être une maille maître.
"""
    ),
    96: _(
        """
Contact.
    -> Les surfaces en contact ont bougé de plus de 1%% depuis la dernière réactualisation.
       Or vous n'avez pas activé la réactualisation géométrique automatique dans la commande DEFI_CONTACT
       (REAC_GEOM='AUTOMATIQUE') ou bien vous utilisez le mode 'CONTROLE'
    -> Risque & conseil :
       Vos résultats risquent d'être faux, les mailles ne seront peut-être pas correctement appariées
       et des interpénétrations pourraient apparaître.
       Si vous avez volontairement négligé la non-linéarité géométrique de contact (pour des raisons
       de performance), nous vous invitons à vérifier visuellement qu'il n'y a effectivement
       pas d'interpénétrations.
"""
    ),
    98: _(
        """
Contact et détection de collision.
    Le contact semble "rasant", c'est-à-dire que le jeu est nul mais que la pression de contact est quasiment
    nulle également.
    -> Risque & Conseil :
      La découpe automatique du pas de temps ne sera probablement pas efficace. Si vous constatez de fortes oscillations
      de la vitesse ou du déplacement au point de contact pour cet instant de calcul malgré l'activation du mode de
      traitement automatique de la collision, vous pouvez essayer d'éviter l'instant pour lequel le contact "rasant"
      apparaît de deux manières différentes:
       - changer la discrétisation en temps initiale
       - changer la valeur de SUBD_INST
"""
    ),
    99: _(
        """
    DEFI_CONTACT/ALGO_CONT=PENALISATION/ADAPTATION=PENE_MAXI.
    Le critère de pénétration PENE_MAXI semble lâche.  Pour information :
        - La plus petite arête de la zone maître est  %(r1).2e.
        - La plus grande arête de la zone maître est  %(r2).2e.
    -> Risque & Conseils :
       Vous risquez d'avoir des résultats inattendus. Vérifiez que le niveau d'interpénétration dans le tableau de convergence est bien inférieure une longueur d'arête dans la zone de contact.
       S'il y a trop d'interpénétrations, il faudrait relancer le calcul avec PENE_MAXI <= %(r1).2e.
"""
    ),
}
