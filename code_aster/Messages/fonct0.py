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

from ..Utilities import _

cata_msg = {
    1: _(
        """
Le fichier %(k1)s existe déjà, on écrit à la suite.
"""
    ),
    2: _(
        """
Il n'y a pas de règles d'interpolation pour LIST_PARA/LIST_RESU,
LIST_PARA/LIST_RESU ne peut donc apparaître qu'une seule fois
et à la première occurrence de COURBE.
"""
    ),
    3: _(
        """
LIST_PARA et LIST_RESU n'ont pas la même taille.
"""
    ),
    4: _(
        """
FONC_X/FONC_Y ne peuvent pas être des nappes !
"""
    ),
    5: _(
        """
Au format 'TABLEAU', FONC_X/FONC_Y ne peut apparaître qu'une seule fois
et à la première occurrence de COURBE
"""
    ),
    6: _(
        """
Il n'y a pas de règles d'interpolation pour ABSCISSE/ORDONNEE,
ABSCISSE/ORDONNEE ne peut donc apparaître qu'une seule fois
et à la première occurrence de COURBE.
"""
    ),
    7: _(
        """
ABSCISSE et ORDONNEE n'ont pas la même taille.
"""
    ),
    8: _(
        """
Format inconnu : %(k1)s
"""
    ),
    9: _(
        """
Erreur lors de l'interpolation de la fonction '%(k1)s'.
"""
    ),
    10: _(
        """sur la maille '%(k1)s'
"""
    ),
    12: _(
        """
Une erreur s'est produite dans la recherche de l'intervalle des abscisses contenant la valeur %(r1)f.
"""
    ),
    13: _(
        """
Le type de la fonction '%(k1)s' est inconnu.
Seules les fonctions, nappes, fonctions constantes peuvent être traitées par %(k3)s.

  -> Débogage :
      le type est '%(k2)s'
"""
    ),
    14: _(
        """
Il n'y a pas assez de paramètres pour évaluer la fonction.
Seulement %(i1)d paramètre(s) sont fourni(s) alors que la fonction en réclame %(i2)d.
"""
    ),
    15: _(
        """
Il y a des doublons dans la liste des paramètres fournis :
   %(ktout)s
"""
    ),
    16: _(
        """
Les paramètres dont dépend la fonction sont :
   %(ktout)s
"""
    ),
    17: _(
        """
Alors que les paramètres fournis pour réaliser l'interpolation sont :
   %(ktout)s

L'interpolation n'est donc pas possible.
Conseil : Modifier les paramètres dont la fonction dépend pour rendre possible
l'interpolation.
"""
    ),
    18: _(
        """
La fonction n'a même pas un point !
"""
    ),
    19: _(
        """
On est hors du domaine de définition de la fonction.
On ne peut pas interpoler la fonction pour cette abscisse car le prolongement à gauche est exclus.
   abscisse demandée              : %(r1)f
   borne inférieure des abscisses : %(r2)f

  -> Risque & Conseil :
    Voir le mot-clé PROL_GAUCHE des commandes qui créent des fonctions.
"""
    ),
    20: _(
        """
On est hors du domaine de définition de la fonction.
On ne peut pas interpoler la fonction pour cette abscisse car le prolongement à droite est exclus.
   abscisse demandée              : %(r1)f
   borne supérieure des abscisses : %(r2)f

  -> Risque & Conseil :
    Voir le mot-clé PROL_DROITE des commandes qui créent des fonctions.
"""
    ),
    21: _(
        """
Erreur de programmation : type d'extrapolation inconnu.

  -> Débogage :
      le type d'extrapolation est '%(k1)s'
"""
    ),
    22: _(
        """
La fonction n'est définie qu'en un point. On ne peut pas l'interpoler en
plus d'un point si le prolongement n'est pas constant des deux cotés.

  -> Risque & Conseil :
    Voir les mots-clés PROL_GAUCHE/PROL_DROITE des commandes qui créent des fonctions.
"""
    ),
    23: _(
        """
La fonction n'est définie qu'en un point. On ne peut pas l'interpoler ailleurs
qu'en ce point si le prolongement n'est pas constant des deux cotés.

  -> Risque & Conseil :
    Voir les mots-clés PROL_GAUCHE/PROL_DROITE des commandes qui créent des fonctions.
"""
    ),
    24: _(
        """
On attend une fonction d'un seul paramètre.
La fonction '%(k1)s' est une fonction de %(i1)d paramètres.
"""
    ),
    25: _(
        """
Le type de la fonction '%(k1)s' est inconnu.
Seules les fonctions, nappes, fonctions constantes et formules sont
traitées par %(k3)s.

  -> Débogage :
      le type est '%(k2)s'
"""
    ),
    27: _(
        """
Un problème d'interpolation a été rencontré.
%(k1)s

  -> Risque & Conseil :
      Vérifier les valeurs fournies derrière le mot-clé 'INTERPOL' lors
      de la création de cette(ces) fonction(s).
"""
    ),
    28: _(
        """
Un problème concernant le nom des abscisses ou ordonnées a été rencontré.
Vous ne pouvez pas faire la transformée de Fourier d'une fonction dont les abscisses sont des fréquences,
   ou si la fonction est a valeurs complexes
Vous ne pouvez pas faire la transformée de Fourier inverse d'une fonction dont les abscisses sont des instants,
   ou si la fonction est a valeur réelle.
%(k1)s

  -> Risque & Conseil :
      Vérifier la valeur fournie derrière les mots-clés 'NOM_PARA'/'NOM_RESU' lors
      de la création de cette(ces) fonction(s).
"""
    ),
    29: _(
        """
Un problème concernant le prolongement de la (des) fonction(s) a été rencontré.
%(k1)s

  -> Risque & Conseil :
      Vérifier la valeur fournie derrière les mots-clés 'PROL_GAUCHE'/'PROL_DROITE'
      lors de la création de cette(ces) fonction(s).
"""
    ),
    30: _(
        """
Une erreur s'est produite lors de l'opération.
%(k1)s

  -> Débogage :
      %(k2)s

Remontée d'erreur (pour aider à l'analyse) :

%(k3)s

"""
    ),
    31: _(
        """
   Génération par défaut de trois amortissements :[%(r1)f, %(r2)f, %(r3)f]
"""
    ),
    32: _(
        """
   Génération par défaut de %(i1)dfréquences :
   %(k1)s
"""
    ),
    33: _(
        """
   SPEC_OSCI, la norme ne peut être nulle.
"""
    ),
    35: _(
        """
   INTERPOL_FFT : Le pas de temps PAS_INST demandé pour le calcul de l'interpolation
   du signal est plus grand que le pas de temps du signal initial.
"""
    ),
    36: _(
        """
   SPEC_OSCI, la méthode choisie suppose des amortissements sous critiques,
   (inférieurs à 1).
"""
    ),
    37: _(
        """
 calcul du MAX, la liste de fonctions n'est pas
 homogène en type (fonctions et nappes)
"""
    ),
    38: _(
        """
 Calcul du MAX, la liste de fonctions n'est pas homogène
 en label NOM_PARA :%(k1)s
"""
    ),
    39: _(
        """
 Calcul du MAX, la liste de fonctions n'est pas homogène
 en label NOM_RESU :%(k1)s
"""
    ),
    40: _(
        """
 Intensité spectrale, avant de calculer l'intensité spectrale,
 il est prudent de vérifier la norme de la nappe sur laquelle
 porte le calcul, ceci peut être une source d'erreurs.
"""
    ),
    41: _(
        """
 Le fichier %(k1)s est introuvable.
"""
    ),
    42: _(
        """
Erreur lors de la lecture des blocs de valeurs :
   %(k1)s
"""
    ),
    43: _(
        """
Les fréquences doivent être strictement positives.
"""
    ),
    44: _(
        """
Les abscisses de la fonction %(k1)s ne sont pas strictement croissantes.
"""
    ),
    45: _(
        """
Les abscisses de la fonction %(k1)s ne sont pas croissantes.
"""
    ),
    46: _(
        """
Les abscisses de la fonction %(k1)s ne sont pas décroissantes.
"""
    ),
    47: _(
        """
Les abscisses de la fonction %(k1)s ne sont pas strictement décroissantes.
"""
    ),
    48: _(
        """
La fonction ou formule ne doit avoir qu'une ou deux variables.
"""
    ),
    49: (
        """
La nappe ou formule a deux paramètres. Il faut renseigner le mot-clé NOM_PARA_FONC
et soit VALE_PARA_FONC, soit LIST_PARA_FONC.
"""
    ),
    50: _(
        """
Seules les formules à une variable peuvent être traitées directement par IMPR_FONCTION.

La formule '%(k1)s' dépend de %(i1)d paramètres.

  -> Risque & Conseil :
      - Si votre formule dépend de 2 paramètres, utilisez CALC_FONC_INTERP pour produire
        une nappe puis appeler IMPR_FONCTION.
      - Si votre formule dépend de 3 paramètres ou plus, vous devez d'abord créer une
        nouvelle formule à un seul paramètre (et appelé IMPR_FONCTION) ou à 2 paramètres
        et passer par CALC_FONC_INTERP puis IMPR_FONCTION.
"""
    ),
    51: _(
        """
L'écart entre le pas de temps demandé (PAS_INST) et celui obtenu est supérieur
à la précision souhaitée :
   Pas de temps du signal de sortie obtenu   : %(r1)f
   Pas de temps du signal de sortie souhaité : %(r2)f
   Écart relatif                             : %(r3).2f %%

"""
    ),
    52: _(
        """
Conseils :
  Si le problème reporté ci-dessus ressemble à 'NameError: 'XXX'...',
  vérifiez que le paramètre 'XXX' fait bien partie des paramètres de définition de
  la formule (mot clé FORMULE / NOM_PARA).
"""
    ),
    53: _(
        """sur le noeud '%(k1)s'
"""
    ),
    54: (
        """
Nombre de paramètres fournis : %(i1)d
Noms des paramètres fournis  : %(ktout)s
"""
    ),
    55: _(
        """
  La liste des bornes de l'intervalle n'est pas cohérente.
  Elle doit comporter un nombre pair de valeurs.
"""
    ),
    56: _(
        """
  La borne inférieurs doit être inférieure à la borne supérieure.
  Veuillez revoir la saisie du mot-clé INTERVALLE.
"""
    ),
    57: _(
        """
Le polynôme est de la forme :
    a[0] x^N + a[1] x^(N-1) + a[2] x^(N-2) + ... + a[N]

avec :
   %(k1)s

"""
    ),
    58: _(
        """
Erreur lors de la vérification des noms des paramètres.
Le nom du premier paramètre de la formule en entrée (%(k1)s) est '%(k2)s'.

Or vous avez demandé à créer une nappe avec NOM_PARA='%(k3)s'.
"""
    ),
    59: _(
        """
Erreur lors de la vérification des noms des paramètres.
Le nom du paramètre de la nappe en entrée (%(k1)s) est '%(k2)s'.

Or vous avez demandé à créer une nappe avec NOM_PARA='%(k3)s'.
"""
    ),
    60: _(
        """
Erreur lors de la vérification des noms des paramètres.
Le nom du deuxième paramètre de la formule en entrée (%(k1)s) est '%(k2)s'.

Or vous avez demandé à créer une nappe avec NOM_PARA_FONC='%(k3)s'
"""
    ),
    61: _(
        """
Erreur lors de la vérification des noms des paramètres.
Le nom du paramètre des fonctions de la nappe en entrée (%(k1)s) est '%(k2)s'.

Or vous avez demandé à créer une nappe avec NOM_PARA_FONC='%(k3)s'
"""
    ),
    62: _(
        """
Création de la fonction '%(k1)s'.
"""
    ),
    63: _(
        """
Création d'une fonction de la nappe '%(k1)s'.
"""
    ),
    64: _(
        """
Les abscisses ne sont pas strictement monotones.
"""
    ),
    65: _(
        """
Les abscisses ont été réordonnées.
"""
    ),
    66: _(
        """
L'ordre des abscisses a été inversé.
"""
    ),
    67: _(
        """
Le nombre de valeurs est différent du nombre de paramètres
"""
    ),
    70: _(
        """
Erreur lors de l'évaluation de la formule.
La remontée d'erreur suivante peut aider à comprendre où se situe l'erreur :
%(k1)s
"""
    ),
    71: _(
        """
Type invalide pour le mot-clé '%(k1)s'.
Il semble que ce ne soit pas une structure de données connue.

Conseil
-------
    Vérifiez que les objets utilisés dans cette commande n'ont pas été
    créé en tant que liste Python. Cela peut se produire si la commande
    précédente se termine par une virgule.
"""
    ),
    72: _(
        """
Le nombre de fréquence imposée est supérieure, on augmente le nombre de points imposé de '%(k1)s' à '%(k2)s'.
"""
    ),
    73: _(
        """
L'algorithme n'a pas converge après %(k1)s itérations , il reste %(k2)s points au lieu de %(k3)s points demandés.
"""
    ),
    74: _(
        """
Les nappes fournit n'ont pas les mêmes valeurs d'amortissements
"""
    ),
    75: _(
        """
En mode CONCEPTION, un seul nombre de fréquence NB_FREQ_LISS est pris en compte. En mode VERIFICATION, les deux premières valeurs de NB_FREQ_LISS sont pris en compte. Les autres valeurs sont ignorés.
"""
    ),
    76: _(
        """
Le nombre de coefficients d'élargissement est différent du nombre de nappes ou tables fournies.
"""
    ),
    77: _(
        """
La fréquence fournie %(k1)s n'appartient pas à la liste des fréquences initiales du spectre.
"""
    ),
    78: _(
        """
Le déplacement maximal %(r1)f souhaité est trop élevé ou trop petit. Le déplacement spectral à la fréquence %(r2)f vaut %(r3)f. On ne peut pas prolonger.
"""
    ),
    80: _(
        """
La fonction %(k1)s n'a pas les bons paramètres d'accès. C'est soit :
   - INST et EPAIS       : Instant et l'épaisseur de la coque/plaque,
   - INST et EXCENT      : Instant et l'épaisseur de la coque/plaque en tenant compte
                           de l'excentrement,
   - INST et X et Y et Z : Instant et coordonnées du centre de gravité de la couche.
"""
    ),
}
