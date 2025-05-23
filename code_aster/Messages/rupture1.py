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
    2: _(
        """
Formalisme de déformation différents au moins pour la maille %(k3)s :
 - type de déformation extrait de la structure de données Résultat   : %(k1)s
 - type de déformation fourni à l'opérateur CALC_G : %(k2)s

--> Risques & conseils :
Vous effectuez un post-traitement de mécanique de la rupture à partir d'un formalisme
de déformation autre que PETIT et différent de celui utilisé pour le calcul mécanique.
Ceci est hors du champ de validation des opérateurs.
La démarche la plus fiable de calcul de G est de réaliser un second calcul en petites
déformations, de calculer G et de faire une équivalence en ouverture : le G en grandes
déformations sera égal au G en petites déformations lorsque les ouvertures de défaut
sont identiques.
"""
    ),
    3: _(
        """
Formalisme de déformation non valide au moins pour la maille %(k2)s :
 - type de déformation extrait de la structure de données Résultat   : %(k1)s


--> Risques & conseils :
Vous effectuez un post-traitement de mécanique de la rupture à partir d'un formalisme
de déformation autre que PETIT ou PETIT_REAC.
Ceci est hors du champ de validation des opérateurs.
La démarche la plus fiable de calcul de G est de réaliser un second calcul en petites
déformations, de calculer G et de faire une équivalence en ouverture : le G en grandes
déformations sera égal au G en petites déformations lorsque les ouvertures de défaut
sont identiques.
"""
    ),
    4: _(
        """
La commande CALC_G ne traite pas le cas des fonds doubles.
"""
    ),
    5: _(
        """
Le paramètre R_INF automatiquement choisi vaut %(r1)f.
Le paramètre R_SUP automatiquement choisi vaut %(r2)f.
"""
    ),
    6: _(
        """
Le rayon R_SUP (ou R_SUP_FO) doit obligatoirement être supérieur au rayon
R_INF (respectivement R_INF_FO).
-> Risque et Conseil :
Veuillez revoir votre mise en données.
"""
    ),
    7: _(
        """
Il n'est pas possible de calculer automatiquement R_INF et R_SUP en cas
de modélisation FEM avec une fissure en configuration décollée.
-> Risque et Conseil :
Veuillez indiquer les mots-clés R_INF et R_SUP (ou R_INF_FO et R_SUP_FO).
"""
    ),
    8: _(
        """
Pour une modélisation 2D, le front de fissure doit être défini par un
groupe contenant un noeud (mot-clé GROUP_NO). Le mot clé GROUP_MA
est réservé aux modélisations 3D.
"""
    ),
    9: _(
        """
Pour une modélisation 3D, le front de fissure doit être défini par un
groupe de mailles (mot-clé GROUP_MA). Le mot clé GROUP_NO est réservé
aux modélisations 2D.
"""
    ),
    10: _(
        """
Le fond de fissure ne doit être défini que par un noeud.
-> Risque et Conseil :
Veuillez revoir le contenu du mot-clé GROUP_NO ou NOEUD ou FOND_FISS.
"""
    ),
    11: _(
        """
Il faut un mot clé parmi FOND_FISS ou FISSURE.
Veuillez le renseigner.
"""
    ),
    12: _(
        """
Le champ de contrainte initiale n'est pas du bon type.
En FEM (fissure représentée dans le maillage), il doit être de type ELNO, NOEU, ou ELGA.
En X-FEM (fissure non représentée dans le maillage), il doit être de type ELGA.
"""
    ),
    13: _(
        """
La recherche du matériau associé au front de fissure a échouée. Cela arrive notamment quand
des noeuds du front du fissure sont situés à la frontière entre plusieurs matériaux.

Solution : Renseigner le mot-clé MATER de POST_K1_K2_K3 avec le matériau souhaité.

"""
    ),
    15: _(
        """
Le résultat n'est pas un EVOL_NOLI.
"""
    ),
    16: _(
        """
La différence entre la taille maximale et la taille minimale des mailles connectées aux
noeuds du fond de fissure est importante.
La taille minimale vaut : %(r1)f
La taille maximale vaut : %(r2)f
-> Risque et Conseil :
Il a été choisi de multiplier par quatre la taille maximale des mailles connectées aux
noeuds du fond de fissure pour calculer le paramètre R_SUP et par deux, pour calculer
le paramètre R_INF. Or, si ces tailles sont importantes, vous risquez de post-traiter
vos résultats sur une zone trop grande autour du fond de fissure et d'obtenir des
valeurs de taux de restitution d'énergie et de facteurs d'intensité moins précises.
Vérifiez que les valeurs de R_INF et de R_SUP calculées sont licites.
Sinon, veuillez spécifier directement les valeurs de R_INF et de R_SUP ou bien revoir
le maillage de manière à rendre les mailles proches du fond de fissure de taille homogène.
"""
    ),
    17: _(
        """
L'association: lissage du champ THETA par les polynômes de Lagrange
               lissage de G autre que par des polynômes de Lagrange
n'est pas possible.
-> Risque et Conseil :
Veuillez consulter la documentation U4.82.03 pour déterminer une
association satisfaisante.
"""
    ),
    18: _(
        """
Les mots-clés R_INF_FO et R_SUP_FO ne peuvent être employés dans le cas 2D.
-> Risque et Conseil :
Veuillez les remplacer par R_INF ou R_SUP, ou bien ne rien indiquer afin que les rayons
de la couronne soient calculés avec les données du maillage.
"""
    ),
    19: _(
        """
L'utilisation du mot-clé facteur %(k1)s n'est compatible qu'avec une modélisation %(k2)s.
Veuillez vérifiez vos données.
"""
    ),
    20: _(
        """
Votre étude comporte une charge de type PRE_EPSI ou une variable de commande EPSA.
Ceci est incompatible avec la présence d'une contrainte initiale dans le calcul de G
(mot clé SIGM_INIT de l'opérateur CALC_G).
-> Risque et Conseil :
On ne peut pas faire de calcul de G en introduisant simultanément une contrainte
initiale ET une déformation initiale. Veuillez revoir les données.
"""
    ),
    21: _(
        """
La liste de taille n'a pas la taille de la liste des groupes mailles.
Vérifiez vos données.
"""
    ),
    22: _(
        """
Les mailles volumiques caractérisant les zones de calcul doivent absolument être des
hexaèdres.
Vérifiez votre maillage ou vos données.
"""
    ),
    23: _(
        """
CALC_G - option CALC_G : détection de chargements non nuls sur l'axe,
le calcul est impossible.
-> Risque et Conseil :
En 2D axisymétrique, le calcul de G n'est pas possible pour les éléments de l'axe de
symétrie si un chargement est imposé sur ceux-ci.
Modifier les couronnes R_INF et R_SUP pour qu'elles soient toutes les deux
plus petites que le rayon du fond de fissure.
"""
    ),
    24: _(
        """
L'option CALC_K_G est incompatible avec les comportements incrémentaux,
avec les comportements non linéaires et avec la déformation GREEN_LAGRANGE.
"""
    ),
    25: _(
        """
Il faut affecter les éléments de bord (E et NU) pour le calcul des fic.
-> Risque et Conseil :
Veuillez revoir la mise en données des opérateurs DEFI_MATERIAU
et AFFE_MATERIAU.
"""
    ),
    26: _(
        """
La masse volumique RHO n'a pas été définie.

-> Risque et Conseil :
En présence de forces de pesanteur ou de rotation, ou pour le calcul de
l'option CALC_K_G avec un résultat de type MODE_MECA, il est indispensable
de fournir la masse volumique du matériau considéré.
RHO doit être défini dans l'opérateur DEFI_MATERIAU.
"""
    ),
    29: _(
        """
Au moins une des mailles caractérisant les zones de calcul a une forme trop
trapézoïdale.
-> Risque et Conseil :
Le calcul de la surface de sa face appartenant au plan de symétrie de
l'entaille risque d'être altéré et par conséquent celui de GP également.
Veuillez vérifier votre maillage.
"""
    ),
    30: _(
        """
Votre étude comporte une charge de type PRE_EPSI et une variable de commande EPSA.
Le calcul de G ne peut être réalisé.
-> Risque et Conseil :
On ne peut pas faire de calcul de G en introduisant simultanément une pré-déformation
ET une déformation initiale. Veuillez revoir les données.
"""
    ),
    31: _(
        """
DELTA_K_EQ contient une ou plusieurs valeurs négatives.
Cela peut provenir d'un calcul erroné de DELTA_K_EQ par l'utilisateur ou bien de l'utilisation
de l'opérande CUMUL='MODE_I' dans POST_RUPTURE dans le cas où K_I aurait des valeurs négatives.
-> Risque et Conseils
Vérifiez le calcul de DELTA_K_EQ.
"""
    ),
    35: _(
        """
Les directions normales au plan de la fissure entre les points %(i1)d et %(i2)d successifs du fond forment un angle
supérieur à 10 degrés.
-> Risque et Conseils
L'interpolation des sauts de déplacements est basée sur les champs singuliers
correspondants à une fissure plane. La fissure utilisée ici est trop irrégulière et
il y a donc un risque d'obtenir des résultats imprécis.
"""
    ),
    38: _(
        """
La fissure contient %(i1)d fond(s) et le calcul est demandé pour le fond numéro %(i2)d.
-> Risque et Conseil :
Vérifier le paramètre défini sous le mot clé NUME_FOND de POST_K1_K2_K3.
"""
    ),
    39: _(
        """
La récupération des contraintes à partir de la structure de données Résultat n'est permise que si les fissures sont maillées.
-> Risque et Conseil :
Veillez à ne pas vous servir de FISSURE avec le mot-clé CALCUL_CONTRAINTE.
"""
    ),
    41: _(
        """
Attention, dans CALC_G vous utilisez le mot clef RELATION.
La relation est normalement récupéré a partir du calcul mécanique dans
MECA_STATIQUE ou STAT_NON_LINE.
"""
    ),
    42: _(
        """
 Lois de comportement différentes au moins pour la maille %(k3)s :
 - loi de comportement extraite de la structure de données Résultat   : %(k1)s
 - loi de comportement fournie à l'opérateur CALC_G : %(k2)s

--> Risques & conseils :
On doit généralement utiliser la même loi de comportement entre le calcul et le
post-traitement. On peut utiliser deux comportements différents, mais alors
l'utilisateur doit être vigilant sur l'interprétation des résultats(cf.U2.05.01).
Si plusieurs comportements sont définis sur la structure, le comportement à
indiquer dans CALC_G est celui du matériau dans lequel la fissure se développe.
Dans ce cas, ce message d'alarme est quand même émis mais le résultat est bien cohérent.
Un post-traitement élastique non-linéaire d'un calcul élastoplastique est
admissible, si le chargement est radial et monotone. La comparaison du G calculé
à partir des contraintes issues de STAT_NON_LINE (option CALC_CONTRAINTE='NON')
ou à partir des contraintes recalculées avec la loi de comportement
(CALC_CONTRAINTES='OUI') peut fournir une indication sur le respect de ces
hypothèses.
"""
    ),
    44: _(
        """
Les paramètres K1 et/ou K2 et/ou G sont absents du tableau des facteurs d'intensité des
contraintes fourni.
-> Risque et Conseil :
Le tableau des facteurs d'intensité des contraintes doit absolument contenir ces trois
paramètres ainsi que K3 en 3D. Veuillez vérifier le contenu de votre tableau.
"""
    ),
    45: _(
        """
 Formalismes de déformations différents au moins pour la maille %(k3)s :
 - type de déformation extrait de la structure de données Résultat   : %(k1)s
 - type de déformation fourni à l'opérateur CALC_G : %(k2)s

--> Risques & conseils :
Vous avez choisi d'effectuer un post-traitement de mécanique de la rupture en petites déformations (seules
valables) alors que le calcul mécanique ne l'était pas. Cela demande d'être vigilant sur les résultats (cf.U2.05.01).
Si plusieurs types de déformation sont définis sur la structure et que la fissure se développe dans une partie en
petites déformations, ce message d'alarme est quand même émis mais le résultat est bien cohérent.
"""
    ),
    46: _(
        """
Le taux de restitution d'énergie G est négatif ou nul sur certains des noeuds du fond de fissure :
le calcul de propagation est impossible.
-> Risque et Conseil :
Vérifier les paramètres du calcul de G (rayons des couronnes ou abscisse curviligne
maximale, type de lissage, ...).
"""
    ),
    47: _(
        """
Vous demandez un calcul de G en post-traitement d'un calcul élastoplastique. Ceci n'est valable que
si votre CHARGEMENT est MONOTONE PROPORTIONNEL.
Si tel est le cas, renseignez, dans CALC_G, l'option RELATION = ELAS_VMIS_XXX pour un calcul de G.
"""
    ),
    49: _(
        """
Vous ne pouvez pas utiliser CALC_G avec comme relation de comportement: %(k2)s ,alors que lors
de la résolution du Problème mécanique vous utilisez une relation de comportement de type: %(k1)s
"""
    ),
    51: _(
        """
PROPA_FISS / METHODE = 'MAILLAGE' : les noeuds définissant la fissure initiale ne sont
pas ordonnés. Vérifiez le maillage donné en entrée (MAIL_ACTUEL).
"""
    ),
    52: _(
        """
PROPA_FISS / METHODE = 'INITIALISATION' : les deux vecteurs VECT_X et VECT_Y
définissant la fissure initiale, de forme elliptique, ne sont pas orthogonaux. Vérifiez
les données d'entrée.
"""
    ),
    53: _(
        """
L'instant %(r1)f n'appartient pas au résultat %(k1)s.
"""
    ),
    54: _(
        """
Les champs de contraintes et de déformations ne sont pas de même taille. Vérifiez que votre
calcul mécanique s'est bien passé.
"""
    ),
    55: _(
        """
Problème dans la liste d'instants du résultats: 2 instants consécutifs sont égaux.
"""
    ),
    56: _(
        """
La contrainte de référence est nulle à l'instant %(r1)f.
"""
    ),
    57: _(
        """
Problème dans la dimension du modèle. POST_BORDET ne supporte pas les raccords 2D-3D
"""
    ),
    58: _(
        """
L'utilisation de POST_BORDET n'est possible qu'avec 1 seul MODELE et 1 seul
CHAM_MATERIAU
"""
    ),
    59: _(
        """
La table %(k1)s ne contient pas le paramètre %(k2)s.
"""
    ),
    60: _(
        """
Le critère 'K2_NUL' donne des mauvais résultats pour des angles supérieurs à 60 degrés.
Il se peut que le signe de l'angle soit faux.
Conseil : utiliser le critère par défaut.
"""
    ),
    61: _(
        """
Impossible de réaliser le comptage sur les quantités demandées car
le nombre de cycles pour chacune d'elles est différent.
Conseil : limiter le comptage des cycles à une seule quantité (K_EQ par exemple).
"""
    ),
    62: _(
        """
Pour l'opération %(k1)s, la table doit être réentrante (reuse obligatoire).
"""
    ),
    63: _(
        """La récupération des contraintes à partir de la structure de données Résultat
en présence d'un état initial n'est pas permise.
Pour l'opération %(k1)s, la table ne doit pas être réentrante (reuse interdit).
"""
    ),
    64: _(
        """
Pour le comptage %(k1)s, la table doit comporter uniquement 1 instant/NUME_ORDRE (ou aucun).
Or la table %(k2)s contient %(i1)d instants/NUME_ORDRE.
Conseil : Vérifier la table en entrée ou utiliser un autre type de comptage des cycles.
"""
    ),
    65: _(
        """
La table %(k1)s ne doit pas contenir le paramètre %(k2)s.
"""
    ),
    66: _(
        """
L'opération %(k1)s nécessite une seule table sous le mot-clé TABLE, or il y en a %(i1)d.
"""
    ),
    67: _(
        """
Les caractéristiques du matériau ne peuvent dépendre que des variables de commande TEMP et NEUT1.
-> Conseil :
Veuillez revoir les données du matériau.
"""
    ),
    68: _(
        """
La macro-commande POST_RUPTURE ne fonctionne pas quand les paramètres matériau ne sont pas constants.
"""
    ),
    69: _(
        """
Comportement incrémental :
La relation de comportement %(k1)s n'est pas prévue dans CALC_G
"""
    ),
    71: _(
        """
La différence relative moyenne entre G et G_IRWIN est: %(r1)f

--> Risques & conseils :
Pour un problème thermo-mécanique, on a pris la solution singulière d'un problème
purement mécanique. Le calcul des facteurs d'intensité des contraintes est approché.
Si la différence relative entre G et G_IRWIN est importante, le calcul des facteurs
d'intensité des contraintes K par CALC_K_G n'est plus valable.

"""
    ),
    73: _(
        """
CALC_G :
Il est interdit d'utiliser NB_POINT_FOND avec les lissages LEGENDRE, LAGRANGE_NO_NO
et MIXTE
-> Conseil :
Veuillez revoir le type de lissage utilisé
"""
    ),
    74: _(
        """
POST_K1_K2_K3 : Cet opérateur est incompatible avec les modélisations incompressibles. En effet, il n'est valide qu'en élasticité linéaire, où les modélisations incompressibles ne sont pas nécessaires.
"""
    ),
    75: _(
        """
CALC_G : On ne peut pas calculer les dérivées des fonctions singulières car on se trouve sur le fond de fissure.
"""
    ),
    76: _(
        """
CALC_G : Le champ fourni par CHAM_THETA_IN n'est pas correct. Les composants DIR_Z, ABSC_CURV et LONG doivent être nuls en 2D.
"""
    ),
    89: _(
        """
L'OPTION NB_COUCHE n'est pas compatible avec des éléments HEXA27 ou PENTA18'
"""
    ),
    90: _(
        """
CALC_G fonctionne uniquement avec des éléments incompressibles de type INCO_UPG. L'utilisation d'autres types d'éléments pourrait entraîner des résultats incorrects (INCO_UP, INCO_UPO, GRAD_INCO).
"""
    ),
}
