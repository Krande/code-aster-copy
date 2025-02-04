# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
 La somme de matrices distribuées n'ayant pas le même profil est impossible
"""
    ),
    2: _(
        """
 La matrice est symétrique.
"""
    ),
    3: _(
        """
 La matrice n'est pas symétrique.
"""
    ),
    4: _(
        """
 erreur LAPACK (ou BLAS) au niveau de la routine  %(k1)s
  le paramètre numéro  %(i1)d
  n'a pas une valeur cohérente %(i2)d
"""
    ),
    6: _(
        """
 Résolution MULTI_FRONTALE :
 problème dans le traitement des résultats
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    10: _(
        """
  le nombre de noeuds de la structure   :  %(i1)d
  la base utilisée est              :  %(k1)s
  les caractéristiques élémentaires :  %(k2)s
  diamètre de la structure          :  %(r1)f
  type de pas                       :  %(i2)d
"""
    ),
    11: _(
        """
  le profil de vitesse de la zone :  %(k1)s
  type de réseau de la zone       :  %(i1)d
"""
    ),
    13: _(
        """
  le noeud d'application            :  %(k1)s
  la base utilisée est              :  %(k2)s
  les caractéristiques élémentaires :  %(k3)s
  diamètre de la structure          :  %(r1)f
  type de configuration             :  %(k4)s
  le coefficient de masse ajoutée   :  %(r2)f
  le profil de masse volumique      :  %(r3)f
"""
    ),
    14: _(
        """
    pas de couplage pris en compte
"""
    ),
    15: _(
        """
   pour le concept  %(k1)s, le mode numéro  %(i1)d
"""
    ),
    16: _(
        """
  de fréquence  %(r1)f
"""
    ),
    17: _(
        """
  de charge critique  %(r1)f
"""
    ),
    18: _(
        """
  a une norme d'erreur de  %(r1)f  supérieure au seuil admis  %(r2)f.
"""
    ),
    20: _(
        """
  est en dehors de l'intervalle de recherche : [ %(r1)f,  %(r2)f ].
"""
    ),
    21: _(
        """
 AMELIORATION='OUI' n'est pas compatible avec les options PROCHE, SEPARE ou AJUSTE.
"""
    ),
    23: _(
        """
   pour le concept  %(k1)s,
"""
    ),
    24: _(
        """
  dans l'intervalle [%(r1)f  ,  %(r2)f]
  il y a théoriquement  %(i1)d fréquence(s) propres()
  et on en a calculé  %(i2)d.
"""
    ),
    25: _(
        """
  dans l'intervalle [%(r1)f  ,  %(r2)f]
  il y a théoriquement  %(i1)d charge(s) critique(s)
  et on en a calculé  %(i2)d.
"""
    ),
    26: _(
        """
 Ce problème peut apparaître lorsqu'il y a des modes multiples (structure avec symétries),
 une forte densité modale ou si vous avez juste choisi d'affiner quelques modes parmi un groupe
 de modes (options PROCHE, SEPARE ou AJUSTE).
"""
    ),
    27: _(
        """
 La valeur du SHIFT %(r1)f coïncide avec une fréquence propre.
"""
    ),
    28: _(
        """
 les nombres de termes des matrices RIGI et MASSE différent
 celui de la matrice MASSE vaut :  %(i1)d
 celui de la matrice RIGI  vaut :  %(i2)d

"""
    ),
    29: _(
        """
 le nombre d'amortissements réduits est trop grand
 le nombre de modes propres vaut  %(i1)d
 et le nombre de coefficients :   %(i2)d
 on ne garde donc que les %(i3)d premiers coefficients

"""
    ),
    30: _(
        """
 le nombre d'amortissements réduits est insuffisant, il en manque :  %(i1)d,
 car le nombre de modes vaut :  %(i2)d
 on rajoute  %(i3)d amortissements réduits avec la valeur du dernier mode propre.
"""
    ),
    31: _(
        """
  incohérence :
     %(i1)d
    %(i2)d
    %(i3)d
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    32: _(
        """
  erreur de type  différent de -1 ou -2  %(i1)d
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    33: _(
        """
 un ddl bloqué a au moins 2 LAMBDA1 ou 2 LAMBDA2
 le ddl bloqué est  %(i1)d
Ce message est un message d'erreur développeur.
Contactez le support technique.

"""
    ),
    34: _(
        """
 incohérence des multiplicateurs de Lagrange
 DDL %(i1)d
 LAMBDA1 %(i2)d
 LAMBDA2 %(i3)d
"""
    ),
    35: _(
        """
 erreur programmeur
 le LAMBDA2  %(i1)d a moins de 2 voisins
 il faut le LAMBDA1 et au moins un DDL
Ce message est un message d'erreur développeur.
Contactez le support technique.

"""
    ),
    36: _(
        """
 Problème dans le calcul des DDL :
 NUM devrait être égal à n1 :
 NUM = %(i1)d , n1 = %(i2)d
 impression des multiplicateurs de Lagrange
"""
    ),
    37: _(
        """
 NUME_DDL incohérence des multiplicateurs de Lagrange
  DDL     %(i1)d
  LAMBDA1 %(i2)d
  LAMBDA2 %(i3)d
"""
    ),
    38: _(
        """
 nombre de relations linéaires %(i1)d
"""
    ),
    39: _(
        """
 LAMBDA1 de R linéaire : %(i1)d
 LAMBDA2 de R linéaire : %(i2)d
"""
    ),
    40: _(
        """
 Données erronées
"""
    ),
    41: _(
        """
 pas de mode statique pour  le noeud :  %(k1)s  et sa composante :  %(k2)s

"""
    ),
    42: _(
        """
 pour les modes statiques :
 on attend un :  %(k1)s
 noeud :  %(k2)s
 composante   :  %(k3)s

"""
    ),
    43: _(
        """
 champ inexistant.
 champ :  %(k1)s
 noeud :  %(k2)s
 composante   :  %(k3)s

"""
    ),
    48: _(
        """
 incohérence de certains paramètres modaux propres à ARPACK
  numéro d'erreur  %(i1)d

"""
    ),
    49: _(
        """
Le nombre de modes propres calculés (%(i1)d) est inférieur au nombre
de fréquences demandées (%(i2)d).
 --> Le calcul continue, on ne prend que les %(i1)d premiers modes.

Conseil
-------
La prochaine fois, relancez en augmentant la taille de l'espace de projection (mots-clés
COEF_DIM_ESPACE ou DIM_SOUS_ESPACE).
Vous pouvez aussi relancer en demandant moins de modes.
Si votre problème est fortement amorti, il est possible que des modes propres non calculés soient
sur amortis, diminuez alors le nombre de modes demandés.
"""
    ),
    50: _(
        """
Le nombre de modes propres calculés (%(i1)d) est supérieur au nombre
de fréquences demandées (%(i2)d).

--> Le calcul continue, on ne prend que les %(i2)d premiers modes.

"""
    ),
    51: _(
        """
 La valeur propre numéro %(i1)d a une partie imaginaire non négligeable.
 Partie réelle     = %(r1)12.5E
 Partie imaginaire = %(r2)12.5E

 Ce phénomène numérique est fréquent sur les valeurs propres en bordure du spectre recherché.

 Conseil
 -------
  En cas de problème, la prochaine fois, relancez en demandant moins de modes ou en augmentant la
  taille de l'espace de projection (mots-clés COEF_DIM_ESPACE ou DIM_SOUS_ESPACE).
"""
    ),
    52: _(
        """
 LAIGLE: Erreur
   - Non convergence à l'itération max : %(i1)d
   - Convergence irrégulière & erreur >   %(r1)f
   - Diminuer la taille d'incrément.
"""
    ),
    53: _(
        """
 Erreur de programmation MULT_FRONT
   * Sur connexité des Lagrange Lambda1
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    54: _(
        """
     ==== Type de maille Aster / Type de maille GMSH ====
"""
    ),
    55: _(
        """
    %(i1)d  éléments %(k1)s découpés en %(i2)d  éléments %(k2)s a %(i3)d noeuds
"""
    ),
    56: _(
        """
Commande FACTORISER :
   Il faut toujours utiliser le mot clé 'reuse'.

   Sauf si la méthode est 'GCPC' et le préconditionnement 'LDLT_INC',
   auquel cas, il est interdit d'utiliser 'reuse'
"""
    ),
    57: _(
        """
    Le préconditionnement d'une matrice assemblée complexe n'est pas permis.
"""
    ),
    58: _(
        """
    La masse du modèle est nulle.
    On ne peut donc pas normer par rapport à la masse.
"""
    ),
    59: _(
        """
 MULT_FRONT: Erreur dans la renumérotation
   - Le super noeud : %(i1)d
   - devrait être le fils de   %(i2)d

 Risques & conseils :
   - Vous devriez rencontrer des problèmes lors de la factorisation.
   - Essayez un autre algorithme pour la renumérotation : 'MD', 'MDA', ...
"""
    ),
    60: _(
        """
    La variante 'QZ_QR' de la méthode 'QZ' ne fonctionne qu'avec une matrice %(k1)s symétrique réelle
    et %(k2)s symétrique réelle définie positive. Donc elle ne traite pas les problèmes de flambement,
    les Lagrange issus de AFFE_CHAR_MECA, des matrices complexes ou non symétriques,
    ni les problèmes modaux quadratiques.
"""
    ),
    61: _(
        """
    Méthode 'QZ' : propriété spectrale non respectée sur la valeur propre numéro %(i1)d.
    Les relations |alpha| < ||A|| et |bêta| < ||B|| ne sont pas vérifiées :
          |alpha|=%(r1)f,  ||A||=%(r2)f
          |bêta| =%(r3)f,  ||B||=%(r4)f
"""
    ),
    62: _(
        """
    Méthode 'QZ' dans CALC_MODES : On trouve un nombre de valeurs propres
    %(i1)d différent du nombre de ddls physiques actifs %(i2)d !
"""
    ),
    63: _(
        """
    Méthode 'QZ' dans CALC_MODES + OPTION='BANDE': On trouve un nombre de
    valeurs propres %(i1)d différent du nombre de valeurs propres détectées
    dans la bande %(i2)d !
"""
    ),
    64: _(
        """
    La méthode de 'JACOBI' n'est pas utilisable pour un problème modal quadratique
    (présence d'une matrice %(k1)s).

    Conseil :
    Utiliser la méthode 'SORENSEN' ou 'TRI_DIAG'.
"""
    ),
    65: _(
        """
    L'option de calcul 'TOUT' (sous le mot-clé facteur %(k1)s)
    est licite seulement avec METHODE='QZ'.
"""
    ),
    66: _(
        """
    Méthode 'QZ' dans CALC_MODES : On souhaite un nombre de valeurs
    propres %(i1)d supérieur au nombre de valeurs propres détectées %(i2)d !
"""
    ),
    68: _(
        """
    Méthode 'QZ' dans CALC_MODES : erreur LAPACK %(i1)d !
"""
    ),
    69: _(
        """
    Au moins une des matrices est non symétrique.
    Pour l'instant, seules les méthodes 'SORENSEN' et 'QZ' peuvent traiter le cas de
    matrices non symétriques.

    Conseils :
    - Si le problème modal est de petite taille (quelques centaines de DDL),
      utiliser la méthode 'QZ'.
    - Sinon, utiliser 'SORENSEN'.
"""
    ),
    70: _(
        """
    Au moins une des matrices est non symétrique, et la matrice %(k1)s est complexe.
    Pour l'instant, ce cas n'a pas été développé dans le code.
"""
    ),
    71: _(
        """
    Cette fonctionnalité requiert un solveur linéaire permettant de
    détecter les éventuelles singularités des matrices.

    Conseil :
    changer de solveur linéaire : sous le mot-clé facteur SOLVEUR,
    utiliser 'MULT_FRONT' ou 'MUMPS'.
"""
    ),
    73: _(
        """
    On a besoin d'effectuer un calcul de déterminant.
    Pour l'instant seuls les solveurs linéaires directs 'MULT_FRONT', 'LDLT'
    et MUMPS (à partir de la version 4.10.0) peuvent effectuer ce type de calcul.

    Conseil :
    Choisir un de ces deux solveurs (mot-clé METHODE sous le mot-clé facteur SOLVEUR).
    De préférence, utiliser 'MUMPS' (à partir de la version 4.10.0) qui est souvent
    le plus efficace pour des gros problèmes et/ou des problèmes difficiles
    (X-FEM, incompressibilité, THM, ...).
"""
    ),
    74: _(
        """
    Vous utilisez une fonctionnalité qui nécessite de connaître le degré de singularité de matrices associées à
    des systèmes linéaires. Or, vous avez désactivé la détection de singularité avec le mot-clé NPREC.

    Conseils :
      - Relancez le calcul avec NPREC > 0 (par exemple 8) sous le mot-clé facteur SOLVEUR.
      - S'il vous est indispensable de désactiver la détection de singularité, essayez d'utiliser un autre solveur linéaire, comme MULT_FRONT par exemple.
"""
    ),
    75: _(
        """
    Le solveur modal n'a pas réussi à capturer tous les modes propres souhaités
    avec le niveau de convergence requis.

    Conseils :
    Pour améliorer la convergence des algorithmes modaux vous pouvez par exemple :
     - Diminuer le nombre de modes recherchés à chaque fois en découpant votre calcul modal en plusieurs bandes.
       Cela améliore aussi souvent grandement les performances des calculs.
     - Avec la méthode de 'SORENSEN', augmenter la taille de l'espace de projection (DIM_SOUS_ESPACE/COEF_DIM_ESPACE)
       ou jouer sur les paramètres numériques qui pilotent la convergence (PREC_SOREN et NMAX_ITER_SOREN).
     - Avec la méthode 'QZ', diminuer NMAX_FREQ ou changer de variante (mot-clé TYPE_QZ).
     Si vous voulez tout de même utiliser les modes ainsi calculés (à vos risques et périls),
     relancer le calcul en augmentant la valeur du seuil de convergence (mot-clé SEUIL)
     ou en utilisant l'option STOP_ERREUR='NON' sous le mot-clé facteur VERI_MODE.
"""
    ),
    76: _(
        """
   Solveur GCPC :
   la création du préconditionneur %(k1)s a échoué car on manque de mémoire.

   Conseil :
   augmenter la valeur du paramètre PCENT_PIVOT sous le mot-clé facteur SOLVEUR.
"""
    ),
    77: _(
        """
Conseils :
Si vous utilisez METHODE='SORENSEN' ou 'TRI_DIAG' ou 'JACOBI', vous pouvez améliorer cette norme :
 - Si la dimension de l'espace réduit est inférieure à (nombre de degrés de liberté actifs - 2), augmenter la valeur de
   COEF_DIM_ESPACE (la valeur par défaut est 4 pour 'TRI_DIAG' et 2 pour 'SORENSEN' et 'JACOBI').
 - Découper le calcul en plusieurs appels de manière à réduire le nombre de modes propres recherchés simultanément
   (%(k1)s ou taille de la BANDE).
"""
    ),
    78: _(
        """
Conseils :
Vous pouvez améliorer cette norme :
 - en augmentant les nombres d'itérations des algorithmes
    (paramètres NMAX_ITER_SEPARE/AJUSTE sous le mot-clé facteur %(k1)s
     et/ou paramètre NMAX_ITER sous le mot-clé facteur CALC_MODE),
 - en augmentant la précision requise
    (paramètres PREC_SEPARE/PREC_AJUSTE sous le mot-clé facteur %(k1)s
     et/ou paramètre PREC sous le mot-clé facteur CALC_MODE),
 - en changeant d'algorithme
    (mot-clé OPTION sous le mot-clé facteur %(k1)s
     et/ou mot-clé OPTION sous le mot-clé facteur CALC_MODE).
"""
    ),
    79: _(
        """
    On souhaite un nombre de valeurs propres NMAX_%(k1)s=%(i1)d
    supérieur au nombre de valeurs propres détectées NUM=%(i2)d !
"""
    ),
    80: _(
        """
    Pour poursuivre le calcul, on impose NMAX_%(k1)s=NUM.
"""
    ),
    81: _(
        """
    Ce problème peut être dû :
    - à un mauvais tri dans les valeurs propres complexes conjuguées.
      Contacter l'équipe de développement.

    - à une mauvaise convergence de la méthode.
      Regarder les paramètres permettant d'améliorer celle-ci.

    - à une action incomplète du SHIFT.
      En diminuant la valeur %(k1)s de l'option 'CENTRE'
      et en augmentant le nombre de valeurs propres retenues (NMAX_%(k1)s),
      on peut souvent capter tous les couples (lambda,conjugué(lambda)) souhaités.

    Sinon utiliser METHODE='QZ' pour les problèmes de petites tailles (<500 ddls).
"""
    ),
    82: _(
        """
L'option 'PLUS_GRANDE' n'est pas utilisable en présence d'une matrice d'amortissement,
d'une matrice de rigidité complexe, ou de matrices non symétriques.
"""
    ),
}
