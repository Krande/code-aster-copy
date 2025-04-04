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
        """On ne peut pas faire propager une interface. Seules les fissures (possédant un fond de fissure) peuvent être propagées."""
    ),
    3: _(
        """On ne peut pas définir du contact X-FEM sur un maillage contenant à la fois des mailles linéaires et quadratiques."""
    ),
    6: _("""Le nombre de fissures est limité à %(i1)d, or vous en avez définies %(i2)d."""),
    7: _(
        """La fissure %(k1)s est de type CZM. On ne peut pas appliquer de chargement de pression ou de condition d'échange thermique sur ses lèvres."""
    ),
    10: _(
        """
  -> Toutes les fissures ne sont pas rattachées au même maillage.
     La fissure %(k1)s est rattachée au maillage %(k2)s alors que
     la fissure %(k3)s est rattachée au maillage %(k4)s.
"""
    ),
    15: _(
        """
  -> Point de FOND_FISS sans maille de surface rattachée.
  -> Risque & Conseil:
     Veuillez revoir la définition des level-sets.
"""
    ),
    17: _(
        """
  -> Segment de FOND_FISS sans maille de surface rattachée
  -> Risque & Conseil:
     Veuillez revoir la définition des level-sets.
"""
    ),
    18: _(
        """
  -> Le mot-clé CRITERE de PROPA_FISS est différent de 'ANGLE_IMPO' ou 'ANGLE_IMPO_BETA_GAMMA' et le tableau
     des facteurs d'intensité de contraintes de la fissure %(k1)s contient
     une colonne 'BETA' ou une colonne 'GAMMA'.
  -> Risque & Conseil:
     Les valeurs de l'angle de bifurcation notées dans ce tableau ne sont
     pas prises en compte. Si vous souhaitez imposer les valeurs de l'angle
     de bifurcation aux points du fonds de fissure, veuillez indiquer
     CRITERE='ANGLE_IMPO' ou 'ANGLE_IMPO_BETA_GAMMA'.
"""
    ),
    19: _(
        """
  -> Le mot-clé CRITERE de PROPA_FISS vaut 'ANGLE_IMPO' et le tableau
     des facteurs d'intensité de contraintes de la fissure %(k1)s ne contient
     pas de colonne 'BETA'.
  -> Risque & Conseil:
     Si vous souhaitez imposer les valeurs de l'angle de bifurcation aux points
     du fonds de fissure, veuillez indiquer CRITERE='ANGLE_IMPO' et ajouter
     une colonne 'BETA' au tableau manuellement ou si le modèle est en 3D,
     en utilisant l'option 'CALC_K_G' de la commande CALC_G_XFEM.
"""
    ),
    20: _(
        """
  -> En 3D, si METHODE_PROPA='MAILLAGE' dans PROPA_FISS il faut absolument une
     colonne 'ABSC_CURV' contenant les abscisses curvilignes des points du fond
     dans le tableau des facteurs d'intensité de contraintes.
  -> Risque & Conseil:
     Veuillez vérifier la présence de cette colonne.
"""
    ),
    30: _(
        """
     Le modèle de visualisation (mot-clé MODELE_VISU) utilisé est un
     modèle X-FEM.

     Risque & Conseil:
     Veuillez créer un modèle FEM sur un maillage de visualisation fissuré
     par l'utilisation des commandes POST_MAIL_XFEM puis AFFE_MODELE

"""
    ),
    50: _(
        """
  -> Le maillage utilisé pour la représentation des level-sets est 2D
     mais il contient des éléments 1D aussi.
  -> La méthode UPWIND sélectionnée dans PROPA_FISS peut gérer des
     grilles 2D définies seulement par des éléments QUAD4.
  -> Risque & Conseil:
     Veuillez donner un maillage défini seulement par des éléments
     QUAD4.
  """
    ),
    51: _(
        """
  -> Il n'y a aucune maille enrichie.
  -> Risque & Conseil:
     Veuillez vérifier les définitions des level-sets.
  """
    ),
    52: _(
        """
  -> Le maillage utilisé pour la représentation des level-sets est 3D
     mais il contient des éléments 2D et/ou 1D aussi.
  -> La méthode UPWIND sélectionnée dans PROPA_FISS peut gérer des
     grilles 3D définies seulement par des éléments HEXA8.
  -> Risque & Conseil:
     Veuillez donner un maillage défini seulement par des éléments
     HEXA8.
  """
    ),
    53: _(
        """
  -> Dans le maillage utilisé pour la représentation des level-sets,
     il y a des éléments qui ne sont pas disponibles pour la méthode
     UPWIND (PROPA_FISS).
  -> Risque & Conseil:
     Veuillez vérifier le maillage et utiliser uniquement des éléments
     QUAD4 en 2D et HEXA8 en 3D.
  """
    ),
    54: _(
        """
  -> Il n'y a pas d'éléments disponibles pour la méthode UPWIND
     (PROPA_FISS) dans le maillage utilisé pour la représentation
     des level-sets.
  -> Risque & Conseil:
     Veuillez vérifier le maillage et utiliser uniquement des éléments
     QUAD4 en 2D et HEXA8 en 3D.
  """
    ),
    55: _(
        """
  -> Dans le maillage utilisé pour la représentation des level-sets
     (PROPA_FISS), il y a des arêtes qui ne sont pas orthogonales aux
     autres arêtes.
  -> Risque & Conseil:
     Risques de résultats faux.
     Veuillez vérifier que toutes les arêtes des éléments du maillage
     soient orthogonales entre elles.
  """
    ),
    57: _(
        """
  -> La définition de un ou plusieurs éléments du maillage utilisé pour
     la représentation des level-sets (PROPA_FISS) n'est pas correcte.
  -> Risque & Conseil:
     Il y a une arête avec une longueur nulle dans le maillage.
     Veuillez vérifier la définition des éléments du maillage (par
     exemple: un noeud est utilisé seulement une fois dans la définition
     d'un élément; il n'y a pas de noeuds doubles...)
  """
    ),
    58: _(
        """
  -> La dimension (2D ou 3D) du modèle physique et la dimension (2D ou
     3D) du modèle utilisé pour la grille auxiliaire ne sont pas égales.
  -> Risque & Conseil:
     Veuillez utiliser deux modèles avec la même dimension (les deux 2D
     ou les deux 3D).
  """
    ),
    59: _(
        """
Il y a au moins une maille qui porte des sous-éléments de peau et qui
ne borde pourtant aucune maille principale.
  """
    ),
    60: _(
        """
  -> L'opérande TEST_MAIL a été utilisée dans l'opérateur PROPA_FISS.
     La même vitesse d'avancée est utilisée pour tous les points du
     fond de fissure et l'angle de propagation est fixé égal à zéro.
  -> Risque & Conseil:
     L'avancée de la fissure n'est pas liée aux contraintes affectant
     la structure et donc les résultats de la propagation n'ont pas
     une signification physique.
     L'opérande TEST_MAIL doit être utilisé uniquement pour vérifier
     si le maillage est suffisamment raffiné pour la représentation
     des level-sets.
  """
    ),
    63: _(
        """
  -> La valeur de l'avancée DA_MAX utilisée est petite par rapport à la
     longueur de la plus petite arrête du maillage utilisé pour
     la représentation des level-sets:
     DA_MAX = %(r1)f
     Longueur minimale arrêt = %(r2)f
  -> Risque & Conseil:
     Risques de résultats faux. Veuillez vérifier les résultats en
     utilisant un maillage plus raffiné pour la représentation des
     level-sets.
  """
    ),
    64: _(
        """
  -> La valeur du RAYON est plus petite que la longueur de la plus petite
     arrête du maillage utilisé pour la représentation des level-sets:
     RAYON = %(r1)f
     LONGUEUR minimale arrêt = %(r2)f
  -> Le calcul du résidu local n'est pas possible.
  -> Risque & Conseil:
     Veuillez utiliser une valeur du RAYON plus grande.
  """
    ),
    70: _(
        """
  -> La macro-commande PROPA_FISS ne peut traiter qu'un seul instant de calcul.
  -> Risque & Conseil:
     Veuillez vérifier que les tableaux des facteurs d'intensité de contraintes
     donnés dans l'opérateur PROPA_FISS ne contiennent qu'un seul instant.
"""
    ),
    73: _(
        """
  -> L'option NB_POINT_FOND a été utilisé dans PROPA_FISS mais le
     modèle est 2D.
  -> Risque & Conseil:
     Cette option n'est utile qu'avec un modèle 3D.
     Ce mot-clé n'est pas pris en compte.
  """
    ),
    74: _(
        """
  -> Aucune fissure du modèle ne se propage.
  -> Risque & Conseil:
     Veuillez vérifier les conditions du chargement du modèle et les
     constantes de la loi de propagation données à PROPA_FISS.
  """
    ),
    75: _(
        """
  -> Une valeur de la liste de NB_POINT_FOND ne correspond pas au nombre de
     lignes du tableau des facteurs d'intensité de contraintes pour
     le fond %(i1)d de la fissure %(k1)s.
  -> Risque & Conseil:
     Veuillez vérifier que la liste NB_POINT_FOND donnée dans PROPA_FISS
     soit identique à celle utilisée pour construire le tableau.
  """
    ),
    78: _(
        """
  -> L'option NB_POINT_FOND a été utilisée dans PROPA_FISS
     mais le nombre de valeurs données n'est pas égale au nombre total
     des morceaux des fissures dans le modèle.

  -> Conseil:
     Veuillez vérifier que l'option NB_POINT_FOND a été utilisée
     correctement dans PROPA_FISS et que les valeurs données pour
     chaque fissure sont correctes.
  """
    ),
    80: _(
        """
  -> Le nombre des valeurs dans un des tableaux des facteurs
     d'intensité de contraintes est supérieur au nombre des
     points du fond de la fissure correspondante.
  -> Risque & Conseil:
     Veuillez vérifier que les tableaux donnés par l'opérateur
     PROPA_FISS sont corrects. Si NB_POINT_FOND a été utilisé, veuillez
     vérifier aussi que la liste donnée pour chaque fissure est correcte.
  """
    ),
    81: _(
        """
  -> Les valeurs de COEF_MULT_MAXI et COEF_MULT_MINI de COMP_LINE sont
     égales à zéro.
  -> Risque & Conseil:
     Au moins une des deux valeurs doit être différente de zéro pour
     avoir un cycle de fatigue. Veuillez vérifier les valeurs données.
  """
    ),
    85: _(
        """
   Les propriétés matériaux dépendent de la température. La température en fond
   de fissure n'étant pas connue, le calcul se poursuit en prenant la température
   de référence du matériau (TEMP = %(r1)f).
"""
    ),
    86: _(
        """
 -> Le maillage/la grille sur lequel/laquelle vous voulez créer le groupe
    n'est pas associé/associée à la fissure donnée.

 -> Risque & Conseil:
    Veuillez vérifier d'avoir spécifié le bon maillage/grille et/ou
    la bonne fissure.
"""
    ),
    87: _(
        """
  -> L'opérande TEST_MAIL a été utilisé dans l'opérateur PROPA_FISS.
  -> Cet opérande n'a de sens que pour un modèle 3D.
  -> Risque & Conseil:
     Ne pas utiliser TEST_MAIL pour un modèle 2D.
  """
    ),
    88: _(
        """
  -> La valeur du rayon du tore de localisation de la zone de mise à
     jour est supérieure à la valeur limite. Cette dernière est
     déterminée par la valeur du rayon du tore utilisée à la propagation
     précédente et la valeur de l'avancée de la fissure (DA_MAX) imposée
     à la propagation courante.

     Rayon actuel = %(r1)f
     Rayon limite = %(r2)f

  -> Risque & Conseil:
     Risques de résultats faux si la fissure ne propage pas en mode I.

     Pour éviter ce risque, vous pouvez utiliser la même avancée de la
     fissure (DA_MAX) et le même rayon (RAYON) que ceux qui ont été
     utilisés à la propagation précédente.

     Si vous ne pouvez pas utiliser les même valeurs, vous pouvez
     choisir une des solutions suivantes:
     - donner une valeur de RAYON_TORE inférieure à la valeur limite
       pour la propagation courante
     - utiliser une valeur de RAYON_TORE plus grande pour les
       propagations précédentes
     - augmenter l'avancée de la fissure (DA_MAX) à la propagation
       courante

     Sinon, même si fortement déconseillé, vous pouvez choisir de ne pas
     utiliser la localisation de la zone de mise à jour
     (ZONE_MAJ='TOUT').
  """
    ),
    89: _(
        """
  -> La fissure à propager n'existe pas dans le modèle:
     FISS_ACTUELLE = %(k1)s
     MODELE        = %(k2)s
  -> Conseil:
     Veuillez vérifier que la fissure et le modèle sont correctement
     définis.
  """
    ),
    91: _(
        """
  -> Le nouveau fond de fissure n'est pas très régulier. Cela signifie
     que le maillage ou la grille auxiliaire utilisés pour la
     représentation de la fissure par level-sets ne sont pas
     suffisamment raffinés pour bien décrire la forme du fond de la
     fissure utilisée.
  -> Risque & Conseil:
     Risques de résultats faux en utilisant le maillage ou la grille
     auxiliaire testés. Veuillez utiliser un maillage ou une grille
     auxiliaire plus raffinés.
  """
    ),
    92: _(
        """
Vous avez demandé la création d'un groupe de noeuds dans un tore construit autour du fond de la fissure suivante:
FISSURE = %(k1)s

Toutefois cette fissure a été calculée par PROPA_FISS en utilisant la localisation du domaine.
Dans ce cas le groupe de noeuds doit être forcement défini en utilisant le tore déjà utilisé par PROPA_FISS.

Le groupe de noeuds sera créé en utilisant le domaine de localisation de la fissure (option TYPE_GROUP='ZONE_MAJ')."""
    ),
    93: _(
        """
  -> Aucune fissure n'est définie sur le modèle spécifié:
     MODELE = %(k1)s
  -> Risque & Conseil:
     Veuillez définir une fissure sur le modèle ci-dessus en utilisant
     les opérateurs DEFI_FISS_XFEM et MODI_MODELE_XFEM avant
     l'utilisation de PROPA_FISS.
  """
    ),
    94: _(
        """
  -> L'avancée donnée (DA_MAX) pour la propagation courante est
     inférieure à la valeur minimale conseillée.

     DA_MAX donnée                     = %(r1)f
     Avancée maximale fissure courante = %(r2)f
     DA_MAX minimal conseillé          = %(r3)f

  -> Risque & Conseil:
     Risque de résultats faux. Dans le cas de propagation 3D en mode
     mixte, on conseille en général d'utiliser une avancée de fissure
     supérieure à celle minimale écrite ci-dessus, même si des bonnes
     résultats peuvent être obtenus en utilisant une avancée inférieure.

     La valeur minimale de DA_MAX dépende de la valeur de l'opérande
     RAYON et de l'angle de propagation de la fissure. Dans le cas où la
     valeur DA_MAX donnée ne peut pas être changée, sa valeur minimale
     conseillée peut être diminuée en agissant sur la valeur de RAYON,
     c'est-à-dire en utilisant une valeur de RAYON plus petite. Cela
     influence l'opérateur CALC_G_XFEM aussi et normalement est faisable en
     utilisant un maillage plus raffiné.
  """
    ),
    95: _(
        """
  -> Le modèle grille donné est défini sur un maillage (%(k1)s)
     et pas sur une grille.

  -> Risque & Conseil:
     Veuillez donner un modèle grille défini sur une grille. Cette
     grille doit être définie par DEFI_GRILLE à partir d'un maillage.
  """
    ),
    96: _(
        """
 -> Le maillage sur lequel vous voulez créer le groupe n'est pas associé à
    la fissure donnée.

 -> Les maillages suivants sont associés à cette fissure:
      maillage physique = %(k1)s
      maillage grille   = %(k2)s

 -> Risque & Conseil:
    Veuillez vérifier d'avoir spécifié le bon maillage et/ou la bonne fissure.
"""
    ),
    97: _(
        """
  -> La localisation de la zone de mise à jour a été utilisé pour la
     détermination de la configuration actuelle des fissures du modèle.
     Par contre, pour la propagation courante, la localisation n'a pas
     été activée.
  -> Risque & Conseil:
     Veuillez utiliser la localisation de la zone de mise à jour
     (ZONE_MAJ='TORE') pour la propagation courante aussi.
  """
    ),
    98: _(
        """
  -> Aucune grille auxiliaire n'est utilisée pour la représentation de
     la fissure donnée.
  -> Risque & Conseil:
     Veuillez vérifier que vous avez demandé les level-sets de la bonne
     fissure.
  """
    ),
    99: _(
        """
  -> La valeur du rayon du tore de localisation de la zone de mise à
     jour est plus petite que celle qui est nécessaire pour la bonne
     mise à jour des level-sets.
     Rayon à utiliser = %(r1)f
     Rayon minimal    = %(r2)f
  -> Risque & Conseil:

  -> Si vous avez utilisé l'opérande RAYON_TORE, veuillez augmenter la
     valeur donné ou diminuer la valeur de DA_MAX ou RAYON.

  -> Si vous n'avez pas utilisé l'opérande RAYON_TORE, cette erreur
     signifie que l'estimation automatique faite par l'opérateur
     PROPA_FISS ne marche pas bien pour la propagation courante. Elle
     peut être utilisée dans les cas où les valeurs de RAYON et DA_MAX
     ne changent pas entre deux propagations successives et la taille
     des éléments dans la zone de propagation est presque constante.
     Veuillez donc donner explicitement une valeur du rayon en utilisant
     l'opérande RAYON_TORE.
     Vous pouvez calculer une estimation de cette valeur en utilisant la
     formule suivante:

     RAYON_TORE=RAYON_max+DA_MAX_max+2*h_max

     où RAYON_max et DA_MAX_max sont les valeurs maximales des opérandes
     RAYON et DA_MAX qu'on va utiliser et h_max est la valeur de la plus
     grande arête des éléments du maillage ou de la grille auxiliaire
     dans la zone de propagation.
  """
    ),
}
