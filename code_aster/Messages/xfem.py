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
  -> Lors de la définition du champ de matériaux %(k1)s, vous avez renseigné le mot-clé
     EVOL de AFFE_MATERIAU / AFFE_VARC avec un résultat thermique X-FEM.
     Cette fonctionnalité n'est pas compatible avec l'opérateur POST_K1_K2_K3 lorsque le
     matériau dépend de la température.
"""
    ),
    2: _(
        """
  -> Le calcul de la distance d'un noeud à l'ellipse n'a pas convergé
     avec le nombre d'itérations maximal fixé (10). Cela est dû à une
     ellipse très allongée.
  -> Conseil:
     Contacter les développeurs.
     Dans la mesure du possible, définissez une ellipse moins allongée.
"""
    ),
    3: _(
        """
  -> Le modèle %(k1)s est incompatible avec la méthode X-FEM.
  -> Risque & Conseil:
     Vérifiez qu'il a bien été créé par l'opérateur MODI_MODELE_XFEM.
"""
    ),
    5: _(
        """
  -> Attention, vous avez défini un enrichissement géométrique sur %(i1)d
     couches d'éléments autour du fond de fissure.
  -> Risque :
     Au delà de 7 couches, il y a des risques de pivots nuls lors de la
     résolution dans STAT_NON_LINE.
  -> Conseils :
     Pour éviter ces risques de pivots nuls, il est conseillé de ne pas
     dépasser NB_COUCHES = 7.
     Vous pouvez aussi laisser NB_COUCHES = %(i1)d, mais il pourra s'avérer
     nécessaire d'augmenter le nombre maximales de décimales perdues dans
     STAT_NON_LINE (mot-clé NPREC de SOLVEUR pour les méthodes LDLT ou MULT_FRONT)
"""
    ),
    6: _(
        """
  -> Le rayon d'enrichissement RAYON_ENRI doit être un réel strictement
     supérieur à 0.
"""
    ),
    7: _(
        """
     Il y a %(i1)d mailles %(k1)s
"""
    ),
    8: _(
        """
     Le nombre de %(k1)s X-FEM est limité à 2.10E9
     Risque & Conseil:
     Veuillez réduire la taille du maillage.
"""
    ),
    9: _(
        """
     Le groupe de mailles donné pour définir la fissure contient des
     mailles qui ne sont pas connectées aux autres. Cela empêche
     d'orienter correctement la normale à la surface de la fissure.

     Risque & Conseil:
       Veuillez vérifier que les mailles données en entrée sont toutes
       connectées entre elles, c'est-à-dire qu'elle forme un groupe de
       mailles contiguës.
       Sinon, il faut définir une fissure pour chacun des groupes
       non connexes.
"""
    ),
    10: _(
        """
     La direction du champ thêta n'a pas été donnée. La direction automatique
     est une direction variable, basée sur le gradient de la level-set tangente.
"""
    ),
    11: _(
        """
  -> On a trouvé plus de 2 points de fond de fissure, ce qui est impossible en 2D.
  -> Risque & Conseil:
     Cela est normalement causé par une mauvaise définition des level-sets.

     Si les level-sets ont été définies par DEFI_FISS_XFEM, veuillez revoir leur
     définition.

     Si les level-sets ont été calculées par PROPA_FISS, vous pouvez essayer
     d'utiliser un maillage plus raffiné dans toute la zone de propagation ou bien
     une grille auxiliaire.
     Si vous avez utilisé la méthode simplexe avec restriction de la zone de mise à
     jour des level-sets (ZONA_MAJ='TORE'), vous pouvez spécifier un rayon plus élevé
     de celui qui a été utilisé (écrit dans le fichier .mess) en utilisant l'opérande
     RAYON_TORE ou vous pouvez déactiver la restriction de la zone de mise à jour
     (ZONE_MAJ='TOUT').
     Sinon vous pouvez changer la méthode utilisée par PROPA_FISS (opérande
     METHODE_PROPA).
"""
    ),
    12: _(
        """
  Le gradient de la level-set tangente est nul au noeud %(k1)s.
  Ceci est certainement du à un point singulier dans la définition de la level-set.
  Il vaut veiller à ce que ce point singulier ne soit pas inclus dans la couronne
  d'intégration du champ thêta.
  Conseil : réduisez la taille de la couronne du champ thêta
"""
    ),
    13: _(
        """
     Dans le modèle, des mailles SEG2 ou SEG3 possèdent des noeuds enrichis par X-FEM.
     Ceci n'est pas encore possible en 3D.
     Conseils : si ces mailles sont importantes pour le calcul (charge linéique...), il faut
     les mettre loin de de la fissure.
     Si ces mailles ne servent pas pour le calcul, il vaut mieux ne pas les affecter dans le modèle,
     ou bien les supprimer du maillage.
"""
    ),
    14: _(
        """
Le calcul de la norme L2 de la pression de contact sur une fissure XFEM n'est pas disponible.
"""
    ),
    15: _(
        """
  -> Cette option n'a pas encore été programmée.
  -> Risque & Conseil:
     Veuillez utiliser un autre chargement (en pression) ou contacter votre
     correspondant.
"""
    ),
    16: _(
        """
  -> Il n'y a aucun élément enrichi.
"""
    ),
    17: _(
        """
     il ne faut qu'un mot-clé parmi RAYON_ENRI et NB_COUCHES.
"""
    ),
    18: _(
        """
     Dimension de l'espace incorrecte.
     Le modèle doit être 2D ou 3D et ne pas comporter de sous-structures.
"""
    ),
    19: _(
        """
     Il y a %(i1)d mailles dans la zone fissure.
"""
    ),
    20: _(
        """
     Vous avez défini plus d'un fond fermé lors d'un appel à DEFI_FISS_XFEM.
     Conseil:
     Veuillez définir chacun des fonds dans des DEFI_FISS_XFEM différents.
"""
    ),
    21: _(
        """
     Vous avez défini au moins un fond ouvert et un fond fermé lors d'un appel à
     DEFI_FISS_XFEM.
     Conseil:
     Veuillez définir chacun des fonds dans des DEFI_FISS_XFEM différents.
"""
    ),
    22: _(
        """
     Une relation cinématique est imposée au niveau du noeud %(k1)s sur lequel une
     fissure passe et cette relation est imposée de part et d'autre de la fissure.
     Deux relations sont donc imposées, bloquant ainsi le ddl Heaviside associé.
     Si vous avez également imposé une relation cinématique sur le degré de liberté
     Heaviside associé du noeud  %(k1)s, la matrice sera non factorisable.
     Conseil :
     Vous pouvez exclure le noeud %(k1)s du groupe affecté par
     cette deuxième relation en utilisant par exemple l'opérande DIFFE.
"""
    ),
    23: _(
        """
     Erreur dans le choix de la méthode de calcul des level-sets.
     Vous souhaitez définir une %(k1)s.
     Or la forme que vous avez sélectionnée < %(k2)s >
     correspond à une %(k3)s.
     Conseil :
     Sélectionnez une forme de %(k1)s.
"""
    ),
    24: _(
        """
     Erreur dans le choix de la méthode de calcul des level-sets.
     Vous souhaitez définir une fissure.
     Pour cela il est nécessaire de définir 2 level-sets : LT et LN.
     Conseil :
     Veuillez renseignez %(k1)s.
"""
    ),
    25: _(
        """
     Erreur dans le choix de la méthode de calcul des level-sets.
     Vous souhaitez définir une interface.
     Pour cela il ne faut est pas définir la level-set normale LT.
     %(k1)s ne sera pas considéré.
     Conseil :
     Pour ne plus obtenir ce message, ne renseignez pas %(k1)s.
"""
    ),
    26: _(
        """
     Numéros des mailles de la zone fissure.
"""
    ),
    27: _(
        """
     L'utilisation de la grille auxiliaire avec la méthode SIMPLEXE n'est pas disponible.
     Utiliser la méthode UPWIND ou effectuer le calcul sur le maillage physique.
"""
    ),
    28: _(
        """
     Les flux de fluide dans les fractures ne sont autorisés que pour les modélisations HM-XFEM.
"""
    ),
    29: _(
        """
     Nombre de mailles contenant le fond de fissure : %(i1)d
"""
    ),
    30: _(
        """
     Nombre de mailles de type Heaviside : %(i1)d
"""
    ),
    31: _(
        """
     Nombre de mailles de type Crack-tip : %(i1)d
"""
    ),
    32: _(
        """
     Nombre de mailles de type Heaviside Crack-tip : %(i1)d
"""
    ),
    33: _(
        """
     Nombre de points du fond de fissure : %(i1)d
"""
    ),
    34: _(
        """
     Nombre de fonds de fissure : %(i1)d
"""
    ),
    35: _(
        """
     Coordonnées des points des fonds de fissure
"""
    ),
    36: _(
        """
     fond de fissure : %(i1)d
"""
    ),
    37: _(
        """
     Nombre de level-sets réajustées : %(i1)d
"""
    ),
    39: _(
        """
     Erreur utilisateur : incohérence entre les mots-clés FISSURE et MODELE_IN.
     Il faut que les (ou la) fissure sous le mot-clé FISSURE soient toutes définies à
     partir du même maillage.
     Or :
     - la fissure %(k1)s est définie à partir du maillage %(k2)s
     - le modèle renseigné sous MODELE_IN est défini à partir du maillage %(k3)s.
     Conseil :
     Veuillez revoir la définition de la fissure %(k1)s ou bien changer MODELE_IN.
"""
    ),
    40: _(
        """
      La maille %(k1)s doit être enrichie avec plus de 4 fonctions Heaviside,
      le multi-Heaviside est limité à 4 fonctions Heaviside, un noeud du maillage
      ne doit pas être connecté à plus de 4 fissures.
      Pour ne pas activer le multi-Heaviside, les fissures doivent être séparées de 2 mailles
      minimum. Veuillez raffiner le maillage entre les fissures (ou écarter les fissures).
"""
    ),
    41: _(
        """
      La maille %(k1)s est quadratique et elle est connectée à 2 fissures,
      le multi-Heaviside n'a été généralisé en quadratique que pour les modélisations HM.
      Pour ne pas activer le multi-Heaviside, les fissures doivent être séparées de 2 mailles
      minimum. Veuillez raffiner le maillage entre les fissures (ou écarter les fissures).
"""
    ),
    42: _(
        """
      La table de facteurs d'intensité des contraintes donnée en entré pour la fissure
      %(k1)s ne contient pas le bon numéro de fonds de fissure et/ou de points.

      Veuillez faire les vérifications suivantes:
      - si la fissure %(k1)s  est formée par un seul fond, veuillez vérifier d'avoir donné
        la bonne table
      - si la fissure %(k1)s est formée par plusieurs fonds, veuillez vérifier que la
        table donnée contient les valeurs des facteurs d'intensité des contraintes de
        chaque fond (voir colonne NUME_FOND de la table)
"""
    ),
    43: _(
        """
      Le contact autre que P1P1 est actif et la maille %(k1)s est connectée à 2 fissures,
      à l'exception des modélisations HM-XFEM, le multi-Heaviside ne peut pas être pris
      en compte si le contact autre que P1P1 est utilisé.
      Pour ne pas activer le multi-Heaviside, les fissures doivent être séparées de 2 mailles
      minimum. Veuillez raffiner le maillage entre les fissures (ou écarter les fissures).
"""
    ),
    44: _(
        """
      La maille %(k1)s est connectée à 2 fissures, or il s'agit d'une maille possédant des
      enrichissements de fond de fissure.
      Le multi-Heaviside ne peut pas être pris en compte en fond de fissure.
      Pour ne pas activer le multi-Heaviside, les fissures doivent être séparées de 2 mailles
      minimum. Veuillez raffiner le maillage entre les fissures (ou écarter les fissures).
"""
    ),
    45: _(
        """
      Jonction X-FEM et contact

      Une facette de contact XFEM doit être redécoupée. Ceci n'est pas implémenté pour l'instant.
      Les efforts de contact ne seront pas prise en compte sur cette facette.
"""
    ),
    46: _(
        """
      Les fissures sont mal ordonnées dans le mot clé FISSURE de MODI_MODELE_XFEM

      L'ordre dans lequel sont définis les fissures avec l'utilisation du mot clé JONCTION impose que %(k1)s doit
      être donné après %(k2)s. Veiller permuter %(k1)s et %(k2)s dans le mot clé FISSURE de MODI_MODELE_XFEM.
"""
    ),
    47: _(
        """
      La fissure %(k1)s est déjà attaché à la fissure %(k2)s, on ne peut pas l'attacher à %(k3)s.

      Il est possible d'attacher globalement %(k1)s à la fois à %(k2)s et %(k3)s,
      mais il ne faut pas qu'un élément soit connecté à la fois à %(k2)s, %(k3)s et %(k1)s.

      Pour résoudre ce problème, soit il faut écarter (ou raffiner le maillage entre) les fissures %(k2)s et %(k3)s.
      Soit il faut lier la fissure %(k3)s à la fissure %(k2)s en ajoutant une ligne du type
      JONCTION=_F(FISSURE=%(k2)s,POINT=...) lorsqu'on appelle DEFI_FISS_XFEM pour définir %(k3)s.
"""
    ),
    48: _(
        """
      Les flux correspondant a PRE2 et TEMP sont interdits pour les
      éléments HM-XFEM.
"""
    ),
    50: _(
        """
     La méthode X-FEM n'est pas disponible avec %(k1)s.
"""
    ),
    51: _(
        """
     La maille %(k1)s possède %(i1)d points de fond de fissure de coordonnées :
"""
    ),
    52: _(
        """
     Une ou des mailles contenant plus de 2 points du fond de fissure ont été détectées.
     Le fond ne peut pas être orienté sous cette condition. Il n'est donc pas possible
     de calculer les abscisses curvilignes du fond et de détecter les fonds multiples.
     Par conséquent le post-traitement avec la commande CALC_G_XFEM n'est pas possible.
"""
    ),
    54: _(
        """
     Attention, deux jonctions de fissure sont présentes dans le même élément. Les facettes
     de contact récupérées dans cet élément risquent de ne pas être conformes à chaque jonction.
  -> Conseil :
     Raffinez le maillage afin de n'avoir qu'une seule jonction de fissure par élément.
"""
    ),
    55: _(
        """
     La découpe des facettes de contact XFEM à rencontrer un problème d'espace mémoire
     lors de l'écriture du mode local %(k1)s.
     La longueur du mode local dans le catalogue est %(i1)d, pourtant on a calculé qu'il
     faudrait %(i2)d.
     Le calcul s'arrête pour prévenir un dépassement mémoire.
  -> Conseils:
     - Vous pouvez augmenter le dimensionnement du mode local %(k1)s.
     - Veuillez aussi contacter l'équipe de développement pour reporter la configuration de
     coupe pour l'amélioration du catalogue.
"""
    ),
    56: _(
        """
  -> L'élément X-FEM de type %(k1)s porté par la maille %(k2)s ne
     possède pas suffisamment de points de Gauss pour stocker les
     informations aux points de Gauss de ses sous-éléments. Le nombre
     de points de Gauss de la famille XFEM pour cet élément est %(i1)d,
     alors que le nombre de points de Gauss nécessaire est %(i2)d.
  -> Conseil:
     Vous pouvez essayer de raffiner le maillage afin d'obtenir une
     découpe mettant en jeu moins de sous-éléments.
"""
    ),
    57: _(
        """
  -> La fissure (ou l'interface) définie dans DEFI_FISS_XFEM ne coupe aucune des mailles
     du maillage. La fissure (ou l'interface) se trouve donc en dehors de la structure ou
     bien coïncide avec un bord de la structure. Cela est interdit. Il y a probablement
     une erreur de mise en données.
  -> Conseil :
     Vérifier la cohérence entre la définition de la géométrie de la fissure et le maillage.
"""
    ),
    58: _(
        """
  -> Aucun point du fond de fissure n'a été trouvé !
     Cela signifie que le fond de fissure se trouve en dehors de la structure.

  -> Risque & Conseil :
     - Si vous souhaitiez définir une interface, il faut choisir TYPE_DISCONTINUITE = 'INTERFACE'
        pour ne plus avoir ce message.
     -  Si vous souhaitiez définir une fissure, il doit y avoir une erreur lors de la définition
        de la level-set tangente Vérifier la définition des level-sets.
     -  S'il s'agit d'un calcul de propagation, cela signifie que la fissure débouche et traverse
        entièrement la structure (risque de pivot nul dans le prochain calcul mécanique pour cause
        de modes de corps rigides non bloqués).
"""
    ),
    59: _(
        """
     Ne pas utiliser le mot-clef RAYON_ENRI lorsque le fond de fissure
     est en dehors de la structure.
"""
    ),
    60: _(
        """
     -> Pour le fond fermé, un point supplémentaire du fond a été ajouté.

"""
    ),
    63: _(
        """
      -> ---Éléments XFEM quadratiques ---
    On a effectué un ajustement géométrique de la fissure, car les arêtes du maillage sain sont
    coupées plusieurs fois par l'isovaleur zéro de la level-set.
    Le nombre d'ajustements effectués est : %(i1)d.
    Cette correction impacte légèrement la localisation de la fissure.
  -> Conseil :
    Veuillez vérifier en post-traitement grâce à la commande POST_MAIL_XFEM que la nouvelle
    géométrie de la fissure respecte la géométrie imposée dans la commande DEFI_FISS_XFEM.
"""
    ),
    64: _(
        """
  -> ---Éléments XFEM quadratiques 2D---
     Le calcul ne peut aboutir pour l'une des raisons suivante :
     - les calculs des coordonnées des points d'intersection entre un élément et la fissure
       se sont mal déroulés
     - l'élément ne peut être découpé selon la configuration de fissure qui le traverse
"""
    ),
    65: _(
        """
  -> ---Éléments XFEM quadratiques 2D---
     Le calcul d'abscisse curviligne sur une arête quadratique n'est pas encore supporté pour
     arête courbe.
"""
    ),
    66: _(
        """
  -> ---Éléments XFEM quadratiques 2D---
     Le calcul d'abscisse curviligne sur une arête quadratique ne peut aboutir pour l'une des
     raisons suivante :
     - les trois points qui définissent l'arête quadratique sont identiques
     - l'arête est "trop" arrondie
"""
    ),
    67: _(
        """
  -> ---Éléments XFEM quadratiques 2D---
     Newton : nombre d'itérations maximal atteint
"""
    ),
    68: _(
        """
  -> Aucune grille n'est associée à la fissure donnée par FISS_GRILLE.

  -> Risque & Conseil:
     Veuillez donner une fissure avec une grille associée.

"""
    ),
    69: _(
        """
  -> La fissure à propager a été définie par DEFI_FISS_XFEM en donnant directement les deux
     champs level-sets (mots-clés CHAMP_NO_LSN et CHAMP_NO_LST).
     Aucune grille auxiliaire n'a été associée à cette fissure.

  -> Risque & Conseil:
     Dans le cas où les deux champs level-sets ont été obtenus d'une fissure propagée par
     PROPA_FISS, les informations sur la localisation du domaine (mot-clé ZONE_MAJ) et sur
     l'utilisation d'une grille auxiliaire ont été perdues, ce qui fait que le calcul de la
     propagation de la fissure pourrait donner des résultats faux.

     Vous pouvez ignorer cette alarme seulement si les deux champs level-sets donnés dans
     DEFI_FISS_XFEM:
     - n'ont pas été calculés par PROPA_FISS, c'est-à-dire qu'ils n'ont pas été extraits
       d'une fissure propagée par PROPA_FISS
     - ont été extraits d'une fissure propagée par PROPA_FISS sans grille auxiliaire associée
       et la mise à jours des level-sets a été faite sur tous les noeuds du maillage
       (mot-clé ZONE_MAJ='TOUT' dans PROPA_FISS)

     Dans tous les autres cas, pour éviter des résultats faux, il faut absolument associer
     une grille auxiliaire à la fissure à propager, grille héritée de la fissure de laquelle
     les deux level-sets ont été extraites (mot-clé FISS_GRILLE de DEFI_FISS_XFEM).

"""
    ),
    70: _(
        """
  -> Un élément convexe a été détecté sur le fond de la fissure %(k1)s. PROPA_FISS ne sait
     pas traiter ce cas. Le calcul de propagation ne peut se poursuivre.
  -> Conseil:
     Veuillez revoir la définition de votre fissure.

"""
    ),
    71: _(
        """
     La jonction de fissures est une fonctionnalité disponible uniquement pour les
     modélisations mécaniques et hydro-mécaniques. Or le modèle %(k1)s est un modèle
     thermique.
  -> Conseil:
     Revoyez la définition de votre modèle, ou celle de la fissure (ou des fissures).
"""
    ),
    72: _(
        """
  -> Vous utilisez le mot-clé FISSURE, or le modèle %(k1)s que vous avez renseigné
     pour le mot clé MODELE n'est pas un modèle X-FEM.
  -> Conseil:
     Veuillez utiliser un autre mot-clé ou revoyez la définition de votre modèle.
"""
    ),
    73: _(
        """
  -> Vous avez renseigné %(k1)s pour le mot-clé FISSURE, or cette fissure est absente
     du modèle %(k2)s que vous avez renseigné pour le mot-clé MODELE.
  -> Conseil:
     Assurez vous de renseigner pour le mot-clé FISSURE une liste de fissures
     présentes dans le modèle ou revoyez la définition de votre modèle.
"""
    ),
    74: _(
        """
     Nombre de points du fond de fissure sur la grille : %(i1)d
"""
    ),
    75: _(
        """
     Coordonnées des points du fond de fissure sur la grille
"""
    ),
    77: _(
        """
  -> Il y a éventuellement des créations de mailles
     supplémentaires de type POI1 lorsque des affectations sont faites sur des nœuds ou des groupes de
     noeuds. Ces mailles ne sont pas accessibles à l’utilisateur. Ceci a crée une maille tardive.

  -> Risque & Conseil:
     Veuillez voir le document [U4.41.01], section 4:  Il est fortement conseillé
     d’utiliser CREA_MAILLAGE [U4.23.02] pour créer des mailles POI1 utilisables dans le fichier de
     commande (pour STAT_NON_LINE par exemple).

"""
    ),
    78: _(
        """
  -> Erreur, pour une modélisation %(k1)s, l'introduction de fissures dans le modèle n'est pas possible.

"""
    ),
    79: _(
        """
  -> Vous devez renseigner le mot-clé MODELE_IN avec un modèle sain (produit par
     l'opérateur AFFE_MODELE), or vous avez renseigné ce mot-clé avec %(k1)s
     qui est un modèle X-FEM (produit par l'opérateur MODI_MODELE_XFEM)
  -> Conseil:
     Revoyez la définition de ce modèle
"""
    ),
    80: _(
        """
  -> En présence du mot-clé MODELE_THER, vous devez renseigner le mot-clé MODELE_IN
     avec un modèle mécanique. Or le modèle %(k1)s n'est pas un modèle mécanique.
  -> Conseil:
     Revoyez la définition de ce modèle.
"""
    ),
    81: _(
        """
  -> Vous devez renseigner le mot-clé MODELE_THER avec un modèle X-FEM (produit par
     l'opérateur MODI_MODELE_XFEM), or vous avez renseigné ce mot-clé avec %(k1)s
     qui est un modèle sain (produit par l'opérateur AFFE_MODELE)
  -> Conseil:
     Revoyez la définition de ce modèle
"""
    ),
    82: _(
        """
  -> Vous devez renseigner le mot-clé MODELE_THER avec un modèle thermique.
     Or le modèle %(k1)s n'est pas un modèle thermique.
  -> Conseil:
     Revoyez la définition de ce modèle
"""
    ),
    83: _(
        """
  -> Les modèles %(k1)s et %(k2)s doivent nécessairement avoir été créés à partir
     du même maillage. Or %(k1)s et %(k2)s ont respectivement été définis à partir
     des maillages %(k3)s et %(k4)s.
  -> Conseil:
     Revoyez la définition de ces modèles.
"""
    ),
    84: _(
        """
  On ne peut pas créer un modèle X-FEM avec contact dans le cas où le mot clé MODELE_THER est présent.
  -> Conseil:
     Vous devez renseigner CONTACT='NON' pour ne pas activer le contact.
"""
    ),
    85: _(
        """
  La maille %(k1)s est affectée par un élément fini thermique dans le modèle thermique enrichi %(k2)s,
  or cette maille n'est affectée par aucun élément fini dans le modèle mécanique sain %(k3)s.
  -> Conseil:
     Revoyez la définition de ces deux modèles.
"""
    ),
    86: _(
        """
  -> Il n'est pas possible de réaliser la propagation d'une fissure en présence de mailles
     quadratiques dans le cadre d'un modèle X-FEM.

  -> Risque & Conseil:
     Veuillez utiliser un maillage linéaire.
"""
    ),
    86: _(
        """
Le CHAM_GD utilisé par l'opération ASSE_DEPL de l'opérateur CREA_CHAMP doit être un champ nodal."
"""
    ),
    87: _(
        """
Le CHAM_GD utilisé par l'opération ASSE_DEPL de l'opérateur CREA_CHAMP doit être de type DEPL_R."
"""
    ),
    88: _(
        """
L'opération ASSE_DEPL de l'opérateur CREA_CHAMP ne doit pas être utilisée si le modèle le comporte pas d'éléments X-FEM."
"""
    ),
    89: _(
        """
L'opération ASSE_DEPL de l'opérateur CREA_CHAMP ne prend pas en charge les éléments HM."
"""
    ),
    95: _(
        """
  -> Le mot-clé CRITERE de PROPA_FISS vaut 'ANGLE_IMPO_GAMMA' ou 'ANGLE_IMPO_BETA_GAMMA' et le tableau
     des facteurs d'intensité de contraintes de la fissure %(k1)s ne contient
     pas de colonne 'GAMMA'.
  -> Risque & Conseil:
     Si vous souhaitez imposer les valeurs de l'angle de déversement aux points
     du fonds de fissure, veuillez indiquer CRITERE='ANGLE_IMPO_GAMMA' ou 'ANGLE_IMPO_BETA_GAMMA' et ajouter
     une colonne 'GAMMA' au tableau des facteurs d'intensité de contraintes manuellement.
"""
    ),
    96: _(
        """
  -> Vous avez renseigné le mot-clé EVOL de AFFE_MATERIAU / AFFE_VARC avec un résultat thermique X-FEM.
     Dans ce cas, vous ne pouvez définir qu'une seule occurrence du mot-clé facteur AFFE_VARC.
"""
    ),
    97: _(
        """
  -> Vous avez renseigné le mot-clé EVOL de AFFE_MATERIAU / AFFE_VARC avec un résultat thermique X-FEM.
     Dans ce cas, le mot-clé MODELE de AFFE_MATERIAU devient obligatoire.
"""
    ),
    98: _(
        """
  -> Les chargements de type FORCE_FACE ne sont pas gérés pour le calcul de l'option CALC_G_XFEM sur les éléments de bord X-FEM.
     Seuls les chargements de type PRES_REP peuvent être pris en compte.
"""
    ),
    99: _(
        """
  -> L'opérateur CALC_G_XFEM ne sait pas traiter le cas d'une pression fonction d'un ou de plusieurs paramètres appliquée
     sur les lèvres d'une fissure X-FEM. Seul le cas d'une pression constante en espace et en temps est prévu.
     Veuillez utiliser l'opérateur POST_K1_K2_K3.
"""
    ),
    100: _(
        """
  -> Le champ matériau est absent de la liste des arguments de l'opérateur
     POST_ERREUR, pour une fissure XFEM, avec fond de fissure. L'assemblage
     du champ de déplacement X-FEM est impossible.
  -> Conseil:
     Rajouter le CHAM_MATER en entrée de la commande POST_ERREUR.
"""
    ),
}
