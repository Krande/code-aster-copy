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
    3: _(
        """
 La conversion d'un champ par éléments aux noeuds en un champ
 aux noeuds a été demandé.
 Or le champ par éléments aux noeuds contient des
 sous-points alors qu'un champ aux noeuds ne peut pas en contenir.
 Les mailles contenant des sous-points vont être exclues.

 Conseil :
   Si vous souhaitez convertir le champ sur les mailles ayant des
   sous-points, vous devez d'abord extraire le champ par éléments
   aux noeuds sur un sous-point avec la commande POST_CHAMP.
"""
    ),
    4: _(
        """
 Vous avez demandé le calcul d'un champ aux noeuds sur des éléments
 de structure. Mais les repères locaux de certaines mailles entourant
 des noeuds sur lesquels vous avez demandés le calcul ne sont pas
 compatibles (Au maximum, on a %(r1)g degrés d'écart entre les angles
 nautiques définissant ces repères).

 Risque & Conseil :
   Il se peut que vous obteniez des résultats incohérents.
   Il est donc recommandé de passer en repère global les champs
   utiles au calcul du champ aux noeuds.
"""
    ),
    7: _(
        """
 Problème norme de axe
"""
    ),
    8: _(
        """
Cette grandeur ne peut pas accepter plus de %(i2)d composantes (%(i1)d fournies).
"""
    ),
    9: _(
        """
 dimension  %(k1)s  inconnue.
"""
    ),
    10: _(
        """
 maillage obligatoire.
"""
    ),
    11: _(
        """
 on ne peut pas créer un champ de VARI_R avec le mot clé facteur AFFE
 (voir u2.01.09)
"""
    ),
    12: _(
        """
 mot clé AFFE/NOEUD interdit ici.
"""
    ),
    13: _(
        """
 mot clé AFFE/GROUP_NO interdit ici.
"""
    ),
    14: _(
        """
 type scalaire non traité :  %(k1)s
"""
    ),
    15: _(
        """
 incohérence entre nombre de composantes et nombre de valeurs
"""
    ),
    16: _(
        """
 il faut donner un champ de fonctions
"""
    ),
    17: _(
        """
 les paramètres doivent être réels
"""
    ),
    18: _(
        """
 maillages différents
"""
    ),
    20: _(
        """
Erreur utilisateur dans la commande CREA_CHAMP.
 le champ  %(k1)s n'est pas de type réel
"""
    ),
    21: _(
        """
 on ne traite que des "CHAM_NO" ou des "CHAM_ELEM".
"""
    ),
    22: _(
        """
 la programmation prévoit que les entiers sont codés sur plus de 32 bits
 ce qui n'est pas le cas sur votre machine
"""
    ),
    23: _(
        """
 on ne trouve aucun champ.
"""
    ),
    24: _(
        """
 Lors de la recopie du champ, certaines composantes ont du être
 mises à zéro.

 En effet, certaines parties attendues dans le champ n'étaient
 pas présentes dans le champ initial, elles ont donc été mises à zéros.

 Ce problème peut survenir lorsque la numérotation du champ initial
 n'est pas intégralement incluse dans celle du champ final
"""
    ),
    25: _(
        """
Erreur utilisateur dans la commande CREA_RESU / TYPE_MAXI= ...
  La commande exige que tous les champs pour lesquels on veut extraire
  le maximum aient la même numérotation pour leurs composantes.
  Ce n'est pas le cas pour cette structure de donnée RESULTAT.

  Cela est certainement du à la présence de mailles tardives à certains instants
  de calcul et pas à d'autres.
"""
    ),
    26: _(
        """
Erreur utilisateur dans la commande CREA_RESU / TYPE_MAXI= ...
  La commande exige que tous les CHAM_NO pour lesquels on veut extraire
  le maximum aient la même numérotation pour leurs composantes.
  Ce n'est pas le cas pour cette structure de donnée RESULTAT.
"""
    ),
    27: _(
        """
 il faut donner un maillage.
"""
    ),
    28: _(
        """
 Le champ : %(k1)s ne peut pas être assemblé en :  %(k2)s
"""
    ),
    29: _(
        """
 La structure de donnée résultat est associée à un maillage différent de la structure de donnée NUME_DDL.
"""
    ),
    30: _(
        """
 Le numéro de la composante %(i1)d est invalide. Il doit être
 compris entre 1 et %(i2)d.
"""
    ),
    31: _(
        """
 NOM_CMP2 et NOM_CMP de longueur différentes.
"""
    ),
    32: _(
        """
 Grandeur incorrecte pour le champ : %(k1)s
 grandeur proposée :  %(k2)s
 grandeur attendue :  %(k3)s
"""
    ),
    33: _(
        """
 le mot-clé 'COEF_C' n'est applicable que pour un champ de type complexe
"""
    ),
    34: _(
        """
 développement non réalisé pour les champs aux éléments. vraiment désolé !
"""
    ),
    35: _(
        """
Erreur utilisateur dans la commande CREA_CHAMP.
 le champ  %(k1)s n'est pas de type complexe
"""
    ),
    36: _(
        """
Erreur dans la commande CREA_CHAMP.
 Le champ de sortie ne contient pas de composante %(k1)s. Les valeurs
 correspondantes ne seront pas prises en compte.
"""
    ),
    37: _(
        """
Erreur utilisateur dans la commande CREA_CHAMP.
Les options 'R2C' et 'C2R' ne sont autorisées que pour les CHAM_NO.
Le champ %(k1)s n'est pas un champ aux noeuds.
"""
    ),
    40: _(
        """
 structure de données inexistante : %(k1)s
"""
    ),
    41: _(
        """
 duplication maillage objet inconnu:  %(k1)s
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    42: _(
        """
CREA_CHAMP/ASSE :

Une ou plusieurs composantes demandées sont inexistantes dans le champ produit.
Si vous ne comprenez pas pourquoi, sachez que cela peut-être dû au champ d'entrée
de la première occurrence qui dans certains cas impose ses composantes.

Par exemple, si vous souhaitez créer un champ de type ELGA_EPSI_R avec
les composantes EPXX et EPYY à partir notamment d'un champ EPVC_ELGA
(extrait d'un résultat), lui aussi de type ELGA_EPSI_R, placé en première
occurrence, aucune des composantes souhaitées ne sera prise en compte
car elles ne font pas partie des composantes du champ créé.

Conseils :

1/ Si plusieurs occurrences sont présentes, vous pouvez en placer une autre
en première position s'il elle ne présente pas le même problème.

2/ Vous pouvez faire les choses en deux temps en passant par l'intermédiaire
d'un champ de type NEUT.

3/ Vous pouvez ajouter artificiellement un champ permettant uniquement d'avoir
les composantes souhaitées. Dans l'exemple ci-dessus, on pourrait extraire
d'un résultat le champ EPSI_ELGA et le donner en première occurrence avec
un COEF_R à 0.
"""
    ),
    43: _(
        """
Erreur utilisateur dans la commande CREA_CHAMP.
Il n'y a pas de composante %(k1)s dans le champ %(k2)s.
"""
    ),
    44: _(
        """
Les instants de la liste ne sont pas uniques selon le critère %(k1)s et précision %(r1)e.
"""
    ),
    46: _(
        """
 la fonction doit s appuyer sur un maillage pour lequel une abscisse curviligne a été définie.
"""
    ),
    47: _(
        """
  le mot clé : %(k1)s n est pas autorise.
"""
    ),
    49: _(
        """
La question " %(k1)s " est inconnue
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    50: _(
        """
 champ de nom symbolique %(k1)s inexistant: %(k2)s
"""
    ),
    51: _(
        """
 il n y a pas de NUME_DDL pour ce CHAM_NO
"""
    ),
    52: _(
        """
 type de charge inconnu
"""
    ),
    54: _(
        """
 trop d objets
"""
    ),
    55: _(
        """
 champ inexistant: %(k1)s
"""
    ),
    56: _(
        """
 Le partitionneur SCOTCH a fait remonter l'erreur %(i1)d. Veuillez contacter l'équipe de
 développement.
"""
    ),
    57: _(
        """
La question n'a pas de réponse sur une grandeur de type matrice
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    59: _(
        """
La question n'a pas de sens sur une grandeur de type matrice
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    60: _(
        """
La question n'a pas de sens sur une grandeur de type composée.
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    61: _(
        """
La question n'a pas de sens, cette grandeur n'existe pas : %(k1)s
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    65: _(
        """
 le type de concept : " %(k1)s " est inconnu
"""
    ),
    66: _(
        """
 le phénomène :  %(k1)s  est inconnu.
"""
    ),
    68: _(
        """
 type de résultat inconnu :  %(k1)s  pour l'objet :  %(k2)s
"""
    ),
    70: _(
        """
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    71: _(
        """
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    72: _(
        """
 on ne traite pas les noeuds tardifs
"""
    ),
    73: _(
        """
 grandeur inexistante
"""
    ),
    74: _(
        """
 composante de grandeur inexistante
"""
    ),
    75: _(
        """
 problème avec la réponse  %(k1)s
"""
    ),
    76: _(
        """
 les conditions aux limites autres que des ddls bloques ne sont pas admises
"""
    ),
    77: _(
        """
 unité logique  %(k1)s , problème lors du close
"""
    ),
    78: _(
        """
 Vous essayez d'allouer dynamiquement un vecteur de dimension %(i1)d <= 0.
 Ce n'est pas possible.
"""
    ),
    82: _(
        """
 méthode d'intégration inexistante.
"""
    ),
    84: _(
        """
 interpolation  %(k1)s  non implantée
"""
    ),
    85: _(
        """
 recherche " %(k1)s " inconnue
"""
    ),
    86: _(
        """
 l'intitule " %(k1)s " n'est pas correct.
"""
    ),
    87: _(
        """
 le noeud " %(k1)s " n'est pas un noeud de choc.
"""
    ),
    88: _(
        """
 nom de sous-structure et d'intitulé incompatible
"""
    ),
    89: _(
        """
 le noeud " %(k1)s " n'est pas un noeud de choc de l'intitule.
"""
    ),
    90: _(
        """
 le noeud " %(k1)s " n'est pas compatible avec le nom de la sous-structure.
"""
    ),
    91: _(
        """
 le paramètre " %(k1)s " n'est pas un paramètre de choc.
"""
    ),
    92: _(
        """
 le noeud " %(k1)s " n'existe pas.
"""
    ),
    93: _(
        """
 la composante " %(k1)s " du noeud " %(k2)s " n'existe pas.
"""
    ),
    94: _(
        """
 type de champ inconnu  %(k1)s
"""
    ),
    95: _(
        """
 "INTERP_NUME" et ("INST" ou "LIST_INST") non compatibles
"""
    ),
    96: _(
        """
 "INTERP_NUME" et ("FREQ" ou "LIST_FREQ") non compatibles
"""
    ),
    99: _(
        """
 objet %(k1)s inexistant
"""
    ),
}
