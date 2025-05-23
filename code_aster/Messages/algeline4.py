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
 Seules les méthodes de résolution LDLT, MUMPS et MULT_FRONT sont autorisées.
"""
    ),
    2: _(
        """
 Solveur modal TRI_DIAG + option de filtrage des modes rigides:
 Le shift utilisé pour poursuivre le calcul modal est quasi-nul: shift=%(r1)g.
 Cela va perturber le calcul modal. Ce problème aurait dû être détecté automatiquement
 et le shift aurait dû être décalé.

 Conseils :
   * Vérifier votre mise  en données concernant les paramètres de décalage (PREC_SHIFT, NMAX_ITER_SHIFT) ou concernant la
      détection de singularité (SOLVEUR/NPREC).
   * Vous pouvez aussi décaler manuellement le shift via les valeurs 'CENTRE' ou 'BANDE' du mot-clé OPTION.
   * Vous pouvez aussi relancer le même calcul en changeant de solveur modal (par exemple avec METHODE='SORENSEN').
"""
    ),
    3: _(
        """
 Solveur GCPC :
 La résolution du système linéaire a échoué car le critère de convergence n'a pu être satisfait avec le nombre d'itérations autorisées.
   norme du résidu relatif après %(i1)d itérations :  %(r1)f
   critère de convergence à satisfaire             :  %(r2)f

 Conseils :
  * Augmentez le nombre d'itérations autorisées (SOLVEUR/NMAX_ITER).
  * Vous pouvez aussi augmenter le niveau de remplissage pour la factorisation incomplète (SOLVEUR/NIVE_REMPLISSAGE).
  * Si vous utilisez une commande non-linéaire (STAT_NON_LINE par exemple), diminuez la précision demandée pour la convergence (SOLVEUR/RESI_RELA).
    Prenez garde cependant car cela peut empêcher la convergence de l'algorithme non-linéaire.
"""
    ),
    4: _(
        """
  Manque de mémoire :
     Mémoire disponible = %(i1)d
     Mémoire nécessaire = %(i2)d
"""
    ),
    5: _(
        """
Erreur utilisateur dans la commande CREA_MAILLAGE :
  Pour créer le nouveau maillage, il faut créer de nouveaux noeuds.
  Ici, on cherche à créer le noeud (%(k1)s) mais il existe déjà.

  Le préfixe peut être modifié par l'utilisateur (mots clés XXXX / PREF_NOEUD).

Risques & conseils :
  Quand on utilise deux fois de suite la commande CREA_MAILLAGE, il est en général nécessaire
  d'utiliser au moins une fois l'un des mots clés PREF_NOEUD
"""
    ),
    6: _(
        """
 Solveur GCPC :
 La résolution du système linéaire a échoué car le critère de convergence n'a pu être satisfait avec le nombre d'itérations autorisées.
   norme du résidu relatif après %(i1)d itérations :  %(r1)f
   critère de convergence à satisfaire             :  %(r2)f

 Conseils :
  * Augmentez le nombre d'itérations autorisées (SOLVEUR/NMAX_ITER).
  * Vous pouvez aussi réactualiser plus souvent le préconditionneur en diminuant la valeur du mot-clé SOLVEUR/REAC_PRECOND.
  * Si vous utilisez une commande non-linéaire (STAT_NON_LINE par exemple), diminuez la précision demandée pour la convergence (SOLVEUR/RESI_RELA).
    Prenez garde cependant car cela peut empêcher la convergence de l'algorithme non-linéaire.
"""
    ),
    9: _(
        """
 Erreur données GROUP_MA déjà existant :  %(k1)s
"""
    ),
    10: _(
        """
 Votre problème modal n'est pas un problème généralisé à matrices réelles symétriques :
 il comporte des matrices non symétriques et/ou complexes, ou bien il s'agit d'un problème quadratique.
 Son spectre n'est donc pas uniquement restreint à l'axe réel, il représente une zone du plan complexe.

 Conseil :
 Il faut donc relancer votre calcul avec le mot-clé TYPE_RESU='MODE_COMPLEXE' et les opérandes associées.
"""
    ),
    11: _(
        """
 erreur données GROUP_NO déjà existant :  %(k1)s
"""
    ),
    13: _(
        """
 L'algorithme APM a atteint le nombre maximal de discrétisations du contour,
 c'est à dire %(i1)d, sans convergence du procédé.

 Conseils :
 Vous pouvez:
  - Changer les dimensions du contour de manière à réduire son périmètre,
  - Changer sa localisation. Il passe peut-être très près de valeurs propres
    ce qui peut induire des perturbations numériques.
"""
    ),
    14: _(
        """
 L'algorithme APM a atteint son nombre maximal d'itérations, c'est à dire %(i1)d,
 sans convergence du procédé.

 Conseils :
  Vous pouvez:
  - Augmenter ce nombre maximal d'itérations via le paramètre NMAX_ITER_CONTOUR,
  - Augmenter la discrétisation initiale du contour via NBPOINT_CONTOUR,
  - Changer les dimensions du contour de manière à réduire son périmètre,
  - Changer sa localisation. Il passe peut-être très près de valeurs propres
    ce qui peut induire des perturbations numériques.
"""
    ),
    15: _(
        """
 L'algorithme APM avec le calcul du polynôme caractéristique via une factorisation
 LDLT a un problème numérique: le point de vérification (%(r1)f +i*%(r2)f)
 est très proche d'une valeur propre ou le solveur linéaire a eu un problème.

 Conseils :
 Vous pouvez:
 - Augmenter les dimensions du contour pour englober cette valeur propre,
 - Changer la discrétisation du contour (plus risqué).
 - Changer le paramétrage du solveur linéaire, ou le solveur linéaire lui-même (expert).
"""
    ),
    19: _(
        """
 Matrice de masse non définie.

 Conseil : essayer un autre algorithme de résolution.
"""
    ),
    20: _(
        """
 Pour l'instant, on est obligé de choisir pour un résultat de type 'DYNAMIQUE' ou
 'FLAMBEMENT', la méthode de comptage 'STURM', et pour 'MODE_COMPLEXE', la méthode
 'APM'.
 Si vos choix ne respectent pas cette règle, on fait le changement pour vous, en
 se référant au type de problème que vous avez choisi.
"""
    ),
    21: _(
        """
 Manque de place mémoire longueur de bloc insuffisante:  %(i1)d
 le super noeud  %(i2)d
  nécessite un bloc de %(i3)d
"""
    ),
    22: _(
        """
 L'algorithme APM a convergé sur un nombre de fréquences aberrant !

 Conseils:
 Vous pouvez:
  - Augmenter la discrétisation initiale du contour via NBPOINT_CONTOUR,
  - Changer les dimensions du contour de manière à réduire son périmètre,
  - Changer sa localisation. Il passe peut-être très près de valeurs propres
    ce qui peut induire des perturbations numériques.
"""
    ),
    25: _(
        """
 combinaison non prévue   type résultat :  %(k1)s    type matrice  :  %(k2)s
    type constante:  %(k3)s
"""
    ),
    27: _(
        """
 combinaison non prévue
 type résultat :  %(k1)s
 type matrice  :  %(k2)s
"""
    ),
    31: _(
        """
 combinaison non prévue
 type résultat :  %(k1)s
"""
    ),
    33: _(
        """
 la normalisation doit se faire en place
 il est impossible d'avoir comme concept produit  %(k1)s et %(k2)s comme concept d'entrée.
"""
    ),
    36: _(
        """
 l'option de normalisation  %(k1)s  n'est pas implantée. %(i1)d
"""
    ),
    37: _(
        """
 problème(s) rencontré(s) lors de la factorisation de la matrice : %(k1)s
"""
    ),
    38: _(
        """
 appel erroné :
 code retour : %(i1)d
 Problème CHAM_NO %(k1)s
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    39: _(
        """
 Au moins une des matrices est non symétrique.
 Dans ce cas, l'option 'BANDE' n'est pas utilisable.
"""
    ),
    40: _(
        """
 Au moins une des matrices est non symétrique.
 Dans ce cas, la détection des modes de corps rigide (OPTION='MODE_RIGIDE')
 n'est pas utilisable.
"""
    ),
    41: _(
        """
 Au moins une des matrices est non symétrique.
 Dans ce cas, le calcul de flambement ne peut pas être mené.
"""
    ),
    42: _(
        """
 pas de produit car les valeurs de la MATRICE sont  %(k1)s
 et celles du CHAM_NO sont  %(k2)s
"""
    ),
    44: _(
        """
 Le solveur itératif GCPC est interdit avec un maillage parallèle.

 Conseils :
  - Utilisez un autre solveur.
"""
    ),
    55: _(
        """
 pas d'extraction pour  %(k1)s
 pour le numéro d'ordre  %(i1)d
"""
    ),
    56: _(
        """
 pas de mode extrait pour  %(k1)s
"""
    ),
    57: _(
        """
 NUME_MODE identique pour le %(i1)d
 mode d'ordre  %(i2)d
"""
    ),
    58: _(
        """
  Problème dans le préconditionnement de la matrice MATAS par LDLT incomplet
  pivot nul à la ligne %(i1)d correspondant au noeud %(k1)s et à la
  composante %(k2)s
"""
    ),
    60: _(
        """
  incohérence sans multiplicateurs de Lagrange %(i1)d reconstitués %(i2)d
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    61: _(
        """
 pas de mode statique pour le noeud :  %(k1)s  et sa composante :  %(k2)s
"""
    ),
    62: _(
        """
 pour les modes statiques, on attend un :  %(k1)s
 noeud :  %(k2)s
 composante   :  %(k3)s
"""
    ),
    63: _(
        """
 Solveur GCPC :
 Le préconditionneur LDLT a été mis à jour.
 """
    ),
    64: _(
        """
 détection d'un terme nul sur la sur diagonale
 valeur de BETA   %(r1)f
 valeur de ALPHA  %(r2)f
"""
    ),
    65: _(
        """
 La  %(i1)d -ème valeur propre du système réduit est complexe.
 Partie imaginaire =  %(r1)f
 et partie imaginaire / partie réelle =  %(r2)f
 """
    ),
    66: _(
        """
 la valeur propre est :   %(r1)f
"""
    ),
    74: _(
        """
 Calcul d'erreur modale :
 une valeur propre réelle est détectée à partir du couple (fréquence, amortissement réduit).
 On ne peut plus la reconstruire.
 Par convention l'erreur modale est fixée à : %(r1)f.
"""
    ),
    76: _(
        """
 la réorthogonalisation diverge après  %(i1)d  itération(s).
"""
    ),
    77: _(
        """
 l'option de normalisation  %(k1)s  n'est pas implantée.
"""
    ),
    80: _(
        """
 type de valeurs inconnu   %(k1)s
"""
    ),
    82: _(
        """
 incohérence de certains paramètres modaux propres à ARPACK
 numéro d'erreur  %(i1)d
"""
    ),
    85: _(
        """
 appel erroné mode numéro %(i1)d position modale %(i2)d
 code retour : %(i3)d
 Problème CHAM_NO %(k1)s
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    86: _(
        """
 la réorthogonalisation diverge après  %(i1)d  itération(s) %(i2)d
       vecteur traité :  %(i3)d
       vecteur testé  :  %(i4)d
 arrêt de la réorthogonalisation %(k1)s
"""
    ),
    87: _(
        """
 pour le problème réduit
 valeur(s) propre(s) réelle(s)                  :  %(i1)d
 valeur(s) propre(s) complexe(s) avec conjuguée :  %(i2)d
 valeur(s) propre(s) complexe(s) sans conjuguée :  %(i3)d

 Indications et Conseils
 -----------------------
 * Une valeur propre est considérée comme réelle dès que le module de sa partie imaginaire (en pulsation) est < 1.d-7.
 * Deux valeurs propres sont considérées comme conjuguées dès que le module de leur écart (en pulsation)
 est inférieur à (2.pi.SEUIL_FREQ)**2.
 Vous pouvez donc augmenter un peu la valeur par défaut de SEUIL_FREQ si vous pensez que la calcul a filtré, à tort,
 des couples de valeurs propres normalement conjuguées.
"""
    ),
    88: _(
        """
 votre problème est fortement amorti.
 valeur(s) propre(s) réelle(s)                  :  %(i1)d
 valeur(s) propre(s) complexe(s) avec conjuguée :  %(i2)d
 valeur(s) propre(s) complexe(s) sans conjuguée :  %(i3)d

 Indications et Conseils
 -----------------------
 * Une valeur propre est considérée comme réelle dès que le module de sa partie imaginaire (en pulsation) est < 1.d-7.
 * Deux valeurs propres sont considérées comme conjuguées dès que le module de leur écart (en pulsation)
 est inférieur à (2.pi.SEUIL_FREQ)**2.
 Vous pouvez donc augmenter un peu la valeur par défaut de SEUIL_FREQ si vous pensez que la calcul a filtré, à tort,
 des couples de valeurs propres normalement conjuguées.
"""
    ),
    93: _(
        """
 Problème généralisé complexe.
"""
    ),
    94: _(
        """
 Problème quadratique complexe.
"""
    ),
    95: _(
        """
 Problème quadratique.
"""
    ),
    96: _(
        """
 Amortissement (réduit) de décalage supérieur en valeur absolue à %(r1)f.
 On le ramène à la valeur : %(r2)f.
"""
    ),
}
