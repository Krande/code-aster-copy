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

"""
Outils d'interpolation en correction de fonctionnalités incomplètes de medcoupling
"""
import math
import numpy as np
import medcoupling as mc

from ...Utilities import logger
from ...Messages import UTMESS


class CELL_TO_POINT:
    """
    GetCellsContainingPoints peut échouer en présence de faces gauches
    On propose ici une stratégie robuste pour passer de valeurs P0 (sur cellules) à un nuage de points
    En chaque point, on retient la valeur du champ dans la cellule à laquelle appartient le point
    Si le point appartient à plusieurs cellules (précision), on fait la moyenne arithmétique des cellules conscernées
    A noter que la précision s'entend comme une valeur relative : multipliée par la taille de la boite englobant la cellule,
    elle fournit le rayon (absolu) de la boule pour la détection de cellules ou relatif
    Si le point n'appartient à aucune cellule, on augmente la précision jusqu'à trouver au moins une cellule
    """

    def __init__(self, mesh_3D, mesh_2D, prec_rel, d_max_3d, ampl=1.5):
        """
        mesh_3D    : maillage de cellules 3D (celui sur lequel seront aussi définis les champs à projeter)
        mesh_2D  : maillage 2D sur lequel les champs seront projetés
        prec_rel : précision relative minimale par rapport à la taille de maille du maillage 2D
        ampl    : facteur d'amplification de la précision lorsqu'on ne trouve pas de cellules correspondant à un point
        """

        self.mesh = mesh_3D
        self.coor = mesh_2D.computeIsoBarycenterOfNodesPerCell()
        self.nbr = len(self.coor)

        self.prec = None
        self.cells = None
        self.pos = None
        self.pt_idx = None

        self.ComputeRelativePrecision(mesh_2D, prec_rel, d_max_3d)
        self.Build(ampl)

    def ComputeRelativePrecision(self, mesh_2D, prec_rel, d_max_3D):
        """
        Estime la précision relative de l'identification des cellules correspondant aux points
        en s'appuyant sur la taille des mailles 2D et la précision relative donnée en entrée
        """
        prec_abs = prec_rel * mesh_2D.computeDiameterField().getArray().getMinValueInArray()
        self.prec = prec_abs / d_max_3D

    def Build(self, ampl):
        """
        Construction du projecteur, à savoir le lien point -> cellules correspondantes
        Identique à getCellsContainingPoints mais en plus robuste vis-à-vis de la précision
        car on augmente progressivement (facteur ampl) la précision jusqu'à ce qu'à chaque point
        corresponde au moins une cellule
        """

        # ------------------------------------------------------------------------------------------
        # 1 - Le maillage support du champ est découpé en tétraèdres
        # ------------------------------------------------------------------------------------------

        meshT4 = self.mesh.deepCopy()

        if not meshT4.getAllGeoTypes() == [mc.NORM_TETRA4]:
            (meshT4, transfer, ibid) = meshT4.tetrahedrize(
                mc.PLANAR_FACE_6
            )  # ne respecte pas la conformité, ce qui explique tout le développement
        else:
            transfer = mc.DataArrayInt(np.arange(meshT4.getNumberOfCells()).tolist())
        meshT4 = meshT4.buildUnstructured()  # Bug medcoupling si on ne fait rien

        # ------------------------------------------------------------------------------------------
        # 2- On repère les cellules T4 contenant les points avec un accroissement progressif de la précision jusqu'à ce que tous les points aient une correspondance
        # ------------------------------------------------------------------------------------------

        prec = self.prec
        pt_idx_inv = mc.DataArrayInt([])
        cells = mc.DataArrayInt([])
        pos = mc.DataArrayInt([0])
        nook_idx = mc.DataArrayInt(np.arange(self.nbr).tolist())

        # boucle sur la précision de recherche
        nb_iter = 0
        while 1:

            logger.info("Projector build ... ")
            logger.info("Precision courante = " + str(prec))

            coor_nook = self.coor[nook_idx].toNumPyArray().ravel().tolist()
            (cells_prec, pos_prec) = meshT4.getCellsContainingPoints(coor_nook, prec)
            nbCells = pos_prec[1:] - pos_prec[:-1]

            logger.info(
                "Nombre maximal de nouvelles cellules (points de Gauss transformés en TETRA) 3D trouvées = "
                + str(nbCells.getMaxValueInArray())
            )
            logger.info(
                "Nombre minimal de nouvelles cellules (points de Gauss transformés en TETRA) 3D trouvées = "
                + str(nbCells.getMinValueInArray())
            )
            logger.info("")

            # Amplification de la précision
            prec = ampl * prec
            nb_iter += 1

            # nb_iter max = 10
            if nb_iter > 10:
                UTMESS("F", "RUPTURE4_23")

            # S'il n'y a pas de nouveaux points OK, on continue avec une précision dégradée
            if nbCells.getMaxValueInArray() == 0:
                if nb_iter > 10:
                    UTMESS("F", "RUPTURE4_23")
                continue

            # Index des nouveaux points pour lesquels des cellules sont identifiées
            found = nbCells.findIdsNotEqual(0)
            ok_idx = nook_idx[found]
            logger.info("On a trouvé " + str(len(ok_idx)) + str(" nouveaux points."))
            logger.info("")
            # Concaténation des points et des cellules
            pt_idx_inv.aggregate(ok_idx)
            pos = pos[:-1]
            pos.aggregate(pos_prec[found] + len(cells))
            cells.aggregate(cells_prec)
            pos.aggregate(mc.DataArrayInt([len(cells)]))

            # Tous les points ont-ils trouvé leurs cellules d'appartenance
            if nbCells.getMinValueInArray() > 0:
                break

            # Points résiduels qui n'ont toujours pas trouvé leurs cellules d'appartenance
            notFound = nbCells.findIdsEqual(0)
            nook_idx = nook_idx[notFound]
            logger.info("Il reste " + str(len(nook_idx)) + str(" points à trouver."))
            logger.info("")
            if nb_iter > 10:
                UTMESS("F", "RUPTURE4_23")

        if nb_iter > 1:
            UTMESS("A", "RUPTURE4_22")

        # Interversion de l'indexation des points : numéro du point utilisateur -> numéro du point du projecteur
        self.pt_idx = mc.DataArrayInt([-1] * self.nbr)
        self.pt_idx[pt_idx_inv] = mc.DataArrayInt(np.arange(self.nbr).tolist())

        # Cellules du maillage d'origine pour chacun des points
        self.cells = transfer[cells]
        self.pos = pos

        # Verifications
        assert (
            self.pos[1:] - self.pos[:-1]
        ).getMinValueInArray() > 0  # on a bien trouvé au moins une cellule pour chaque point
        assert len(self.pos) == self.nbr + 1

    def GetCellsContainingPoints(self):
        """
        On renvoie les éléments caractéristiques du projecteur (cells, pos, pt_idx)
        cells et pos ont la même signification que la méthode éponyme de medcoupling
        En revanche, on a une indirection supplémentaire : cells et pos sont repérés par rapport
        à la numérotation des points interne au projecteur
        pt_idx[numero_point_utilisateur] = numero_point_projecteur
        """
        return (self.cells.deepCopy(), self.pos.deepCopy(), self.pt_idx.deepCopy())

    def Eval(self, fieldValues):
        """
        Evaluation du champ scalaire P0 sur les points d'échantillonnage (un DataArrayDouble)
        """

        # Coefficient de pondération propre à chaque cellule
        nbCells = self.pos[1:] - self.pos[:-1]
        weights = mc.DataArrayDouble((1.0 / nbCells.toNumPyArray()).tolist())

        # Calcul des moyennes arithmétiques dans la numérotation des points propre au projecteur
        offset = 0
        values = mc.DataArrayDouble(self.nbr, fieldValues.getNumberOfComponents())
        values.fillWithValue(0.0)

        while nbCells.getMaxValueInArray() > 0:
            # Points pour lesquels il reste des cellules à calculer
            active_idx = nbCells.findIdsGreaterOrEqualTo(1)

            # Contribution des cellules courante à la moyenne
            cells = self.cells[self.pos[active_idx] + offset]
            values[active_idx] = values[active_idx] + weights[active_idx] * fieldValues[cells]

            # Préparation des pointeurs pour les termes suivants
            nbCells = nbCells - 1
            offset += 1

        # indexation des valeurs dans la numérotation de l'utilisateur
        values = values[self.pt_idx]

        return values
