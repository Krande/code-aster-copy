# coding: utf-8

# Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.

# person_in_charge: nicolas.sellenet@edf.fr

from ..Cata.Language.SyntaxObjects import FactorKeyword
from ..Objects import (MechanicalLoadReal, ParallelMechanicalLoadReal,
                       ConnectionMesh, Model)
from ..Supervis import ExecuteCommand
from ..Utilities import deprecate, force_list


class MechanicalLoadDefinition(ExecuteCommand):

    """Command that defines :class:`~code_aster.Objects.MechanicalLoadReal`.
    """
    command_name = "AFFE_CHAR_MECA"
    neumannLoads = ['PESANTEUR', 'ROTATION', 'FORCE_FACE', 'FORCE_ARETE', 'FORCE_CONTOUR', 'FORCE_INTERNE', 'PRE_SIGM', 'PRES_REP', 'EFFE_FOND', 'PRE_EPSI', 'FORCE_POUTRE', \
                    'FORCE_TUYAU', 'FORCE_COQUE', 'FORCE_ELEC', 'INTE_ELEC', 'VITE_FACE', 'ONDE_FLUI', 'FLUX_THM_REP', 'FORCE_SOL', 'FORCE_NODALE', 'EVOL_CHAR', 'VECT_ASSE']
    dirichletLoads = ['DDL_IMPO', 'DDL_POUTRE', 'FACE_IMPO', 'CHAMNO_IMPO', 'ARETE_IMPO', 'LIAISON_DDL', 'LIAISON_OBLIQUE', 'LIAISON_GROUP', 'LIAISON_MAIL', 'LIAISON_PROJ', \
                      'LIAISON_CYCL', 'LIAISON_SOLIDE', 'LIAISON_ELEM', 'LIAISON_UNIF', 'LIAISON_CHAMNO', 'LIAISON_RBE3', 'LIAISON_INTERF', 'LIAISON_COQUE', \
                      'RELA_CINE_BP', 'IMPE_FACE']

    @staticmethod
    def compat_syntax(keywords):
        """Replace LIAISON='ENCASTRE' by BLOCAGE.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords, changed
                in place.
        """
        if not keywords.get("DDL_IMPO"):
            return
        # replace DDL_IMPO/LIAISON=ENCASTRE by DDL_IMPO/BLOCAGE
        keywords["DDL_IMPO"] = force_list(keywords.get("DDL_IMPO", []))
        for fact in keywords["DDL_IMPO"]:
            block = fact.pop("LIAISON", None)
            if block == 'ENCASTRE':
                deprecate("DLL_IMPO/LIAISON='ENCASTRE'", case=3, level=5,
                          help="Use BLOCAGE = ('DEPLACEMENT', 'ROTATION')")
                fact["BLOCAGE"] = ('DEPLACEMENT', 'ROTATION')

    def _addGroup(self, mcf, groups, keys):
        for name in keys:
            mc = mcf.get(name)
            if mc:
                groups.update(force_list(mc))

    def _excludeGroup(self, mcf, keys):
        for name in keys:
            if mcf.get(name):
                raise RuntimeError("Keyword %s not accepted in parallel AFFE_CHAR_MECA"%name)

    def _getGroups(self, keywords):
        """for parallel load, return all node and cells groups present in AFFE_CHAR_MECA, in order to define the connection mesh
        """
        load_types = [key for key in list(self._cata.keywords.keys()) if isinstance(self._cata.keywords[key], FactorKeyword)]
        nodeGroups = set()
        cellGroups = set()
        for key in list(keywords.keys()):
            if key in ("LIAISON_DDL", "DDL_IMPO", "LIAISON_OBLIQUE", "LIAISON_UNIF", \
                            "LIAISON_SOLIDE", "DDL_POUTRE", "FACE_IMPO", "ARETE_IMPO"):
                for mcf in keywords[key]:
                    self._addGroup(mcf, nodeGroups, ("GROUP_NO", "SANS_GROUP_NO"))
                    self._addGroup(mcf, cellGroups, ("GROUP_MA", "SANS_GROUP_MA"))
                    self._excludeGroup(mcf, ("NOEUD", "SANS_NOEUD", "MAILLE", "SANS_MAILLE"))
            elif key in ("LIAISON_GROUP", "LIAISON_COQUE", "LIAISON_ELEM"):
                for mcf in keywords[key]:
                    self._addGroup(mcf, nodeGroups, ("GROUP_NO_1", "GROUP_NO_2", "SANS_GROUP_NO"))
                    self._addGroup(mcf, cellGroups, ("GROUP_MA_1", "GROUP_MA_2"))
                    self._excludeGroup(mcf, ("NOEUD_1", "NOEUD_2", "SANS_NOEUD", "MAILLE_1", "MAILLE_2"))
            elif key in ("LIAISON_RBE3", "LIAISON_MAIL"):
                for mcf in keywords[key]:
                    self._addGroup(mcf, nodeGroups, ("GROUP_NO_MAIT", "GROUP_NO_ESCL"))
                    self._addGroup(mcf, cellGroups, ("GROUP_MA_MAIT", "GROUP_MA_ESCL"))
                    self._excludeGroup(mcf, ("NOEUD_ESCL", "NOEUD_MAIT", "MAILLE_MAIT", "MAILLE_ESCL"))
            elif key in load_types:
                raise NotImplementedError("Type of load {0!r} not yet "
                                      "implemented in parallel".format(key))
        # must be sorted to be identical on all procs
        return sorted(list(nodeGroups)), sorted(list(cellGroups))

    def _hasDirichletLoadings(self, keywords):
        """return True if instance has Dirichlet loadings"""
        for key in self.dirichletLoads:
            if keywords.get(key):
                return True
        return False

    def _hasNeumannLoadings(self, keywords):
        """return True if instance has Neumann loadings"""
        for key in self.neumannLoads:
            if keywords.get(key):
                return True
        return False

    def _hasOnlyDDL_IMPO(self, keywords):
        """ No need to create ConnectionMesh for DDL_IMPO"""
        for key in self.dirichletLoads:
            mc = keywords.get(key)
            if mc and key != "DDL_IMPO" :
                return False
        return False
        # Not Enable for the moment
        # return True

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        model = keywords["MODELE"]
        l_neum = self._hasNeumannLoadings(keywords)
        l_diri = self._hasDirichletLoadings(keywords)
        if not model.getMesh().isParallel():
            self._result = MechanicalLoadReal(model)
        else :
            if l_neum:
                if l_diri:
                    raise TypeError("Not allowed to mix up Dirichlet and Neumann loadings in the same parallel AFFE_CHAR_MECA")
                else:
                    self._result = MechanicalLoadReal(model)
            if self._hasOnlyDDL_IMPO(keywords):
                self._result = MechanicalLoadReal(model)


    def exec_(self, keywords):
        """Override default _exec in case of parallel load
        """
        if isinstance(self._result, MechanicalLoadReal):
            super(MechanicalLoadDefinition, self).exec_(keywords)
        else:
            model = keywords.pop("MODELE")
            nodeGroups, cellGroups = self._getGroups(keywords)
            connectionMesh = ConnectionMesh(model.getMesh(), nodeGroups, cellGroups)

            connectionModel = Model( connectionMesh )
            connectionModel.setFrom( model )

            keywords["MODELE"] = connectionModel
            partialMechanicalLoad = AFFE_CHAR_MECA(**keywords)
            keywords["MODELE"] = model
            self._result = ParallelMechanicalLoadReal(partialMechanicalLoad, model)

AFFE_CHAR_MECA = MechanicalLoadDefinition.run
