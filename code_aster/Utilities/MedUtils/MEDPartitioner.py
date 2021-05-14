# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
# ----------------------------------- ---------------------------------

# person_in_charge: nicolas.pignet at edf.fr

from mpi4py import MPI
import medcoupling as mc
import string
import random

from .myMEDSplitter_mpi import BuildPartNameFromOrig, MakeThePartition, GetGraphPartitioner, setVerbose
from .deal_with_pt1_post_process import add_pt1

from ...Messages import UTMESS
from ..logger import logger


class MEDPartitioner:

    """This class allows to read, partitioning and write a MED mesh in parallel

    """

    def __init__(self, filename, meshname=None):
        """Initialization
        :param filename: name of MED file
        :type filename: string
        :param meshname: name of mesh
        :type meshname: string
        """
        self._filename = filename
        self._meshname = meshname
        self._meshPartitioned = None
        self._writedFilename = None

        med_vers_min = 300 # minimal med version 3.0.0
        med_vers_num = int(mc.MEDFileVersionOfFileStr(self._filename).replace(".", ""))
        if med_vers_num < med_vers_min:
            UTMESS("F", "MED_20", valk=("3.0.0", mc.MEDFileVersionOfFileStr(self._filename)))

    def filename(self):
        """ Return the name of MED file"""
        return self._filename

    def writedFilename(self):
        """ Return the name of the writed MED file"""
        return self._writedFilename

    def meshname(self):
        """ Return the name of the mesh"""
        return self.meshname

    def getPartition(self):
        """ Return the partition of the mesh"""
        return self._meshPartitioned

    def partitionMesh(self, verbose=0):
        """ Partition the mesh

        Arguments:
            verbose (int) : 0 - warnings
                            1 - informations about main steps
                            2 - informations about all steps
        """
        level = logger.getEffectiveLevel()
        setVerbose(verbose, True)

        self._meshPartitioned = MakeThePartition(self._filename, self._meshname, GetGraphPartitioner(None))

        logger.setLevel(level)

    def writeMesh(self, path=None):
        """ Write the partitioning mesh file in MED format

         Arguments:
            path (str): path where to write the file
        """

        if self._meshPartitioned is not None:
            if self.filename()[0:4] == "fort":
                new_name = ''.join(random.choice(string.ascii_lowercase) for i in range(8))
            else:
                new_name = self.filename()

            full_path = new_name
            if (path is not None):
                if path.endswith("/"):
                    full_path = path + new_name
                else:
                    full_path = path + "/" + new_name

            self._writedFilename = BuildPartNameFromOrig(full_path, MPI.COMM_WORLD.rank)
            self._meshPartitioned.write33( self.writedFilename(), 2)

    def addPO1(self):
        """ Add PO1 in the partitionning meshes
            Be carefull all the meshes have to be in the same folder and the initial mesh
            is loaded so it can be memory expensive

            This is a temporary function until the partitionner do it itself
        """

        if self._writedFilename is not None:
            if MPI.COMM_WORLD.rank == 0:
                name_files = self._writedFilename.replace("0.med", "*.med")
                args_dict = {"verbosity": 2, "nb_dest_par": MPI.COMM_WORLD.size, \
                    "origin_seq": self._filename, "dest_par": name_files, "force": 1}
                add_pt1(args_dict)

            MPI.COMM_WORLD.barrier()
