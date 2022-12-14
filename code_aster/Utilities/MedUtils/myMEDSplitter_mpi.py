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

# person_in_charge: nicolas.pignet@edf.fr

import abc
import argparse
import logging
import os
import sys
from datetime import datetime
from distutils.version import StrictVersion
from resource import RUSAGE_SELF, getrusage

import med
import medcoupling as mc
import numpy as np

STANDALONE = True
try:
    from .. import MPI
    from ..logger import logger

    STANDALONE = False
except:
    logger = logging.getLogger()
    from mpi4py import MPI


def comm_world():
    """Dynamically returns the current COMM_WORLD.

    Returns:
        *mpi4py.MPI.Comm*: MPI Communicator.
    """
    return MPI.COMM_WORLD if STANDALONE else MPI.ASTER_COMM_WORLD


class ColoredFormatter(logging.Formatter):
    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(30, 38)
    COLORS = {"WARNING": YELLOW, "INFO": WHITE, "DEBUG": BLUE, "CRITICAL": YELLOW, "ERROR": RED}

    def __init__(self, *args, **kwargs):
        logging.Formatter.__init__(self, *args, **kwargs)

    def format(self, record):
        RESET_SEQ = "\033[0m"
        COLOR_SEQ = "\033[1;%dm"
        record.levelname = (
            COLOR_SEQ % ColoredFormatter.COLORS[record.levelname] + record.levelname + RESET_SEQ
        )
        return logging.Formatter.format(self, record)


def setVerbose(verbose=1):
    """Set verbosity level"""
    if STANDALONE:
        formatter = ColoredFormatter(
            "%(levelname)s #{} : %(asctime)s : %(message)s ".format(comm_world().rank), style="%"
        )
        formatter.default_time_format = "%H:%M:%S"
        formatter.default_msec_format = "%s.%03d"
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    verbose_map = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}
    logger.setLevel(verbose_map[verbose])


def memory_peak(mess=None, debug=False):
    # memory peak in Mo
    mem_used = int(getrusage(RUSAGE_SELF).ru_maxrss / 1024)
    if mess is not None:
        pat = 'Memory peak of "{}" : {} Mo'.format(mess, mem_used)
    else:
        pat = "Memory peak : {} Mo".format(mem_used)

    if debug:
        logger.debug(pat)
    else:
        logger.info(pat)


class ChronoCtxMgGen:
    def __init__(self, what):
        self._what = what

    def __enter__(self):
        pat = 'Start de "{}"'.format(self._what)
        self(pat)
        self.start = datetime.now()

    def __exit__(self, exctype, exc, tb):
        self.stop = datetime.now()
        pat = 'Fin de "{}"'.format(self._what)
        self(pat)
        delta = self.stop - self.start
        self('Time spent for "{}" : {}'.format(self._what, delta))
        memory_peak(self._what)


class ChronoCtxMg(ChronoCtxMgGen):
    def __init__(self, what):
        ChronoCtxMgGen.__init__(self, what)

    def __call__(self, pat):
        global logger
        logger.info(pat)


class ChronoCtxMgDbg(ChronoCtxMgGen):
    def __init__(self, what):
        ChronoCtxMgGen.__init__(self, what)

    def __call__(self, pat):
        logger.debug(pat)


class MasterChronoCtxMg:
    def __init__(self, what):
        self._what = what

    def __enter__(self):
        if comm_world().rank == 0:
            pat = 'Start de "{}"'.format(self._what)
            logger.info(pat)
            self.start = datetime.now()
            memory_peak(pat)

    def __exit__(self, exctype, exc, tb):
        if comm_world().rank == 0:
            self.stop = datetime.now()
            pat = 'Fin de "{0}"'.format(self._what)
            delta = self.stop - self.start
            logger.info(pat)
            logger.info('Time spent for "{}" : {}'.format(self._what, delta))
            memory_peak(pat)


def BuildPartNameFromOrig(fn, rank):
    return "{}_new_{}.med".format(os.path.splitext(fn)[0], rank)


def RetrieveFamGrpsMapInternal(fn, mn):
    famsPy = {}
    grpsPy = {}
    fid = med.MEDfileOpen(fn, med.MED_ACC_RDONLY)
    med.MEDnFamily(fid, mn)
    s = set()
    for elt in range(med.MEDnFamily(fid, mn)):
        nbGrps = med.MEDnFamilyGroup(fid, mn, elt + 1)
        gro = med.MEDCHAR(med.MED_LNAME_SIZE * nbGrps + 1)
        famName, famId, grps = med.MEDfamilyInfo(fid, mn, elt + 1, gro)
        famName = "FAM_" + str(famId)
        famsPy[famName] = famId
        grps2 = [
            "".join(gro[i * med.MED_LNAME_SIZE : (i + 1) * med.MED_LNAME_SIZE]).rstrip()
            for i in range(nbGrps)
        ]
        grps2 = [elt.rstrip("\x00") for elt in grps2]
        grpsPy[famName] = grps2
        s.add(famName)
        pass
    med.MEDfileClose(fid)
    return famsPy, grpsPy


def RetrieveFamGrpsMap(medFileContxt):
    return RetrieveFamGrpsMapInternal(medFileContxt.getFileName(), medFileContxt.getMeshName())
    # return mc.GetFamiliesGroupsInfo(medFileContxt.getFileName(),medFileContxt.getMeshName())


class MedFileContext:
    """Ease to get context of meshName contained in fileName"""

    def __init__(self, fileName, meshName):
        self._fileName = fileName
        self._meshName = meshName
        self.getMeshName()

    def getFileName(self):
        """return the fileName"""
        return self._fileName

    def getMeshName(self):
        """return the meshName"""
        if self._meshName is not None:
            return self._meshName
        if comm_world().rank == 0:
            meshNames = mc.GetMeshNames(self._fileName)
            if self._meshName is None or self._meshName not in meshNames:
                self._meshName = meshNames[0]
        else:
            self._meshName = None
        self._meshName = comm_world().bcast(self._meshName)
        return self._meshName

    def getNumberOfCells(self):
        if hasattr(self, "_nbCells"):
            return self._nbCells
        # 0 -> cell et 0 -> maxlevel
        tps = mc.GetUMeshGlobalInfo(self._fileName, self._meshName)[0][0]
        self._nbCells = sum([nbCells for typ, nbCells in tps])
        return self._nbCells

    def getGlobalNumberOfNodes(self):
        return mc.GetUMeshGlobalInfo(self._fileName, self._meshName)[-1]


class PartialMedFileUMesh:
    """Hold a partial umesh corresponding to current rank"""

    def __init__(self, medFileDataContext):
        self._medFileDataContext = medFileDataContext
        self.getMEDFileUMesh()

    def cellsGen(self, lev):
        rank = comm_world().rank
        size = comm_world().size
        tps = mc.GetUMeshGlobalInfo(
            self._medFileDataContext.getFileName(), self._medFileDataContext.getMeshName()
        )[0][abs(lev)]
        params = []
        acc = 0
        for cellType, nbCellsType in tps:
            slc = mc.DataArray.GetSlice(slice(0, nbCellsType, 1), rank, size)
            params.append(mc.DataArrayInt.Range(acc + slc.start, acc + slc.stop, slc.step))
            acc += nbCellsType
            pass
        return mc.DataArrayInt.Aggregate(params)

    @property
    def cells(self):
        return self.cellsGen(0)

    def storeOrphanNodesIfAny(self, arrayOfGlobalIdsOfOrphanNodes):
        """
        Args
            arrayOfGlobalIdsOfOrphanNodes : tous les noeuds a charger
        """
        if arrayOfGlobalIdsOfOrphanNodes.empty():
            self._orphan_coords = mc.DataArrayDouble([])
            self._orphan_global_ids = mc.DataArrayInt([])
            self._orphan_fam_ids = mc.DataArrayInt([])
            self._orphan_num_ids = mc.DataArrayInt([])
        else:
            dt, it, _ = self._medFileUMesh.getTime()
            vmin = arrayOfGlobalIdsOfOrphanNodes.getMinValueInArray()
            vmax = arrayOfGlobalIdsOfOrphanNodes.getMaxAbsValueInArray()
            dataFromFile = mc.MEDFileUMesh.LoadPartCoords(
                self._medFileDataContext.getFileName(),
                self._medFileDataContext.getMeshName(),
                dt,
                it,
                self._medFileUMesh.getCoords().getInfoOnComponents(),
                vmin,
                vmax + 1,
            )
            idsInLoadedArrays = arrayOfGlobalIdsOfOrphanNodes - vmin
            self._orphan_coords = dataFromFile[0][idsInLoadedArrays]
            self._orphan_global_ids = arrayOfGlobalIdsOfOrphanNodes
            self._orphan_fam_ids = None
            if dataFromFile[2]:
                self._orphan_fam_ids = dataFromFile[2][idsInLoadedArrays]
            self._orphan_num_ids = None
            if dataFromFile[3]:
                self._orphan_num_ids = dataFromFile[3][idsInLoadedArrays]

    def dealWithOrphanNodes(self, meshResult):
        """
        meshResult : MEDFileUMesh instance to be enriched by orphan nodes if any
        """
        if self._orphan_coords.empty():
            return
        logger.debug("ajout des noeuds orphelins")
        coo = meshResult.getCoords()
        fam_coo = meshResult.getFamilyFieldAtLevel(1)
        glob_coo = meshResult.getGlobalNumFieldAtLevel(1)
        coo = mc.DataArrayDouble.Aggregate([coo, self._orphan_coords])
        if fam_coo:
            fam_coo = mc.DataArrayInt.Aggregate([fam_coo, self._orphan_fam_ids])
        glob_coo = mc.DataArrayInt.Aggregate([glob_coo, self._orphan_global_ids])
        meshResult.setCoords(coo)
        if fam_coo:
            meshResult.setFamilyFieldArr(1, fam_coo)
        meshResult.setGlobalNumFieldAtLevel(1, glob_coo)

    def getMeshDimension(self):
        return self.getMEDFileUMesh().getMeshDimension()

    def presenceOfPoint1(self):
        return self._presence_of_point1 is not None

    def getFamilyFieldAtLevelZero(self):
        return self._presence_of_point1.getFamilyFieldAtLevel(0)

    def getParaUMeshAtLevelZero(self):
        zeMesh = self._presence_of_point1[0]
        zeMesh.zipCoords()
        zeCells = self.cellsGen(-self.getMeshDimension())
        paramesh = mc.ParaUMesh(zeMesh, zeCells, self._fni_lev_0)
        return paramesh, self._presence_of_point1_fam_node_ids

    def getMEDFileUMesh(self):
        """return a MEDFileUMesh by loading a part of mesh named self._medFileDataContext.getMeshName()
        contained in file self._medFileDataContext.getFileName()
        Cells of the mesh are retrieved depending to the rank and size of MPI comm"""
        if hasattr(self, "_medFileUMesh"):
            return self._medFileUMesh
        rank = comm_world().rank
        size = comm_world().size
        ttps = mc.GetUMeshGlobalInfo(
            self._medFileDataContext.getFileName(), self._medFileDataContext.getMeshName()
        )[0]
        params = []
        cts = []
        for tps in ttps:
            for cellType, nbCellsType in tps:
                slc = mc.DataArray.GetSlice(slice(0, nbCellsType, 1), rank, size)
                params += [slc.start, slc.stop, slc.step]
                cts.append(cellType)
                pass
        mrs = mc.MEDFileMeshReadSelector()
        mrs.setNumberOfCoordsLoadSessions(10)
        self._medFileUMesh = mc.MEDFileUMesh.LoadPartOf(
            self._medFileDataContext.getFileName(),
            self._medFileDataContext.getMeshName(),
            cts,
            params,
            -1,
            -1,
            mrs,
        )

        self._presence_of_point1 = None
        # detect orphan nodes in parallel
        # si on a des MED_POINT1 on les retire.
        if -self._medFileUMesh.getMeshDimension() in self._medFileUMesh.getNonEmptyLevels():
            self._presence_of_point1 = mc.MEDFileUMesh()
            mc_mesh_lev_0 = self._medFileUMesh[-self.getMeshDimension()]
            fni_tmp = mc_mesh_lev_0.computeFetchedNodeIds()
            self._fni_lev_0 = fni_tmp.deepCopy()
            self._fni_lev_0.transformWithIndArr(self.getMEDFileUMesh().getPartDefAtLevel(1).toDAI())
            self._presence_of_point1[0] = mc_mesh_lev_0
            famFieldLev0 = self._medFileUMesh.getFamilyFieldAtLevel(-self.getMeshDimension())
            if famFieldLev0:
                self._presence_of_point1.setFamilyFieldArr(0, famFieldLev0)
            self._presence_of_point1_fam_node_ids = None
            if self._medFileUMesh.getFamilyFieldAtLevel(1):
                self._presence_of_point1_fam_node_ids = self._medFileUMesh.getFamilyFieldAtLevel(1)[
                    fni_tmp
                ]
            del self._medFileUMesh[-self.getMeshDimension()]
            self._medFileUMesh.zipCoords()
        self._fni = self.getMEDFileUMesh().getPartDefAtLevel(1).toDAI()

        # detect orphan nodes in parallel
        pda = mc.ParaDataArrayInt(self._fni)
        res = pda.buildComplement(self._medFileDataContext.getGlobalNumberOfNodes())
        if rank == 0:
            self.storeOrphanNodesIfAny(res)
        return self._medFileUMesh

    def getFetchedNodeIds(self):
        return self._fni

    def umeshAtLevel(self, lev):
        return self.getMEDFileUMesh()[lev]

    @property
    def umesh(self):
        """return the first MedCouplingUMesh contained in self.__getMEDFileUMesh()"""
        return self.umeshAtLevel(0)

    def getParaUMeshAtLevel(self, lev):
        zeMesh = self.umeshAtLevel(lev).deepCopy()
        zeCells = self.cellsGen(lev)
        zeNodes = self.getMEDFileUMesh().getPartDefAtLevel(1).toDAI()
        paramesh = mc.ParaUMesh(zeMesh, zeCells, zeNodes)
        nodeFamIds = None
        if self.getMEDFileUMesh().getFamilyFieldAtLevel(1):
            nodeFamIds = self.getMEDFileUMesh().getFamilyFieldAtLevel(1)
        return paramesh, nodeFamIds


class NeighborsOfNodes(metaclass=abc.ABCMeta):
    @abc.abstractproperty
    def nodes(self):
        """return nodes"""
        raise NotImplementedError("pure virtual, should be reimplemented")

    @abc.abstractproperty
    def neighborsSkyLineArray(self):
        """return node connectivities"""
        raise NotImplementedError("pure virtual, should be reimplemented")


def computeEnlargedNeighborsOfNodesLowMem(m):
    """
    Surcharge de mc.MEDCouplingUMesh.computeEnlargedNeighborsOfNodes qui reduit la conso m??moire
    """
    # n,ni = m.computeEnlargedNeighborsOfNodes()
    m_dup = m.deepCopyConnectivityOnly()
    fni_m = m_dup.computeFetchedNodeIds()
    ni3 = mc.DataArrayInt(m_dup.getNumberOfNodes())
    ni3[:] = 0
    m_dup.zipCoords()
    n2, ni2 = m_dup.computeEnlargedNeighborsOfNodes()
    if len(ni2) == 1:
        ni2 = mc.DataArrayInt(m.getNumberOfNodes() + 1)
        ni2[:] = 0
        # assert(n.isEqual(n2))
        # assert(ni.isEqual(ni2))
        return n2, ni2
    dsi = ni2.deltaShiftIndex()
    ni3[fni_m] = dsi
    ni3.computeOffsetsFull()
    n3 = fni_m[n2]
    # assert(n.isEqual(n3))
    # assert(ni.isEqual(ni3))
    return n3, ni3


class IncompleteNeighborsOfNodes(NeighborsOfNodes):
    """Retrieve neighbors of nodes (global numbering) for partialMedFileUMesh"""

    def __init__(self, partialMedFileUMesh):
        # nodes as fetched nodes
        self._nodes = partialMedFileUMesh.getFetchedNodeIds()
        sks = []
        for lev in partialMedFileUMesh._medFileUMesh.getNonEmptyLevels():
            neighbors, neighborsIdx = computeEnlargedNeighborsOfNodesLowMem(
                partialMedFileUMesh._medFileUMesh[lev]
            )
            neighbors.transformWithIndArr(self._nodes)
            sks.append(mc.MEDCouplingSkyLineArray(neighborsIdx, neighbors))
        sk = mc.MEDCouplingSkyLineArray.AggregatePacks(sks)
        sk.uniqueNotSortedByPack()
        neighborsIdx = sk.getIndexArray().buildUnique()
        """self._nodes_point1 = mc.DataArrayInt([])
        # on traite le cas particulier des MED_POINT1 connect?? a aucun autre node le buildUnique ci-dessus le vire du neighborhood et il est pourtant dans self._nodes
        if len(neighborsIdx)-1 != len(self._nodes):
            if -partialMedFileUMesh._medFileUMesh.getMeshDimension() not in partialMedFileUMesh._medFileUMesh.getNonEmptyLevels():
                raise RuntimeError("Mismatch between len(_nodes) and len(neighborsIdx)-1 in context extra MED_POINT1 ! Presence of cell with multiple same node id in nodal connectivity ?")
            tmp_umesh = partialMedFileUMesh._medFileUMesh.deepCopy()
            del tmp_umesh[-tmp_umesh.getMeshDimension()]
            nodesEff = tmp_umesh.computeFetchedNodeIds() + partialMedFileUMesh.nodalOffset
            self._nodes_point1 = self._nodes.buildSubstraction(nodesEff)
            self._nodes = nodesEff"""

        self._neighborsSkyLineArray = mc.MEDCouplingSkyLineArray(neighborsIdx, sk.getValuesArray())

    def getIdxByNode(self):
        """return mapping for self.nodes"""
        if hasattr(self, "_idxByNode"):
            return self._idxByNode
        self._idxByNode = self.nodes.invertArrayN2O2O2NOptimized()
        return self._idxByNode

    @property
    def nodes(self):
        return self._nodes

    @property
    def neighborsSkyLineArray(self):
        return self._neighborsSkyLineArray


class CompleteNeighborsOfNodes(NeighborsOfNodes):
    """Retrieve neighbors of nodes (global numbering) for partialMedFileUMesh + update neighbors from other procs"""

    def __init__(self, partialMedFileUMesh, globalNumberOfNodes):
        self.incomplete = IncompleteNeighborsOfNodes(partialMedFileUMesh)
        self._pska = mc.ParaSkyLineArray(
            self.incomplete.neighborsSkyLineArray, self.incomplete.nodes
        )
        self._pska = self._pska.equiRedistribute(globalNumberOfNodes)

    @property
    def nodes(self):
        return self._pska.getGlobalIdsArray()

    @property
    def neighborsSkyLineArray(self):
        return self._pska.getSkyLineArray()


class NodesPartition:
    def __init__(self, neighborsOfNodes, MyPartitioner):
        self._graph = MyPartitioner(neighborsOfNodes)
        self._graphnodes = neighborsOfNodes.nodes
        self.nodes

    def __getProcIdByIdx(self):
        """return a DataArrayInt containing all procIds by index"""
        if hasattr(self, "_procIdByIdx"):
            return self._procIdByIdx
        self._procIdByIdx = self._graph.get()
        return self._procIdByIdx

    def __getNodesByProc(self):
        """return a list of DataArrayInt containing all nodes for each proc"""
        if hasattr(self, "_nodesByProc"):
            return self._nodesByProc
        size = comm_world().size
        self._nodesByProc = [None] * size
        for procId in range(size):
            idx = self.__getProcIdByIdx().findIdsEqual(procId)
            self._nodesByProc[procId] = self._graphnodes[idx]
        return self._nodesByProc

    def __sendrecvNodesByProc(self):
        """receive and send nodes from/to other procs to update
        self._nodesByProc according to partitionning"""
        size = comm_world().size
        rank = comm_world().rank
        for procId in range(size):
            if procId == rank:
                continue
            self.__getNodesByProc()[procId] = comm_world().sendrecv(
                self.__getNodesByProc()[procId], source=procId, dest=procId
            )

    @property
    def nodes(self):
        if hasattr(self, "_nodes"):
            return self._nodes
        size = comm_world().size
        self.__sendrecvNodesByProc()
        self._nodes = mc.DataArrayInt([])
        for procId in range(size):
            self._nodes.aggregate(self.__getNodesByProc()[procId])
        return self._nodes


def AddGlobalFields(mm):
    mfd = mc.MEDFileData()
    mfd.setFields(mc.MEDFileFields())
    if mc.MEDCouplingSizeOfIDs() == 32:
        intField = mc.MEDCouplingFieldInt(mc.ON_NODES)
        intField.setMesh(mm[0])
        intField.setName("Numerotation Globale")
        intField.setArray(mm.getGlobalNumFieldAtLevel(1))
        fileIntFieldMultiTS = mc.MEDFileIntFieldMultiTS()
        fileIntFieldMultiTS = mc.MEDFileIntFieldMultiTS()
        fileIntFieldMultiTS.appendFieldNoProfileSBT(intField)
        mfd.getFields().pushField(fileIntFieldMultiTS)
    mfd.setMeshes(mc.MEDFileMeshes())
    mfd.getMeshes().pushMesh(mm)
    return mfd

    def __init__(self, medFileContext):
        self._medFileContext = medFileContext

    def __getMedCouplingUmesh(self):
        if hasattr(self, "_medCouplingUmesh"):
            return self._medCouplingUmesh
        self._medCouplingUmesh = self._medFileContext.getMEDFileMesh()[0]
        return self._medCouplingUmesh

    def __add(self, dataArrayInt, name):
        """add the globalField dataArrayInt with name into self.__getMedCouplingUmesh()"""
        intField = mc.MEDCouplingFieldInt(mc.ON_NODES)
        intField.setMesh(self.__getMedCouplingUmesh())
        intField.setName(name)
        intField.setArray(dataArrayInt)
        fileIntFieldMultiTS = mc.MEDFileIntFieldMultiTS()
        fileIntFieldMultiTS.appendFieldNoProfileSBT(intField)
        self._medFileContext.getMEDFileData().getFields().pushField(fileIntFieldMultiTS)

    def addNumGlo(self):
        nbNodes = self.__getMedCouplingUmesh().getNumberOfNodes()
        numGlob = mc.DataArrayInt(nbNodes)
        numGlob.iota()
        self.__add(numGlob, "Numerotation Globale")


class MedJoints:
    def __init__(self, mesh, intNodes):
        self._medFileData = mesh
        self._intNodes = intNodes

    def __getMeshName(self):
        """return self._meshName"""
        if hasattr(self, "_meshName"):
            return self._meshName
        self._meshName = self._medFileData.getName()
        return self._meshName

    def __getMedFileUmesh(self):
        """return self._medFileUmesh"""
        if hasattr(self, "_medFileUmesh"):
            return self._medFileUmesh
        return self._medFileData

    def __getMedCouplingUMesh(self):
        """return self._medCouplingUMesh"""
        if hasattr(self, "_medCouplingUMesh"):
            return self._medCouplingUMesh
        self._medCouplingUMesh = self.__getMedFileUmeshgetMeshes()[0]
        return self._medCouplingUMesh

    def __getIntNodes(self):
        """return a DataArrayInt containing internal nodes (global numbering)"""
        return self._intNodes

    def __getNodes(self):
        """return a DataArrayInt containing all nodes (global numbering)"""
        if hasattr(self, "_nodes"):
            return self._nodes
        self._nodes = self._medFileData.getGlobalNumFieldAtLevel(1)
        return self._nodes

    def __getExtNodes(self):
        """return a DataArrayInt containing external nodes (global numbering)"""
        if hasattr(self, "_extNodes"):
            return self._extNodes
        self._extNodes = self.__getNodes().buildSubstraction(self.__getIntNodes())
        return self._extNodes

    def __getNodesGloToLoc(self):
        """return the global to local mapping for self._nodes"""
        if hasattr(self, "_nodeGloToLoc"):
            return self._nodeGloToLoc
        self._nodeGloToLoc = self.__getNodes().invertArrayN2O2O2NOptimized()
        return self._nodeGloToLoc

    def __findLocNodesFromGlo(self, nodes):
        """find local nodes corresponding to global nodes"""
        nodes = nodes[:]
        nodes.transformWithIndArr(self.__getNodesGloToLoc())
        return nodes

    def __bcastExtNodes(self):
        """Communicate self.__getExtNodes() with other procs.
        Fill in self._extNodesByProc from other procs (global numbering)"""
        rank = comm_world().rank
        size = comm_world().size
        self._extNodesByProc = [None] * size
        for procId in range(size):
            if procId == rank:
                self._extNodesByProc[procId] = comm_world().bcast(self.__getExtNodes(), procId)
            else:
                self._extNodesByProc[procId] = comm_world().bcast(None, procId)

    def __computeIntNodesToSendByProc(self):
        """compute intersections between self._extNodesByProc and self.__getIntNodes()
        fill in self._intNodesToSendByProc (global numbering)"""
        self.__bcastExtNodes()
        size = comm_world().size
        rank = comm_world().rank
        self._intNodesToSendByProc = [None] * size
        for procId in range(size):
            if procId == rank:
                continue
            self._intNodesToSendByProc[procId] = self._extNodesByProc[procId].buildIntersection(
                self.__getIntNodes()
            )
        return self._intNodesToSendByProc

    def __getIntNodesToSendByProc(self):
        """return self._intNodesToSendByProc (global numbering)"""
        if hasattr(self, "_intNodesToSendByProc"):
            return self._intNodesToSendByProc
        self.__computeIntNodesToSendByProc()
        return self._intNodesToSendByProc

    def __getIntLocNodesToSendByProc(self):
        """return self._intLocNodesToSendByProc (local numbering) corresponding to self._intNodesToSendByProc"""
        if hasattr(self, "_intLocNodesToSendByProc"):
            return self._intLocNodesToSendByProc
        size = comm_world().size
        rank = comm_world().rank
        self._intLocNodesToSendByProc = [None] * size
        for procId in range(size):
            if procId == rank:
                continue
            self._intLocNodesToSendByProc[procId] = self.__findLocNodesFromGlo(
                self.__getIntNodesToSendByProc()[procId]
            )
        return self._intNodesToSendByProc

    def __sendIntNodesAndRecvExtNodesByProc(self):
        """Send int nodes contained in self._intNodesToSendByProc to other procs
        Recv ext nodes from other procs and store then into self._extNodesRecvByProc"""
        size = comm_world().size
        rank = comm_world().rank
        self._extNodesRecvByProc = [None] * size
        for procId in range(size):
            if procId == rank:
                continue
            self._extNodesRecvByProc[procId] = comm_world().sendrecv(
                self.__getIntNodesToSendByProc()[procId], source=procId, dest=procId
            )

    def __sendIntLocNodesAndRecvExtLocNodesByProc(self):
        """Send int nodes contained in self._intLocNodesToSendByProc to other procs
        Recv ext nodes from other procs and store them into self._extLocNodesRecvByProc"""
        size = comm_world().size
        rank = comm_world().rank
        self._extLocNodesRecvByProc = [None] * size
        for procId in range(size):
            if procId == rank:
                continue
            self._extLocNodesRecvByProc[procId] = comm_world().sendrecv(
                self.__getIntLocNodesToSendByProc()[procId], source=procId, dest=procId
            )

    def __getExtNodesRecvByProc(self):
        """return self._extNodesRecvByProc (global numbering)"""
        if hasattr(self, "_extNodesRecvByProc"):
            return self._extNodesRecvByProc
        self.__getIntNodesToSendByProc()
        self.__sendIntNodesAndRecvExtNodesByProc()
        return self._extNodesRecvByProc

    def __getExtLocNodesRecvByProc(self):
        """return self._extLocNodesRecvByProc (local numbering of other procs)"""
        if hasattr(self, "_extLocNodesRecvByProc"):
            return self._extLocNodesRecvByProc
        self.__getIntLocNodesToSendByProc()
        self.__sendIntLocNodesAndRecvExtLocNodesByProc()
        return self._extLocNodesRecvByProc

    def __getExtLocNodesToSendByProc(self):
        """return self._extLocNodesToSendByProc (local numbering of this procs)"""
        if hasattr(self, "_extLocNodesToSendByProc"):
            return self._extLocNodesToSendByProc
        self.__getExtNodesRecvByProc()
        size = comm_world().size
        rank = comm_world().rank
        self._extLocNodesToSendByProc = [None] * size
        for procId in range(size):
            if procId == rank:
                continue
            self._extLocNodesToSendByProc[procId] = self.__findLocNodesFromGlo(
                self.__getExtNodesRecvByProc()[procId]
            )
        return self._extLocNodesToSendByProc

    def __sendExtLocNodesAndRecvIntLocNodesByProc(self):
        """Send ext nodes contained in self._extLocNodesToSendByProc to other procs
        Recv int nodes from other procs and store them into self._intLocNodesRecvByProc"""
        size = comm_world().size
        rank = comm_world().rank
        self._intLocNodesRecvByProc = [None] * size
        for procId in range(size):
            if procId == rank:
                continue
            self._intLocNodesRecvByProc[procId] = comm_world().sendrecv(
                self.__getExtLocNodesToSendByProc()[procId], source=procId, dest=procId
            )

    def __getIntLocNodesRecvByProc(self):
        """return self._intLocNodesRecvByProc (local numbering of other procs)"""
        if hasattr(self, "_intLocNodesRecvByProc"):
            return self._intLocNodesRecvByProc
        self.__getExtLocNodesToSendByProc()
        self.__sendExtLocNodesAndRecvIntLocNodesByProc()
        return self._intLocNodesRecvByProc

    def __push(self, jointName, procId, thisLocNodes, otherLocNodes):
        """add a joint into self._medFileJoints"""
        if not hasattr(self, "_medFileJoints"):
            self._medFileJoints = mc.MEDFileJoints()
        correspondence = mc.DataArrayInt(0, 2)
        correspondence.aggregate(mc.DataArrayInt.Meld(thisLocNodes + 1, otherLocNodes + 1))
        medFileJoint = mc.MEDFileJoint(
            jointName, self.__getMeshName(), self.__getMeshName(), procId
        )
        medFileJointOneStep = mc.MEDFileJointOneStep()
        medFileJointOneStep.pushCorrespondence(mc.MEDFileJointCorrespondence(correspondence[:]))
        medFileJoint.pushStep(medFileJointOneStep)
        self._medFileJoints.pushJoint(medFileJoint)

    def __set(self):
        if hasattr(self, "_medFileJoints"):
            self.__getMedFileUmesh().setJoints(self._medFileJoints)

    def write(self):
        size = comm_world().size
        rank = comm_world().rank
        for procId in range(size):
            if procId == rank:
                continue
            otherLocNodes = self.__getExtLocNodesRecvByProc()[procId]
            thisLocNodes = self.__getExtLocNodesToSendByProc()[procId]
            if not otherLocNodes.empty():
                jointName = (str(rank) + " " + str(procId)).strip()
                self.__push(jointName, procId, thisLocNodes, otherLocNodes)
            otherLocNodes = self.__getIntLocNodesRecvByProc()[procId]
            thisLocNodes = self.__getIntLocNodesToSendByProc()[procId]
            if not otherLocNodes.empty():
                jointName = (str(procId) + " " + str(rank)).strip()
                self.__push(jointName, procId, thisLocNodes, otherLocNodes)
        self.__set()

    def setGlobalNumField(self):
        self.__getMedFileUmesh().setGlobalNumFieldAtLevel(1, self.__getNodes())


def GetUMeshFrom(paramesh, cellFamIds, nodeFamIds, cellsPartition):
    with ChronoCtxMgDbg("paramesh.redistributeCells"):
        pm2 = paramesh.redistributeCells(cellsPartition)
    mergeMesh = mc.MEDFileUMesh()
    mergeMesh[0] = pm2.getMesh()
    if cellFamIds is not None:
        with ChronoCtxMgDbg("paramesh.redistributeCellField"):
            arr = paramesh.redistributeCellField(cellsPartition, cellFamIds)
        mergeMesh.setFamilyFieldArr(0, arr)
    if nodeFamIds is not None:
        with ChronoCtxMgDbg("paramesh.redistributeNodeField"):
            arrn = paramesh.redistributeNodeField(cellsPartition, nodeFamIds)
        mergeMesh.setFamilyFieldArr(1, arrn)
    mergeMesh.setGlobalNumFieldAtLevel(1, pm2.getGlobalNodeIds())
    return mergeMesh


def FuseLevels(dico):
    """
    dico : dictionnaire dont les clefs sont les niveaux et les valeurs les instance de MEDFileUMesh associ??es
    Cette methode prend tous MEDFileUMesh de tous les niveaux et le merge dans un MEDFileUMesh qui sera retourn??.
    Particularit?? de cette m??thode : les coordonn??es des niveaux autres que 0 ne sont pas inclus forc??ment dans le level 0.
    """
    globalNodeIds = []
    famNodes = []
    coords = []
    for _, mg in dico.items():
        famNode = mg.getFamilyFieldAtLevel(1)
        if famNode:
            famNodes.append(famNode)
        globalNodeIds.append(mg.getGlobalNumFieldAtLevel(1))
        coords.append(mg.getCoords())

    aggregatedNodeIds = mc.DataArrayInt.Aggregate(globalNodeIds)
    aggregatedNodeIdsSorted = aggregatedNodeIds.copySorted()
    nodeIdsIntoAggregatedIds = mc.DataArrayInt.FindPermutationFromFirstToSecondDuplicate(
        aggregatedNodeIdsSorted, aggregatedNodeIds
    )
    idxOfSameNodeIds = aggregatedNodeIdsSorted.indexOfSameConsecutiveValueGroups()
    n2o_nodes = nodeIdsIntoAggregatedIds[idxOfSameNodeIds[:-1]]
    finalGlobalNodeIds = aggregatedNodeIdsSorted[idxOfSameNodeIds[:-1]]

    finalCoords = mc.DataArrayDouble.Aggregate(coords)[n2o_nodes]

    ret = mc.MEDFileUMesh()
    for lev, mg in dico.items():
        current = mg.getGlobalNumFieldAtLevel(1)
        aa = finalGlobalNodeIds.findIdForEach(current)
        m = mg[0].deepCopy()
        m.renumberNodesInConn(aa)
        m.setCoords(finalCoords)
        ret[lev] = m
        famCellField = mg.getFamilyFieldAtLevel(0)
        if famCellField:
            ret.setFamilyFieldArr(lev, famCellField)

    if any(famNodes):
        ret.setFamilyFieldArr(1, mc.DataArrayInt.Aggregate(famNodes)[n2o_nodes])

    ret.setGlobalNumFieldAtLevel(1, finalGlobalNodeIds)
    return ret


def MakeThePartition(fileName, meshName, MyPartitioner):
    medFileContext = MedFileContext(fileName, meshName)

    with ChronoCtxMg("lecture partielle du maillage"):
        partialMedFileUMesh = PartialMedFileUMesh(medFileContext)

    with ChronoCtxMg("r??cup??ration des voisins du maillage partiel"):
        neighborsOfNodes = CompleteNeighborsOfNodes(
            partialMedFileUMesh, medFileContext.getGlobalNumberOfNodes()
        )

    with ChronoCtxMg("cr??ation de la partition via ptscotch"):
        nodesPartition = NodesPartition(neighborsOfNodes, MyPartitioner)

    logger.debug("r??cup??ration des elements de la partition")

    # on rajoute tous les levels
    mm_levs = {}
    for lev in partialMedFileUMesh.getMEDFileUMesh().getNonEmptyLevels():
        pm_lev, fam_field_node = partialMedFileUMesh.getParaUMeshAtLevel(lev)
        if lev == 0:
            with ChronoCtxMgDbg("getCellIdsLyingOnNodes False"):
                cellsPartition_lev = pm_lev.getCellIdsLyingOnNodes(nodesPartition.nodes, False)
        else:
            with ChronoCtxMgDbg("getCellIdsLyingOnNodes True"):
                cellsPartition_lev = pm_lev.getCellIdsLyingOnNodes(globalNodeIdsOfRk0, True)
        #  connaissant les cells du niveau lev dont j ai besoin j appelle la communaut?? pour me donner mon morceau
        mm_lev = GetUMeshFrom(
            pm_lev,
            partialMedFileUMesh.getMEDFileUMesh().getFamilyFieldAtLevel(lev),
            fam_field_node,
            cellsPartition_lev,
        )
        #
        if lev == 0:
            globalNodeIdsOfRk0 = mm_lev.getGlobalNumFieldAtLevel(1)[
                mm_lev[0].computeFetchedNodeIds()
            ]
            globalNodeIdsOfRk0 = globalNodeIdsOfRk0.buildUnion(nodesPartition.nodes)
        #
        mm_levs[lev] = mm_lev
        pass

    if partialMedFileUMesh.presenceOfPoint1():
        pm_lev, fam_field_node = partialMedFileUMesh.getParaUMeshAtLevelZero()
        with ChronoCtxMgDbg("getCellIdsLyingOnNodes False : presence of point 1"):
            cellsPartition_lev = pm_lev.getCellIdsLyingOnNodes(nodesPartition.nodes, False)
        mm_lev = GetUMeshFrom(
            pm_lev,
            partialMedFileUMesh.getFamilyFieldAtLevelZero(),
            fam_field_node,
            cellsPartition_lev,
        )
        mm_levs[-partialMedFileUMesh.getMeshDimension()] = mm_lev

    meshResult = FuseLevels(mm_levs)
    # on rajoute les noms des familles/groupes attach??s aux ids

    famsPy, grpsPy = RetrieveFamGrpsMap(medFileContext)
    for fam, famid in famsPy.items():
        meshResult.addFamily(fam, famid)
        meshResult.setGroupsOnFamily(fam, grpsPy[fam])
        pass

    logger.debug("ajout des raccords")
    medJoints = MedJoints(meshResult, nodesPartition.nodes)
    medJoints.setGlobalNumField()
    medJoints.write()

    # on rajoute sur le noeud 0 les noeuds orphelins s'il y en avait

    if comm_world().rank == 0:
        partialMedFileUMesh.dealWithOrphanNodes(meshResult)

    # on ajoute les champs globaux aux maillage et on met le tout dans un objet MEDFileData

    mfd = AddGlobalFields(meshResult)
    return mfd


class GraphPartitionerCompute:
    def __init__(self, neighborsOfNodes):
        self._graph = mc.MEDPartitioner.Graph(
            neighborsOfNodes.neighborsSkyLineArray, mc.Graph.PTSCOTCH, None, neighborsOfNodes.nodes
        )

    def get(self):
        self._graph.partGraph(comm_world().size)
        return self._graph.getPartition().getValuesArray()


def GraphPartitionerReload(graphName):
    class GraphPartitionerReload:
        """Pour tester avec une decomposition donn??e"""

        def __init__(self, neighborsOfNodes):
            self.nodes = neighborsOfNodes.nodes

        def get(self):
            ids = mc.DataArrayInt(
                np.array(
                    np.load(graphName)["arr_0"],
                    dtype=eval("np.int{}".format(mc.MEDCouplingSizeOfIDs())),
                )
            )
            return ids[self.nodes]

    return GraphPartitionerReload


def GetGraphPartitioner(graph_scotch):
    if not graph_scotch:
        return GraphPartitionerCompute
    else:
        return GraphPartitionerReload(graph_scotch)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", dest="fileName", required=True, help="med file path")
    parser.add_argument(
        "-m",
        "--meshname",
        dest="meshName",
        default=None,
        help="mesh contained in filepath to be considered",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbosity",
        type=int,
        default=1,
        help="verbosity. Default 1 (print info). 0 only errors reported. 2 and higher all debug messages",
    )
    parser.add_argument(
        "-g",
        "--graph_scotch",
        dest="graphScotch",
        default=None,
        help="Partitioning file ( numpy format giving for each node its target rank ). It allows to shortcut PTSCOTCH computation",
    )

    args = parser.parse_args()

    setVerbose(args.verbosity)

    if comm_world().size == 1:
        print("using 1 cpu, nothing to do")
        sys.exit(0)

    with MasterChronoCtxMg("Decoupage"):

        if StrictVersion(mc.MEDFileVersionOfFileStr(args.fileName)) < StrictVersion("3.0.0"):
            raise RuntimeError(
                'File "{}" is < 3.0.0 : Not compatible with // load. To use this // decoupeur, just convert your old MED file into MED file version > 3.0.0, to do that ConvertMEDFileTo33.py MEDCoupling tool may help you.'
            )

        mfd = MakeThePartition(args.fileName, args.meshName, GetGraphPartitioner(args.graphScotch))

        with ChronoCtxMg("ecriture sur disque"):
            mfd.write33(BuildPartNameFromOrig(args.fileName, comm_world().rank), 2)

        comm_world().barrier()

    pass
