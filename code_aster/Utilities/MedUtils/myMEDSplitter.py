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
# --------------------------------------------------------------------

import argparse
import os
import sys
from collections import defaultdict, deque
from sys import stdout
from time import *

import numpy as np

from medcoupling import *
from MEDLoaderSplitter import MEDLoaderSplitter

GLO_NUM_FIELD_NAME = "Numerotation Globale"

def versiontuple(v):
    return tuple(map(int, (v.split("."))))

def imprimerTemps():
    """
    Utilitaire d'impression du temps
    """
    tmp = localtime()
    heures = str(tmp[3])
    if len(heures) == 1: heures = "0"+heures
    minutes = str(tmp[4])
    if len(minutes) == 1: minutes = "0"+minutes
    secondes = str(tmp[5])
    if len(secondes) == 1: secondes = "0"+secondes
    txt = heures+" h "+minutes+" min "+secondes+" s"
    print(txt)

a = versiontuple(str(MEDCouplingVersionStr()))
b = versiontuple("6.6.0")
if not a > b:
    print("Mauvaise version de MEDCoupling", MEDCouplingVersionStr())
    assert False

class MyMedSplitter:
    def __ecritureMaillages(self, fichierMED, maillage, listPartitionsMailles, grpsNoeuds):

        imprimerTemps()
        print("Création d'une numérotation globale")

        champNumGlobal = MEDCouplingFieldID(ON_NODES)
        champNumGlobal.setMesh(maillage)
        champNumGlobal.setName(GLO_NUM_FIELD_NAME)
        nbNodes = maillage.getNumberOfNodes()
        valGlob = DataArrayInt(nbNodes) ; valGlob.iota()
        champNumGlobal.setArray(valGlob)
        f = MEDFileIDFieldMultiTS()
        f.appendFieldNoProfileSBT(champNumGlobal)
        fichierMED.getFields().pushField(f)

        imprimerTemps()
        print("Écriture des maillages découpés")
        listeNomsGrpNoeuds = fichierMED.getMeshes()[0].getGroupsOnSpecifiedLev(1)
        ancienGrps = []
        for nom in listeNomsGrpNoeuds:
            groupe = fichierMED.getMeshes()[0].getNodeGroupArr(nom, False)
            ancienGrps.append(groupe)
        fichierMED.getMeshes()[0].setGroupsAtLevel(1, ancienGrps + grpsNoeuds, False)

        loaderAndSplitter = MEDLoaderSplitter(fichierMED, listPartitionsMailles)
        fichierDecoupes = loaderAndSplitter.getSplittedInstances()
        noeudsFrontiere = []
        maillages = []

        for numero, mfdp in enumerate(fichierDecoupes):
            # MEDCouplingUMesh
            maillageCourant = mfdp.getMeshes()[0][0]
            # MEDFileUmesh
            maillageCourant2 = mfdp.getMeshes()[0]

            maillages.append(maillageCourant)

            noeudsInt = maillageCourant2.getGroupArr(1, "EXT_" + str(numero))
            nbNoeudsTot = maillageCourant.getNumberOfNodes()

            frontiere = noeudsInt.buildComplement(nbNoeudsTot)
            noeudsFrontiere.append(frontiere)

        return maillages, noeudsFrontiere, fichierDecoupes

    def __ecritureRaccords(self, maillages, noeudsFrontiere, procIdOnNodes, fichierDecoupes):
        nbPart = len(maillages)

        imprimerTemps()
        print("Création des correspondances numérotations locales/globales")

        nomMaillageLocalEtDistant = maillages[0].getName()

        # nsellenet
        #import cProfile, pstats, StringIO
        #pr = cProfile.Profile()
        #pr.enable()
        # ... do something ...
        listLocGlo = []
        listGloLoc = []
        for proc in range(nbPart):
            print(proc)
            fichierMED = fichierDecoupes[proc]
            f1ts=fichierMED.getFields()[GLO_NUM_FIELD_NAME][(-1,-1)]
            curField = f1ts.field(fichierMED.getMeshes()[f1ts.getMeshName()])
            valeur = curField.getArray()

            #maillageCourant = maillages[proc]
            #nbNoeuds = maillageCourant.getNumberOfNodes()
            dictGloLoc = valeur.invertArrayN2O2O2NOptimized()
            mm=fichierMED.getMeshes()[0]
            mm.setGlobalNumFieldAtLevel(1,valeur)
            listLocGlo.append(valeur[:])
            listGloLoc.append(dictGloLoc)
        #
        #pr.disable()
        #s = StringIO.StringIO()
        #sortby = 'cumulative'
        #ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        #ps.print_stats()
        #print s.getvalue()
        # nsellenet

        imprimerTemps()
        print("Détermination des raccords")
        correspondances = [DataArrayInt(0,2) for i in range(nbPart*nbPart)]

        for proc1 in range(nbPart):
            # maillageCourant1 = maillages[proc1]
            listeNoeuds1 = noeudsFrontiere[proc1]
            dictLocGlo = listLocGlo[proc1]
            #
            numGlo = dictLocGlo[listeNoeuds1]
            procs2 = procIdOnNodes[numGlo]
            for proc2 in procs2.getDifferentValues().getValues():
                ids=procs2.findIdsEqual(proc2)
                globLoc=listGloLoc[proc2]
                numLoc=numGlo[ids] ; numLoc.transformWithIndArr(globLoc)
                idNoeud=listeNoeuds1[ids]
                correspondances[nbPart*proc1 + proc2].aggregate(DataArrayInt.Meld(idNoeud+1,numLoc+1))
                pass
            stdout.write("\r\t%d processeurs terminés" % (proc1 + 1))
            stdout.flush()
        stdout.write("\n")

        imprimerTemps()
        print("Écriture des raccords")
        js=[MEDFileJoints() for i in range(nbPart)]
        format = "%-4d"
        for proc1 in range(nbPart):
            for proc2 in range(nbPart):
                raccord = correspondances[nbPart*proc1 + proc2]
                if not raccord.empty():
                    assert proc1 != proc2
                    raccord2 = raccord[:]
                    raccord3 = raccord[:,[1,0]]

                    nomRaccord = format%proc1 + " " + format%proc2
                    nomRaccord = nomRaccord.strip()
                    j1=MEDFileJoint(nomRaccord,nomMaillageLocalEtDistant,nomMaillageLocalEtDistant,proc2)
                    j1_p=MEDFileJointOneStep()
                    j1_p.pushCorrespondence(MEDFileJointCorrespondence(raccord2))
                    j1.pushStep(j1_p)
                    js[proc1].pushJoint(j1)
                    #
                    j2=MEDFileJoint(nomRaccord,nomMaillageLocalEtDistant,nomMaillageLocalEtDistant,proc1)
                    j2_p=MEDFileJointOneStep()
                    j2_p.pushCorrespondence(MEDFileJointCorrespondence(raccord3))
                    j2.pushStep(j2_p)
                    js[proc2].pushJoint(j2)
            stdout.write("\r\t%d processeurs terminés" % (proc1 + 1))
            stdout.flush()
        stdout.write("\n")
        for proc in range(nbPart):
            fichierDecoupes[proc].getMeshes()[nomMaillageLocalEtDistant].setJoints(js[proc])
            pass

    def __partitionnementMaillage(self, maillage, nbPart, nomFichierGraph):
        imprimerTemps()
        print("Création de la connectivité")

        # nbNoeudsTot = maillage.getNumberOfNodes()

        if not nomFichierGraph:
            edgetabArray,verttabArray=maillage.computeEnlargedNeighborsOfNodes()
            #
            sk=MEDCouplingSkyLineArray()
            sk.set(verttabArray,edgetabArray)
            g=MEDPartitioner.Graph(sk,1)
            g.partGraph(nbPart)
            procIdOnNodes=g.getPartition().getValuesArray()
        else:
            procIdOnNodes=DataArrayInt( np.array(np.load(nomFichierGraph)["arr_0"],dtype=eval("np.int{}".format(MEDCouplingSizeOfIDs()))) )
        #
        print("Création des partitions")
        #######
        # nbMailles = maillage.getNumberOfCells()
        grpsNoeuds=[procIdOnNodes.findIdsEqual(i) for i in range(nbPart)]
        for i,tn in zip(range(nbPart),grpsNoeuds):
            tn.setName("EXT_" + str(i))
        stdout.write("\r\t%d %%" % 0)
        stdout.flush()
        listPartitionsMailles=[maillage.getCellIdsLyingOnNodes(gn,False) for gn in grpsNoeuds]
        return listPartitionsMailles, grpsNoeuds, procIdOnNodes

    def decoupageMaillageEnMemoire(self, fichierMED, maillage, nbPartition, nomScotch = None):
        listPartitionsMailles, grpsNoeuds, procIdOnNodes = self.__partitionnementMaillage(maillage, nbPartition, nomScotch)
        maillages, noeudsFrontiere,fichierDecoupes       = self.__ecritureMaillages(fichierMED, maillage,
                                                                                    listPartitionsMailles,
                                                                                    grpsNoeuds)
        self.__ecritureRaccords(maillages, noeudsFrontiere, procIdOnNodes, fichierDecoupes)
        return fichierDecoupes

    @classmethod
    def BuildFileNameOfPart(cls,nomFichier,procId):
        nomSansExtension = os.path.splitext(nomFichier)[0]
        return "%s_%i.med"%(nomSansExtension, procId)
        #return str(procId)+".med"

    def decoupageMaillage(self, nomFichier, nomMaillage, nbPartition, nomScotch = None):
        print("Lecture du maillage", nomMaillage, "dans le fichier", nomFichier, "")
        fichierMED = MEDFileData(nomFichier)
        tmpMEDFileMesh = fichierMED.getMeshes()[nomMaillage]
        tmpCouplingMesh = tmpMEDFileMesh[0]
        # nsellenet
        maillage = tmpCouplingMesh.buildUnstructured()
        ######
        imprimerTemps()
        fichierDecoupes=self.decoupageMaillageEnMemoire(fichierMED,maillage,nbPartition,nomScotch)
        ###### Writing
        for proc,mfd in enumerate(fichierDecoupes):
            mfd.write(MyMedSplitter.BuildFileNameOfPart(nomFichier,proc),2)
        imprimerTemps()

if __name__ == '__main__':

    assert sys.version[:3] > '3.0'

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fichier", help = "Nom du fichier")
    parser.add_argument("-m", "--maillage", help = "Nom du maillage")
    parser.add_argument("-p", "--partitions", help = "Nombre de domaines", type = int)
    parser.add_argument("--graph_scotch", help = "Nom du fichier contenant le partitionement (format numpy qui pour chaque noeud donne le rank). Permet de court cicuiter SCOTCH")

    mesOpts = parser.parse_args()
    nomFichierMED = mesOpts.fichier
    nomMaillage = mesOpts.maillage
    nbPart = mesOpts.partitions
    nomScotch = mesOpts.graph_scotch

    assert nomFichierMED != None and nomFichierMED != ""
    assert nomMaillage != None and nomMaillage != ""
    assert nbPart > 1

    monMedSplitter = MyMedSplitter()
    monMedSplitter.decoupageMaillage(nomFichierMED, nomMaillage, nbPart, nomScotch)
    pass
