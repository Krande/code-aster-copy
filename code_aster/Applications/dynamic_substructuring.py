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

import numpy as np
import scipy.sparse
import scipy.sparse.linalg
import os
from ..Objects.user_extensions import WithEmbeddedObjects
try :
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from ..Commands import (CREA_CHAMP, RESOUDRE, CREA_RESU, FACTORISER)



class Interface(WithEmbeddedObjects):
    """ Interface object for dynamic substructring.

    The interface is a 2D surface, modeled by a group of cells.
    It needs not to be conforming  : the mesh of each substructure on the interface
    needs not to be the same.

    The interface is fully defined by :
    - 2 substructures
    - the name of the group of cells modeling the interface
    => Each substructure must have a group of cells with the same name (the name of the interface)
    """

    # Mandatory class attribute in order to manage aster serialization process
    aster_embedded = ["sub1", "sub2"]

    def __init__(self, sub1, sub2, interfaceName: str):
        self.sub1 = sub1
        self.sub2 = sub2
        self.name = interfaceName

        if not (self.sub1.mesh.hasGroupOfCells(interfaceName) and
                self.sub2.mesh.hasGroupOfCells(interfaceName)):
            raise RuntimeError(
                "Each substructure must have a group of cells named {}".format(interfaceName))

    def __getstate__(self):
        """Method for aster serialization process
        """
        return [self.sub1, self.sub2, self.name]

    def __setstate__(self, state):
        """Method for aster serialization process
        """
        assert len(state) == 3, state
        (self.sub1, self.sub2, self.name) = state

    def computeInterfaceDofs(self, dofsType: str):
        """
        Compute interface displacements / pressure used for static corrections
        Usefull for component mode synthesis (see Craig & Chang)

        Input :
        ------
        dofsType : either 'Struct', 'Fluid' or 'IFS'
        """
        # check the value of the argument
        VALID_DOFS = {'Structure', 'Fluid', 'IFS'}
        if dofsType not in VALID_DOFS:
            raise ValueError("dofsType must be one of %r." % VALID_DOFS)
        # we begin by merging all nodes of the interface in common objects
        nodes = np.hstack((self.sub1.nodes, self.sub2.nodes))
        coords = np.vstack((self.sub1.coords, self.sub2.coords))
        if not self.sub1.mesh.hasGroupOfNodes(self.name):
            connect = self.sub1.mesh.getConnectivity()
            lCells = np.array(self.sub1.mesh.getCells(self.name)) - 1
            lINodes = []
            for cell in lCells:
                lINodes += connect[cell]
            interfaceNodes1 = list(np.unique(np.array(lINodes)))
        else:
            interfaceNodes1 = self.sub1.mesh.getNodes(self.name)
        if not self.sub2.mesh.hasGroupOfNodes(self.name):
            connect = self.sub2.mesh.getConnectivity()
            lCells = np.array(self.sub2.mesh.getCells(self.name)) - 1
            lINodes = []
            for cell in lCells:
                lINodes += connect[cell]
            interfaceNodes2 = list(np.unique(np.array(lINodes)))
        else:
            interfaceNodes2 = self.sub2.mesh.getNodes(self.name)
        interfaceNodes = np.array(interfaceNodes1 + interfaceNodes2)
        nInterfaceNodes = len(interfaceNodes)

        # we get the indices of the nodes of the interface in the common list
        ind_no1 = np.where(np.in1d(nodes[:len(self.sub1.nodes)], interfaceNodes1))[0]
        ind_no2 = np.where(np.in1d(nodes[len(self.sub1.nodes):], interfaceNodes2))[0]
        # since the indices are for a sub-list, we must shift them for use in the common list
        ind_no2 += len(self.sub1.nodes)
        ind_no = np.hstack((ind_no1, ind_no2))

        # evaluate the mean plane
        Coord_t = coords[ind_no, :]
        for i1 in range(3):
            Coord_t[:, i1] -= np.mean(Coord_t[:, i1])

        # if there is an ovious mean plane, we use it. Otherwise we compute it with a SVD
        CoordSum = np.abs(Coord_t).sum(axis=0)
        if np.any(CoordSum<1.e-12):
            indx = np.where(CoordSum<1.e-12)[0][0]
            U = np.eye(3)
            U[:,[2, indx]] = U[:,[indx,2]]  # swap the columns
        else:
            U, _, _ = np.linalg.svd(Coord_t.T)

        # project onto the mean plane
        NewCoord = np.matmul(U.T, Coord_t.T)

        # place the centre in [0,0,0]
        Xmax = np.max(NewCoord[0, :])
        Ymax = np.max(NewCoord[1, :])
        Xmin = np.min(NewCoord[0, :])
        Ymin = np.min(NewCoord[1, :])

        NewCoord[0, :] -= (Xmax + Xmin) / 2.
        NewCoord[1, :] -= (Ymax + Ymin) / 2.

        R = np.sqrt(NewCoord[0, :]**2 + NewCoord[1, :]**2)
        Th = np.arctan2(NewCoord[1, :], NewCoord[0, :])

        Rmax = np.max(R)
        Interf = []

        #  section for fluid dofs
        if dofsType == 'Fluid' or dofsType == 'IFS':
            """
            Interface pressure fields
            * constant pressure
            * linear pressure along local X vector
            * linear pressure along local Y vector
            * quadratic field - 1 at the center / 0 on edges

            """
            Pres = np.zeros((nInterfaceNodes, 3))
            Pres[:, 0] = np.ones(nInterfaceNodes)  #  constant pressure
            Pres[:, 1] = NewCoord[0, :] / Rmax  # linear pressure
            Pres[:, 2] = NewCoord[1, :] / Rmax  # linear pressure
            # Pres[:,3] = 1.- (R/Rmax)**2  # quadratic field - not so usefull

            Pres *= 1.e5  #  Normalize

            for i1 in range(3):
                Interf.append(np.zeros((nInterfaceNodes, 4)))
                Interf[-1][:, 3] += Pres[:, i1]

#  section for structure dofs
        if dofsType == 'Structure' or dofsType == 'IFS':
            """
            Interface displacement fields
            * 6 rigid body modes
            * 1 inflating mode
            * 2 ovalisation modes

            """

            for i1 in range(11):
                Interf.append(np.zeros((nInterfaceNodes, 4)))
                Disp = np.zeros((nInterfaceNodes, 3))

                if i1 < 3:  # -- translation rigid body modes
                    Disp[:, i1] += 1.

                elif (i1 >= 3) & (i1 < 5):  # -- rotation along X / Y local vectors
                    Disp = np.zeros((nInterfaceNodes, 3))
                    Disp[:, 2] = R*np.cos(Th + (i1-3)*np.pi/2)/Rmax

                elif i1 == 5:  # -- rotation along Z local vectors
                    Disp[:, 0] = NewCoord[1, :] / Rmax
                    Disp[:, 1] = -NewCoord[0, :] / Rmax

                elif i1 == 6:  # -- inflation
                    Disp[:, 0] = np.cos(Th)
                    Disp[:, 1] = np.sin(Th)

                elif (i1 < 11):  # -- ovalisation
                    Disp[:, 0] = np.cos(Th)*np.sin(2*Th + (i1-7)*np.pi/2)
                    Disp[:, 1] = np.sin(Th)*np.sin(2*Th + (i1-7)*np.pi/2)
                # elif (i1 == 7):  # -- ovalisation1
                #     Disp[:, 0] = np.cos(Th)*1.5
                #     Disp[:, 1] = np.sin(Th)*0.2

                # elif (i1 == 8):  # -- ovalisation2
                #     Disp[:, 0] = np.cos(Th+np.pi/2)*1.5
                #     Disp[:, 1] = np.sin(Th+np.pi/2)*0.2
                else:
                    print('not done yet...')

                Disp = np.matmul(Disp, U.T)
                for k1 in range(3):
                    Interf[-1][:, k1] += Disp[:, k1]

        # transfer the interface data to the specific substructure
        self.sub1.setInterfaceData(self.name, interfaceNodes1, [
                                   a[:len(ind_no1), :] for a in Interf])
        self.sub2.setInterfaceData(self.name, interfaceNodes2, [
                                   a[len(ind_no1):, :] for a in Interf])


class SubStructure(WithEmbeddedObjects):
    """ SubStructure object for dynamic substructring.

    The subStructure is fully defined by :
    - its stiffness and mass matrices
    - its modes, computed with blocked interfaces
    """

    # Mandatory class attribute in order to manage aster serialization process
    aster_embedded = ["mass", "mesh", "modes", "nodes", "coords",
                      "dofNumbering", "stiffness", "lINodes", "lIDispl", "iModes", "lIName"]

    def __init__(self, stiffness, mass, modes):

        self.mass = mass
        self.mesh = stiffness.getMesh()
        self.modes = modes
        coords = np.array(self.mesh.getCoordinates().getValues())
        self.nodes = np.array(self.mesh.getNodes())
        coords = np.array(coords)
        Nb_no = len(self.nodes)
        self.coords = coords.reshape((Nb_no, 3))
        self.dofNumbering = stiffness.getDOFNumbering()
        self.lINodes = []
        self.lIDispl = []
        self.iModes = None
        self.lIName = []
        self.stiffness = stiffness

    def __getstate__(self):
        """Method for aster serialization process
        """
        return [self.mass, self.mesh, self.modes, self.nodes, self.coords,
                self.dofNumbering, self.stiffness, self.lINodes, self.lIDispl, self.iModes,
                self.lIName]

    def __setstate__(self, state):
        """Method for aster serialization process
        """
        assert len(state) == 11, state
        (self.mass, self.mesh, self.modes, self.nodes, self.coords, self.dofNumbering,
         self.stiffness, self.lINodes, self.lIDispl, self.iModes, self.lIName) = state

    def setInterfaceData(self, iName: str, iNodes: list, iDispl: list):
        """
        Utility method for the substructure to recieve information
        for its interface
        """
        self.lIName.append(iName)
        self.lINodes.append(iNodes)
        self.lIDispl.append(iDispl)

    def computeInterfaceModes(self, resi_rela=1e-5):
        """
        Compute interface modes - both fluid and structure

        Output :
        --------
        interfModes : interface modes packaged in a code_aster's result
        """

        model = self.stiffness.getModel()
        dofNumbering = self.dofNumbering
        FACTORISER(reuse=self.stiffness, PCENT_PIVOT=150, MATR_ASSE=self.stiffness)

        createOutputFile = True

        for iName, iDispl, iNodes in zip(self.lIName, self.lIDispl, self.lINodes):

            for displ in iDispl:

                Cham = CREA_CHAMP(MODELE=model,
                                  NUME_DDL=dofNumbering,
                                  PROL_ZERO='OUI',
                                  OPERATION='AFFE',
                                  TYPE_CHAM='NOEU_DEPL_R',
                                  AFFE=(_F(TOUT='OUI',
                                           NOM_CMP=('DX',),
                                           VALE=(0.)),
                                        ),)

                for i1, no in enumerate(iNodes):
                    dictDofValues = {}
                    dictDofValues['NOEUD'] = int(no)
                    for j1, comp in enumerate(['DX', 'DY', 'DZ', 'PRES']):
                        dictDofValues[comp] = displ[i1, j1]
                    try:
                        Cham.setDirichletBC(**dictDofValues)
                    except Exception as e:
                        print("Erreur avec noeud {} comp {} et valeur {}".format(no, comp,
                                                                                 displ[i1, j1]))
                        print(e)
                        pass

                resu = RESOUDRE(MATR=self.stiffness, CHAM_NO=Cham, RESI_RELA=-1.)

                if createOutputFile:
                    interfModes = CREA_RESU(OPERATION='AFFE',
                                            TYPE_RESU='MULT_ELAS',
                                            NOM_CHAM='DEPL',
                                            AFFE=_F(CHAM_GD=resu, MODELE=model),
                                            )
                    createOutputFile = False
                else:
                    interfModes = CREA_RESU(reuse=interfModes,
                                            OPERATION='AFFE',
                                            TYPE_RESU='MULT_ELAS',
                                            NOM_CHAM='DEPL',
                                            AFFE=_F(CHAM_GD=resu, MODELE=model),
                                            )
                interfModes.userName = "Resu_{}".format(iName)

            self.iModes = interfModes

        return interfModes

    def exportModesToAsterResult(self, omArray, evArray, precision=None):
        """
        Export eigenvalues and eigenvectors to MODE_MECA

        Input :
        -------
        omArray : array of eigenvalues
        evArray : array of eigenvectors

        Output :
        --------
        MODE_MECA DataStructure
        """

        model = self.stiffness.getModel()
        dofNumbering = self.dofNumbering

        prevEig = -1
        if precision:
            omArrayRound = np.around(omArray.real, precision)
        else:
            omArrayRound = omArray.real
        for idx, eig in enumerate(omArrayRound):
            # pour ne pas avoir deux valeurs propres strictement egal
            if precision:
                if eig == prevEig:
                    eig += 10**-precision
                prevEig = eig
            print('frequence = %s' % eig)
            Cham = CREA_CHAMP(MODELE=model,
                              NUME_DDL=dofNumbering,
                              PROL_ZERO='OUI',
                              OPERATION='AFFE',
                              TYPE_CHAM='NOEU_DEPL_R',
                              AFFE=(_F(TOUT='OUI',
                                        NOM_CMP='DX',
                                        VALE=0.),),)

            Cham.updateValuePointers()
            for i in range(evArray.shape[0]):
                Cham[i] = evArray[i, idx].real  # keep only real part

            if idx == 0:
                Resu = CREA_RESU(OPERATION='AFFE',
                                 TYPE_RESU='MODE_MECA',
                                 NOM_CHAM='DEPL',
                                 AFFE=_F(CHAM_GD=Cham, NUME_MODE=idx, FREQ=eig),
                                 )
            else:
                Resu = CREA_RESU(reuse=Resu,
                                 OPERATION='AFFE',
                                 TYPE_RESU='MODE_MECA',
                                 NOM_CHAM='DEPL',
                                 AFFE=_F(CHAM_GD=Cham, NUME_MODE=idx, FREQ=eig),
                                 )
            Resu.userName = "Modes_{}".format(model.getName())

        return Resu


class Structure(WithEmbeddedObjects):
    """ Structure object for dynamic substructring.

    The subStructure is fully defined by the list of substructures and interfaces
    """

    # Mandatory class attribute in order to manage aster serialization process
    aster_embedded = ["lSubS", "lInterfaces"]

    def __init__(self, lSubS: list, lInterfaces: list):
        self.lSubS = lSubS
        self.lInterfaces = lInterfaces

    def __getstate__(self):
        """Method for aster serialization process
        """
        return [self.lSubS, self.lInterfaces]

    def __setstate__(self, state):
        """Method for aster serialization process
        """
        assert len(state) == 2, state
        (self.lSubS, self.lInterfaces) = state

    def computeGlobalModes(self, nmodes=None, precision=None):
        # donne l'index des modes propres et des modes d'interface
        # de la i-ème sous structure
        lIndexOfInterfModes = []
        lNumberOfInterfModes = []
        lNumberOfPhysicalEqs = []

        # Compute Reduced Matrices
        Kred, Mred, allModes = self._computeReducedOperators(
            lNumberOfPhysicalEqs, lIndexOfInterfModes, lNumberOfInterfModes)
        totalNumberOfModes = allModes.shape[1]

        # Kernel of the constraint matrix
        T = self._computeKernelBasis(
            lNumberOfInterfModes, totalNumberOfModes, lIndexOfInterfModes)

        # Projection on the kernel
        Kredred = T.T.dot(Kred.dot(T))
        Mredred = T.T.dot(Mred.dot(T))

        if False:
            im = plt.imshow(Kredred)
            plt.colorbar(im)
            plt.show()

        # Definition of the global modal problem
        rightVP = NonSymEVP([Kredred, Mredred])

        # Solve
        nmodes = nmodes or rightVP.shape[0] - 2
        omredred, evredred = scipy.sparse.linalg.eigs(rightVP, k=nmodes, which='LR',
                                                      maxiter=1000*nmodes)

        # tri des valeurs propres
        omredred = np.sqrt(1. / omredred) / 2. / np.pi
        tri = np.real(omredred).argsort()
        # print(omredred.shape,evredred.shape)
        omredred = omredred[tri]
        evredred = evredred[:, tri]
        # print(omredred.shape,evredred.shape)

        # Switch back to physical space
        evred = T.dot(evredred)
        ev = allModes.dot(evred)

        decalage = 0
        resu = []
        for isub, sub in enumerate(self.lSubS):
            neq = lNumberOfPhysicalEqs[isub]
            evSub = ev[decalage:decalage + neq, :]
            decalage += neq
            resuSub = sub.exportModesToAsterResult(omredred, evSub, precision)
            resu.append(resuSub)

        return omredred, resu

    def _computeKernelBasis(self, number_interfModes, totalNumberOfModes, index_interfModes):
        # each interface generates n interface modes on the left structure
        # and n interface modes on the right structure
        # there are n constraints to glue it all together
        nConstraints = np.add.reduce(np.array(number_interfModes)) // 2
        iConstraints = 0
        idx1 = np.zeros(nConstraints, dtype=int)
        jdx1 = np.zeros(nConstraints, dtype=int)
        values1 = np.zeros(nConstraints, dtype=float)
        idx2 = np.zeros(nConstraints, dtype=int)
        jdx2 = np.zeros(nConstraints, dtype=int)
        values2 = np.zeros(nConstraints, dtype=float)

        # build the constraint matrix to impose equal components of the interface modes
        for interface in self.lInterfaces:
            sub1 = interface.sub1
            sub2 = interface.sub2
            for isub, sub in enumerate(self.lSubS):
                if sub == sub1:
                    # number of local modes per interface
                    nRanks = len(sub.iModes.getRanks()) // len(sub.lIName)
                    # name of the interface where the interface modes are built
                    inames = [n for n in sub.lIName for _ in range(nRanks)]
                    # keep only the interface modes of the current interface
                    index = np.where(np.array(inames) == interface.name)[0]
                    nConstraintSub1 = len(index)
                    indexConstraintSub = index_interfModes[isub][index]
                    idx1[iConstraints : iConstraints + nConstraintSub1] = np.arange(
                        iConstraints, iConstraints + nConstraintSub1
                    )
                    jdx1[iConstraints : iConstraints + nConstraintSub1] = indexConstraintSub
                    values1[iConstraints : iConstraints + nConstraintSub1] = 1
                if sub == sub2:
                    # number of local modes per interface
                    nRanks = len(sub.iModes.getRanks()) // len(sub.lIName)
                    # name of the interface where the interface modes are built
                    inames = [n for n in sub.lIName for _ in range(nRanks)]
                    index = np.where(np.array(inames) == interface.name)[0]
                    nConstraintSub2 = len(index)
                    indexConstraintSub = index_interfModes[isub][index]
                    idx2[iConstraints : iConstraints + nConstraintSub2] = np.arange(
                        iConstraints, iConstraints + nConstraintSub2
                    )
                    jdx2[iConstraints : iConstraints + nConstraintSub2] = indexConstraintSub
                    values2[iConstraints : iConstraints + nConstraintSub2] = -1
            assert nConstraintSub1 == nConstraintSub2
            iConstraints += nConstraintSub1
        values = np.hstack((values1, values2))
        idx = np.hstack((idx1, idx2))
        jdx = np.hstack((jdx1, jdx2))
        constraintMatrix = scipy.sparse.coo_matrix(
            (values, (idx, jdx)), shape=(nConstraints, totalNumberOfModes)
        )

        # compute the kernel of the constraint matrix
        _, s, v = np.linalg.svd(constraintMatrix.toarray())

        s = np.concatenate((s, np.zeros(constraintMatrix.shape[1] - s.shape[0])))

        # pick almost zero singular values
        inz = np.where(s < 1.0e-15)[0]
        # a basis of ker(constraint matrix)
        T = v[inz, :].transpose()
        # Visualize the constraint matrix
        if False:
            try:
                plt.switch_backend("TkAgg")
                plt.spy(constraintMatrix)
                plt.show()
            except ImportError:
                print("Cannot swith to interactive backend")
        # Verification que l'on a construit une base du noyau
        if False:
            constraintMatrix.dot(T)

        return T

    def _computeReducedOperators(self, lNumberOfPhysicalEqs: list, lIndexOfInterfModes: list,
                                 lNumberOfInterfModes: list):
        Kred = []
        Mred = []
        allModes = []
        decalage = 0

        for sub in self.lSubS:
            stiff = sub.stiffness
            mass = sub.mass
            values, idx, jdx, neq = stiff.EXTR_MATR(sparse=True)
            K = scipy.sparse.coo_matrix((values, (idx, jdx)), shape=(neq, neq))
            values, idx, jdx, neq = mass.EXTR_MATR(sparse=True)
            M = scipy.sparse.coo_matrix((values, (idx, jdx)), shape=(neq, neq))
            lNumberOfPhysicalEqs.append(neq)

            numb = sub.dofNumbering
            # the method is 1-index based
            lag = np.array(numb.getRowsAssociatedToLagrangeMultipliers()) - 1

            modes = sub.modes
            em = [modes.getFieldOnNodesReal('DEPL', r-1).getValues() for r in modes.getRanks()]
            decalage += len(em)
            em = np.vstack(em).T
            em[lag, :] = 0.  # mise a zero des Lagrange

            imodes = sub.iModes
            im = [imodes.getFieldOnNodesReal('DEPL', r-1).getValues() for r in imodes.getRanks()]
            lIndexOfInterfModes.append(np.arange(len(im)) + decalage)
            lNumberOfInterfModes.append(len(im))
            decalage += len(im)
            im = np.vstack(im).T
            im[lag, :] = 0.  # mise a zero des Lagrange

            allModesSub = np.hstack((em, im))
            Kred.append(allModesSub.T.dot(K.dot(allModesSub)))
            Mred.append(allModesSub.T.dot(M.dot(allModesSub)))
            allModes.append(allModesSub)

        allModes = scipy.sparse.block_diag(allModes)
        Kred = scipy.sparse.block_diag(Kred)
        Mred = scipy.sparse.block_diag(Mred)
        return Kred, Mred, allModes


class NonSymEVP(scipy.sparse.linalg.LinearOperator):

    def __init__(self, lMatrices, side='Right'):

        self.K = scipy.sparse.csc_matrix(lMatrices[0])
        self.M = lMatrices[1]
        self.shape = self.K.shape
        self.Kfact = scipy.sparse.linalg.factorized(self.K)

        if side[0].lower() == 'r':
            print('\nFactorizing - Looking for right eigenvectors \n')
        elif side[0].lower() == 'l':
            print('\nFactorizing - Looking for left eigenvectors \n')

    def _matvec(self, x):
        Mx = np.real(self.M.dot(x))
        sol = self.Kfact(Mx)
        return np.real(sol)

    def _rmatvec(self, x):
        Mx = np.real(self.M.dot(x))
        sol = self.Kfact(Mx)
        return np.real(sol)


def macPlot(lres1, lres2, lmass, fluid_material=None, massprod=True, normalize=True,
            name1=None, name2=None, list1=None, list2=None, dof=None, interactive_plot=False,
            save_plot_filename=""):
    """Compute and plot the MAC for 2 data-structures computed for a given modal problem.
    The function can also handle lists of data-structures. This can be handy in the case where
    dynamic substructuring is used. In this case,
    lres1 can be the list of the modal response of each substructure and lres2 can be the global
    response projected on the meshes of the
    aformentioned substructures.
    The fields in each data-structure must be of type displacement (aka DEPL) - displacement and
    pressure components are supported.

    Input :
    res1, res2 (displ data-structure): modal problem results
    mass (matrix) : mass matrix of the modal problem
    massprod (bool) : indicate that the MAC is computed in the norm of the mass matrix
    normalize (bool) : indicate if the MAC must be normalized by its abs max diagonal value
    name1, name2 (str) : names of the data-structures in the figure
    list1, list2 (list[ints]) : indicate which eigenvalues to keep to plot the MAC
    dof (list[str]) : components of the fields kept to computed the MAC
    interactive_plot (bool) : display a plot of the MAC
    save_plot_filename (str) : name of the file where to save the MAC plot
    """

    interactive_is_possible = True
    try:
        plt.switch_backend("TkAgg")  # switch to interactive plot
    except ImportError:
        interactive_is_possible = False
    if HAS_MATPLOTLIB and os.getenv("DISPLAY"):
        plt.figure()
        plt.jet()
    # build lists of modes for each substructure
    if not isinstance(lres1, (list, tuple)):
        lres1 = [lres1]
    if not isinstance(lres2, (list, tuple)):
        lres2 = [lres2]
    if not isinstance(lmass, (list, tuple)):
        lmass = [lmass]
    if len(lres1) != len(lres2) or len(lres1) != len(lmass):
        raise KeyError("res1, res2, and mass must be of same length")
    nStruct = len(lres1)
    # selection of the eigenvalues to be used
    nres1 = np.min(np.array([_res.getNumberOfRanks() for _res in lres1]))
    nres2 = np.min(np.array([_res.getNumberOfRanks() for _res in lres2]))
    if list1 and not all([i in range(nres1) for i in list1]):
        raise KeyError("list1 is out of bound")
    if list2 and not all([i in range(nres2) for i in list2]):
        raise KeyError("list2 is out of bound")
    nModes1 = len(list1) if list1 else nres1
    nModes2 = len(list2) if list2 else nres2
    lMode1 = list1 if list1 else range(nModes1)
    lMode2 = list2 if list2 else range(nModes2)
    lFreq1 = lres1[0].getAccessParameters()["FREQ"]
    lFreq2 = lres2[0].getAccessParameters()["FREQ"]
    rhof = fluid_material.RCVALE("FLUIDE", nomres=("RHO"), stop=2)[0][0] if fluid_material else 1.
    # start MAC computation : MAC = MAC1 **2 / MAC2 / MAC3
    MAC1 = np.zeros((nModes2, nModes1))
    MAC2 = np.zeros((nModes2, nModes1))
    MAC3 = np.zeros((nModes2, nModes1))
    for istru in range(nStruct):
        mass = lmass[istru]
        res1 = lres1[istru]
        res2 = lres2[istru]
        # selection of the dofs
        lPhysical = (
            np.array(mass.getDOFNumbering().getRowsAssociatedToPhysicalDofs()) - 1
        )  # 0-based index
        lDOF = lPhysical
        if dof:
            dict_dof = {}
            for row in lPhysical:
                dd = mass.getDOFNumbering().getComponentAssociatedToRow(
                    int(row) + 1
                )  # 1-based index
                dict_dof.setdefault(dd, []).append(row)
            lDOF = sum([dict_dof[d] for d in dof], [])
        # extract mass matrix in the form a 3 arrays
        Mp = mass.EXTR_MATR(sparse=True)
        # turn it into a sparse matrix
        M = scipy.sparse.coo_matrix((Mp[0], (Mp[1], Mp[2])), shape=(Mp[3], Mp[3]))
        # extract considered dofs (only csr format support it)
        M = M.tocsr()[lDOF, :][:, lDOF]
        # function to retrieve left and right modes from a modal result

        def getLeftAndRightModes(_res, _idx, _imode, _dof, _lFreq):
            vectot = _res.getFieldOnNodesReal("DEPL", _imode).EXTR_COMP().valeurs
            v0_right = (
                np.concatenate(
                    [
                        _res.getFieldOnNodesReal("DEPL", _imode).EXTR_COMP(d).valeurs
                        for d in _dof
                    ]
                )
                if _dof
                else None
            )
            v_right = v0_right if _dof else vectot
            v0_left = (
                np.concatenate(
                    [
                        rhof
                        * (2 * np.pi * _lFreq[_idx]) ** 2
                        * _res.getFieldOnNodesReal("DEPL", _imode).EXTR_COMP(d).valeurs
                        if d in ["DX", "DY", "DZ"]
                        else _res.getFieldOnNodesReal("DEPL", _imode).EXTR_COMP(d).valeurs
                        for d in _dof
                    ]
                )
                if _dof
                else None
            )
            v_left = (
                v0_left
                if _dof
                else np.array(
                    [
                        rhof * (2 * np.pi * _lFreq[_idx]) ** 2 * vectot[id]
                        if d in ["DX", "DY", "DZ"]
                        else vectot[id]
                        for id, d in enumerate(
                            _res.getFieldOnNodesReal("DEPL", _imode).EXTR_COMP(topo=1).comp
                        )
                    ]
                )
            )
            return v_left, v_right

        for idx1, m1 in enumerate(lMode1):
            v1_left, v1_right = getLeftAndRightModes(res1, idx1, m1, dof, lFreq1)
            Mv1 = M.dot(v1_right) if massprod else v1_right
            for idx2, m2 in enumerate(lMode2):
                v2_left, v2_right = getLeftAndRightModes(res2, idx2, m2, dof, lFreq2)
                Mv2 = M.dot(v2_right) if massprod else v2_right
                print(
                    "{:.1f}%% achieved".format(
                        (float(m1 * (nModes2 - 1) + m2) / nModes1 / nModes2) * 100
                    ),
                    end="\r",
                )
                MAC1[idx2, idx1] += np.real(np.dot(v1_left, Mv2))
                MAC2[idx2, idx1] += np.real(np.dot(v1_left, Mv1))
                MAC3[idx2, idx1] += np.real(np.dot(v2_left, Mv2))
    mac = MAC1 ** 2 / MAC2 / MAC3
    print(" " * 100, end="\r")  # in order to clean the progress print
    if normalize:
        mac = np.dot(mac, np.linalg.inv(np.diag(np.amax(mac, axis=0))))
    if HAS_MATPLOTLIB and os.getenv("DISPLAY"):
        plt.imshow(mac, interpolation="nearest")
        plt.grid(False)
        ax = plt.axes()
        label1 = name1 or res1.getName()
        label2 = name2 or res2.getName()
        ax.set_xlabel(label1)
        ax.set_ylabel(label2)
        plt.xticks(range(nModes1), ["{:.1f}".format(lFreq1[idx]) for idx in lMode1], rotation=45)
        plt.yticks(range(nModes2), ["{:.1f}".format(lFreq2[idx]) for idx in lMode2])
        cbar = plt.colorbar()
        cbar.ax.set_ylabel("MAC")
        for (x_val, y_val), val in np.ndenumerate(np.transpose(mac)):
            ax.text(x_val, y_val, "{:.1f}".format(val) if (val > 0.2) else "", va="center",
                    ha="center", color="white",)
        if interactive_plot and interactive_is_possible:
            plt.show()
        if save_plot_filename:
            plt.savefig(save_plot_filename)
    return mac
