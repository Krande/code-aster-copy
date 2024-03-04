! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
! ==================================================================================================
!
! Module for pairing in contact
!
! ==================================================================================================
!
module ContactPairing_module
! ==================================================================================================
    use contact_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: getContactCellGeom, linearizeCell, calcNormExte, calcNorm
    public :: projPointVector, isPointInCell
    private :: getContactCellType, compLocalBase, getCoorNodePara, getCoorNodeGlob
    private :: projPointVector_
    public :: setProjAlgoPara, fromParaToGLob
! ==================================================================================================
    private
#include "asterc/r8gaem.h"
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/teattr.h"
#include "contact_module.h"
#include "Contact_type.h"
#include "jeveux.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! getContactCellGeom
!
! Get geometry of contact cell
!
! In  nomte            : name of element type
! In  jvGeom           : pointer to geometric field
! IO  contCellGeom     : properties of contact cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getContactCellGeom(nomte, jvGeom, contCellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: nomte
        integer, intent(in) :: jvGeom
        type(Contact_CellGeom), intent(inout) :: contCellGeom
! ----- Local
        integer :: shiftGeom
!   ------------------------------------------------------------------------------------------------
!
        call getContactCellType(nomte, contCellGeom)
        shiftGeom = 0
        call getCoorNodeGlob(jvGeom, contCellGeom%cellDime, contCellGeom%slav, shiftGeom)
        call getCoorNodePara(contCellGeom%slav)
        call getCoorNodeGlob(jvGeom, contCellGeom%cellDime, contCellGeom%mast, shiftGeom)
        call getCoorNodePara(contCellGeom%mast)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getCoorNodeGlob
!
! Get coordinates of cell
!
! In  jvGeom           : pointer to geometric field
! In  modelDime        : global dimension of space (2 or 3)
! IO  cellGeom         : general geometric properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getCoorNodeGlob(jvGeom, modelDime, cellGeom, shiftGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: jvGeom, modelDime
        type(Cell_Geom), intent(inout) :: cellGeom
        integer, intent(inout) :: shiftGeom
! ----- Local
        integer :: iNode, iDime
!   ------------------------------------------------------------------------------------------------
!
        cellGeom%coorNodeGlob = 0.d0
        do iNode = 1, cellGeom%nbNode
            do iDime = 1, modelDime
                cellGeom%coorNodeGlob(iDime, iNode) = zr(jvGeom-1+shiftGeom+iDime)
            end do
            shiftGeom = shiftGeom+modelDime
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getCoorNodePara
!
! Get coordinates of cell in parametric frame
!
! IO  cellGeom         : general geometric properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getCoorNodePara(cellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Cell_Geom), intent(inout) :: cellGeom
!   ------------------------------------------------------------------------------------------------
!
        cellGeom%coorNodePara = 0.d0
        if (cellGeom%cellCode .eq. 'SE3') then
            cellGeom%coorNodePara(1, 1) = -1.d0
            cellGeom%coorNodePara(1, 2) = 1.d0
            cellGeom%coorNodePara(1, 3) = 0.d0
        elseif (cellGeom%cellCode .eq. 'POI1') then
            cellGeom%coorNodePara(1, 1) = 1.d0
        elseif (cellGeom%cellCode .eq. 'XXXX') then
        elseif (cellGeom%cellCode .eq. 'SE2') then
            cellGeom%coorNodePara(1, 1) = -1.d0
            cellGeom%coorNodePara(1, 2) = 1.d0
        elseif (cellGeom%cellCode .eq. 'TR3') then
            cellGeom%coorNodePara(1, 1) = 0.d0
            cellGeom%coorNodePara(2, 1) = 0.d0
            cellGeom%coorNodePara(1, 2) = 1.d0
            cellGeom%coorNodePara(2, 2) = 0.d0
            cellGeom%coorNodePara(1, 3) = 0.d0
            cellGeom%coorNodePara(2, 3) = 1.d0
        elseif (cellGeom%cellCode .eq. 'TR6') then
            cellGeom%coorNodePara(1, 1) = 0.d0
            cellGeom%coorNodePara(2, 1) = 0.d0
            cellGeom%coorNodePara(1, 2) = 1.d0
            cellGeom%coorNodePara(2, 2) = 0.d0
            cellGeom%coorNodePara(1, 3) = 0.d0
            cellGeom%coorNodePara(2, 3) = 1.d0
            cellGeom%coorNodePara(1, 4) = 0.5d0
            cellGeom%coorNodePara(2, 4) = 0.d0
            cellGeom%coorNodePara(1, 5) = 0.5d0
            cellGeom%coorNodePara(2, 5) = 0.5d0
            cellGeom%coorNodePara(1, 6) = 0.d0
            cellGeom%coorNodePara(2, 6) = 0.5d0
        else if (cellGeom%cellCode .eq. 'QU4') then
            cellGeom%coorNodePara(1, 1) = -1.0d0
            cellGeom%coorNodePara(2, 1) = -1.0d0
            cellGeom%coorNodePara(1, 2) = 1.d0
            cellGeom%coorNodePara(2, 2) = -1.d0
            cellGeom%coorNodePara(1, 3) = 1.d0
            cellGeom%coorNodePara(2, 3) = 1.d0
            cellGeom%coorNodePara(1, 4) = -1.d0
            cellGeom%coorNodePara(2, 4) = 1.d0
        else if (cellGeom%cellCode .eq. 'QU8') then
            cellGeom%coorNodePara(1, 1) = -1.0d0
            cellGeom%coorNodePara(2, 1) = -1.0d0
            cellGeom%coorNodePara(1, 2) = 1.d0
            cellGeom%coorNodePara(2, 2) = -1.d0
            cellGeom%coorNodePara(1, 3) = 1.d0
            cellGeom%coorNodePara(2, 3) = 1.d0
            cellGeom%coorNodePara(1, 4) = -1.d0
            cellGeom%coorNodePara(2, 4) = 1.d0
            cellGeom%coorNodePara(1, 5) = 0.d0
            cellGeom%coorNodePara(2, 5) = -1.d0
            cellGeom%coorNodePara(1, 6) = 1.d0
            cellGeom%coorNodePara(2, 6) = 0.d0
            cellGeom%coorNodePara(1, 7) = 0.d0
            cellGeom%coorNodePara(2, 7) = 1.d0
            cellGeom%coorNodePara(1, 8) = -1.d0
            cellGeom%coorNodePara(2, 8) = 0.d0
        else if (cellGeom%cellCode .eq. 'QU9') then
            cellGeom%coorNodePara(1, 1) = -1.0d0
            cellGeom%coorNodePara(2, 1) = -1.0d0
            cellGeom%coorNodePara(1, 2) = 1.d0
            cellGeom%coorNodePara(2, 2) = -1.d0
            cellGeom%coorNodePara(1, 3) = 1.d0
            cellGeom%coorNodePara(2, 3) = 1.d0
            cellGeom%coorNodePara(1, 4) = -1.d0
            cellGeom%coorNodePara(2, 4) = 1.d0
            cellGeom%coorNodePara(1, 5) = 0.d0
            cellGeom%coorNodePara(2, 5) = -1.d0
            cellGeom%coorNodePara(1, 6) = 1.d0
            cellGeom%coorNodePara(2, 6) = 0.d0
            cellGeom%coorNodePara(1, 7) = 0.d0
            cellGeom%coorNodePara(2, 7) = 1.d0
            cellGeom%coorNodePara(1, 8) = -1.d0
            cellGeom%coorNodePara(2, 8) = 0.d0
            cellGeom%coorNodePara(1, 9) = 0.d0
            cellGeom%coorNodePara(2, 9) = 0.d0
        else
            WRITE (6, *) "cellGeom%cellCode: ", cellGeom%cellCode
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getContactCellType
!
! Get type of contact cell
!
! In  nomte            : name of element type
! IO  cellGeom         : general geometric properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getContactCellType(nomte, cellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: nomte
        type(Contact_CellGeom), intent(inout) :: cellGeom
! ----- Local
        integer :: cellDime, nbNode, ier
        character(len=8) :: cellCode, attrib
!   ------------------------------------------------------------------------------------------------
!
        nbNode = 0
        cellDime = 0

! ----- Slave side
        select case (nomte(3:4))
        case ("S2")
            cellDime = 1
            cellCode = 'SE2'
            nbNode = 2
        case ("S3")
            cellDime = 1
            cellCode = 'SE3'
            nbNode = 3
        case ("Q4")
            cellDime = 2
            cellCode = 'QU4'
            nbNode = 4
        case ("Q8")
            cellDime = 2
            cellCode = 'QU8'
            nbNode = 8
        case ("Q9")
            cellDime = 2
            cellCode = 'QU9'
            nbNode = 9
        case ("T3")
            cellDime = 2
            cellCode = 'TR3'
            nbNode = 3
        case ("T6")
            cellDime = 2
            cellCode = 'TR6'
            nbNode = 6
        case ("P1")
            cellCode = 'POI1'
            nbNode = 1
            if (nomte(1:6) .eq. 'CMP1L2' .or. nomte(1:6) .eq. 'CMP1N2') then
                cellDime = 2
            elseif (nomte(1:6) .eq. 'CMP1L3' .or. nomte(1:6) .eq. 'CMP1N3') then
                cellDime = 3
            else
                cellDime = 0
            end if
        case default
            ASSERT(ASTER_FALSE)
        end select
        cellGeom%slav%cellDime = cellDime
        cellGeom%slav%nbNode = nbNode
        cellGeom%slav%cellCode = cellCode
        ASSERT(nbNode .le. 9)

! ----- Master side
        select case (nomte(5:6))
        case ("S2")
            cellDime = 1
            cellCode = 'SE2'
            nbNode = 2
        case ("S3")
            cellDime = 1
            cellCode = 'SE3'
            nbNode = 3
        case ("Q4")
            cellDime = 2
            cellCode = 'QU4'
            nbNode = 4
        case ("Q8")
            cellDime = 2
            cellCode = 'QU8'
            nbNode = 8
        case ("Q9")
            cellDime = 2
            cellCode = 'QU9'
            nbNode = 9
        case ("T3")
            cellDime = 3
            cellCode = 'TR3'
            nbNode = 3
        case ("T6")
            cellDime = 2
            cellCode = 'TR6'
            nbNode = 6
        case ("L2")
            cellCode = 'XXXX'
            nbNode = 0
            cellDime = 0
        case ("L3")
            cellCode = 'XXXX'
            nbNode = 0
            cellDime = 0
        case ("N2")
            cellCode = 'XXXX'
            nbNode = 0
            cellDime = 0
        case ("N3")
            cellCode = 'XXXX'
            nbNode = 0
            cellDime = 0
        case default
            ASSERT(ASTER_FALSE)
        end select
        ASSERT(nbNode .le. 9)
        cellGeom%mast%cellDime = cellDime
        cellGeom%mast%nbNode = nbNode
        cellGeom%mast%cellCode = cellCode
        call teattr('S', 'DIM_COOR_MODELI', attrib, ier)
        read (attrib, '(I8)') cellGeom%cellDime
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! linearizeCell
!
! Linearization of cell
!
! In  cellGeom         : geometry of cell
! Out cellGeomLine     : geometry of cell after linearization
!
! --------------------------------------------------------------------------------------------------
    subroutine linearizeCell(cellGeom, cellGeomLine)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Cell_Geom), intent(in) :: cellGeom
        type(Cell_Geom), intent(out) :: cellGeomLine
!   ------------------------------------------------------------------------------------------------
!
        cellGeomLine%cellDime = cellGeom%cellDime
        cellGeomLine%coorNodeGlob = cellGeom%coorNodeGlob
        cellGeomLine%coorNodePara = cellGeom%coorNodePara
        if (cellGeom%cellCode .eq. 'SE2' .or. cellGeom%cellCode .eq. 'SE3') then
            cellGeomLine%cellCode = 'SE2'
            cellGeomLine%nbNode = 2
        elseif (cellGeom%cellCode .eq. 'TR3' .or. cellGeom%cellCode .eq. 'TR6') then
            cellGeomLine%cellCode = 'TR3'
            cellGeomLine%nbNode = 3
        else if (cellGeom%cellCode .eq. 'QU4' .or. &
                 cellGeom%cellCode .eq. 'QU8' .or. &
                 cellGeom%cellCode .eq. 'QU9') then
            cellGeomLine%cellCode = 'QU4'
            cellGeomLine%nbNode = 4
        else
            ASSERT(ASTER_FALSE)
        end if
        cellGeomLine%coorNodeGlob(:, cellGeomLine%nbNode+1:9) = 0.d0
        cellGeomLine%coorNodePara(:, cellGeomLine%nbNode+1:9) = 0.d0
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! projPointVector
!
! Projection of a point on surface by given vector
!
! In  projAlgoPara     : parameters of projection algorithm
! In  modelDime        : global dimension of space (2 or 3)
! In  cellTarget       : target cell
! IO  projPara         : parameters of projection
!
! --------------------------------------------------------------------------------------------------
    subroutine projPointVector_(projAlgoPara, modelDime, cellTarget, projPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Contact_ProjAlgoPara), intent(in) :: projAlgoPara
        integer, intent(in) :: modelDime
        type(Cell_Geom), intent(in) :: cellTarget
        type(Contact_ProjPara), intent(inout) :: projPara
! ----- Local
        real(kind=8), parameter :: zero = 0.d0, one = 1.d0
        integer :: iNode, iDime
        real(kind=8) :: ff(9), dff(3, 9)
        real(kind=8) :: vect_posi(3), dist
        real(kind=8) :: residu(3), matrix(3, 3), det
        real(kind=8) :: dksi(2), dbeta, beta
        integer :: iterNewt
        real(kind=8) :: toleAbso, toleRela, toleNewt
        real(kind=8) :: distMini, ksiMini(2), betaMini
        real(kind=8) :: refe, test
!   ------------------------------------------------------------------------------------------------
!
        projPara%errorCode = 0
        projPara%ksi = zero
        beta = one
        iterNewt = 0
        toleAbso = projAlgoPara%newtTole/100.d0
        toleRela = projAlgoPara%newtTole
        distMini = r8gaem()

! ----- Newton loop
20      continue
        vect_posi = zero
        projPara%tau1 = zero
        projPara%tau2 = zero
        matrix = zero
        residu = zero
        dksi = zero
        dbeta = zero

! ----- Get shape functions and derivated shape functions
        call elrfvf(cellTarget%cellCode, projPara%ksi, ff)
        call elrfdf(cellTarget%cellCode, projPara%ksi, dff)

! ----- Position vector of current point
        do iDime = 1, 3
            do iNode = 1, cellTarget%nbNode
                vect_posi(iDime) = &
                    cellTarget%coorNodeGlob(iDime, iNode)*ff(iNode)+vect_posi(iDime)
            end do
        end do

! ----- Compute local base
        call compLocalBase(cellTarget, modelDime, dff, projPara%tau1, projPara%tau2)

! ----- Quantity to minimize
        do iDime = 1, 3
            vect_posi(iDime) = projPara%pointCoor(iDime)-vect_posi(iDime)
        end do
        dist = sqrt(vect_posi(1)*vect_posi(1)+vect_posi(2)*vect_posi(2)+vect_posi(3)*vect_posi(3))

! ----- Newton residual
        do iDime = 1, 3
            residu(iDime) = vect_posi(iDime)-beta*projPara%projVect(iDime)
        end do

! ----- Tangent matrix (Newton)
        do iDime = 1, 3
            matrix(iDime, 1) = projPara%tau1(iDime)
            if (modelDime .eq. 2) then
                matrix(iDime, 2) = projPara%projVect(iDime)
            elseif (modelDime .eq. 3) then
                matrix(iDime, 2) = projPara%tau2(iDime)
                matrix(iDime, 3) = projPara%projVect(iDime)
            else
                ASSERT(ASTER_FALSE)
            end if
        end do

! ----- System determinant
        if (modelDime .eq. 2) then
            det = matrix(1, 1)*matrix(2, 2)-matrix(1, 2)*matrix(2, 1)
        else if (modelDime .eq. 3) then
            det = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3)-matrix(3, 2)*matrix(2, 3))- &
                  matrix(2, 1)*(matrix(1, 2)*matrix(3, 3)-matrix(3, 2)*matrix(1, 3))+ &
                  matrix(3, 1)*(matrix(1, 2)*matrix(2, 3)-matrix(2, 2)*matrix(1, 3))
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (abs(det) .le. r8prem()) then
            projPara%errorCode = 1
            goto 999
        end if

! ----- Solve system
        if (modelDime .eq. 2) then
            dksi(1) = (residu(1)*matrix(2, 2)-residu(2)*matrix(1, 2))/det
            dksi(2) = 0.d0
            dbeta = (residu(2)*matrix(1, 1)-residu(1)*matrix(2, 1))/det
        else if (modelDime .eq. 3) then
            dksi(1) = (residu(1)*(matrix(2, 2)*matrix(3, 3)-matrix(3, 2)*matrix(2, 3))+ &
                       residu(2)*(matrix(3, 2)*matrix(1, 3)-matrix(1, 2)*matrix(3, 3))+ &
                       residu(3)*(matrix(1, 2)*matrix(2, 3)-matrix(2, 2)*matrix(1, 3)))/det
            dksi(2) = (residu(1)*(matrix(3, 1)*matrix(2, 3)-matrix(2, 1)*matrix(3, 3))+ &
                       residu(2)*(matrix(1, 1)*matrix(3, 3)-matrix(3, 1)*matrix(1, 3))+ &
                       residu(3)*(matrix(2, 1)*matrix(1, 3)-matrix(2, 3)*matrix(1, 1)))/det
            dbeta = (residu(1)*(matrix(2, 1)*matrix(3, 2)-matrix(3, 1)*matrix(2, 2))+ &
                     residu(2)*(matrix(3, 1)*matrix(1, 2)-matrix(1, 1)*matrix(3, 2))+ &
                     residu(3)*(matrix(1, 1)*matrix(2, 2)-matrix(2, 1)*matrix(1, 2)))/det
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Update
        projPara%ksi = projPara%ksi+dksi
        beta = beta+dbeta

! ----- Save values if Newton avoids
        if (dist .le. distMini) then
            distMini = dist
            ksiMini = projPara%ksi
            betaMini = beta
        end if

! ----- Convergence
        refe = (projPara%ksi(1)*projPara%ksi(1)+projPara%ksi(2)*projPara%ksi(2)+beta*beta)
        if (refe .le. toleRela) then
            toleNewt = toleAbso
            test = sqrt(dksi(1)*dksi(1)+dksi(2)*dksi(2)+dbeta*dbeta)
        else
            toleNewt = toleRela
            test = sqrt(dksi(1)*dksi(1)+dksi(2)*dksi(2)+dbeta*dbeta)/sqrt(refe)
        end if

! ----- Continue or not ?
        if ((test .gt. toleNewt) .and. (iterNewt .lt. projAlgoPara%newtIterMaxi)) then
            iterNewt = iterNewt+1
            goto 20
        else if ((iterNewt .ge. projAlgoPara%newtIterMaxi) .and. (test .gt. toleNewt)) then
            projPara%ksi = ksiMini
            call elrfvf(cellTarget%cellCode, projPara%ksi, ff)
            call elrfdf(cellTarget%cellCode, projPara%ksi, dff)
            call compLocalBase(cellTarget, modelDime, dff, projPara%tau1, projPara%tau2)
            projPara%errorCode = 1
        end if

! ----- End of loop
999     continue
!
        if (projAlgoPara%newtDebug) then
            write (6, *) "Dimension de l'espace: ", modelDime
            write (6, *) "Coordonnées du point à projeter: ", &
                projPara%pointCoor(1), &
                projPara%pointCoor(2), &
                projPara%pointCoor(3)
            write (6, *) "Direction de projection: ", &
                projPara%projVect(1), &
                projPara%projVect(2), &
                projPara%projVect(3)
            write (6, *) 'Maille cible de la projection: ', &
                cellTarget%cellCode, &
                cellTarget%nbNode
            do iNode = 1, cellTarget%nbNode
                write (6, *) '  Noeud ', iNode
                write (6, *) '   (X,Y,Z)', &
                    cellTarget%coorNodeGlob(1, iNode), &
                    cellTarget%coorNodeGlob(2, iNode), &
                    cellTarget%coorNodeGlob(3, iNode)
            end do
            write (6, *) 'KSI   : ', projPara%ksi(1), projPara%ksi(2)
            write (6, *) 'BETA  : ', beta
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLocalBase
!
! Compute local base at current point
!
! In  cellGeom         : geometry of cell
! In  modelDime        : global dimension of space (2 or 3)
! In  dff              : values of derivative of shape functions at current point
! Out tau1             : first tangent at current point
! Out tau2             : second tangent at current point
!
! --------------------------------------------------------------------------------------------------
    subroutine compLocalBase(cellGeom, modelDime, dff, tau1, tau2)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Cell_Geom), intent(in) :: cellGeom
        integer, intent(in) :: modelDime
        real(kind=8), intent(in) :: dff(3, 9)
        real(kind=8), intent(out) :: tau1(3), tau2(3)
! ----- Local
        real(kind=8), parameter :: zero = 0.d0
        integer :: iNode, iDime
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(cellGeom%nbNode .le. 9)
        tau1 = zero
        tau2 = zero
        do iDime = 1, 3
            do iNode = 1, cellGeom%nbNode
                tau1(iDime) = cellGeom%coorNodeGlob(iDime, iNode)*dff(1, iNode)+tau1(iDime)
                if (modelDime .eq. 3) then
                    tau2(iDime) = cellGeom%coorNodeGlob(iDime, iNode)*dff(2, iNode)+tau2(iDime)
                end if
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setProjAlgoPara
!
! Set parameters of projection algorithm
!
! In  projPhase        : phase of projection
! In  debug            : flag for debug
! Out projAlgoPara     : parameters of projection algorithm
!
! --------------------------------------------------------------------------------------------------
    subroutine setProjAlgoPara(projPhase, debug, projAlgoPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: projPhase
        aster_logical, intent(in) :: debug
        type(Contact_ProjAlgoPara), intent(out) :: projAlgoPara
!   ------------------------------------------------------------------------------------------------
!
        if (projPhase .eq. PROJ_ALGO_GAPI) then
            projAlgoPara%newtDebug = debug
            projAlgoPara%newtIterMaxi = 75
            projAlgoPara%newtTole = PROJ_TOLE
            projAlgoPara%withPrepLine = ASTER_TRUE
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! calcNormExte
!
! Compute external norm to cell
!
! In  cellGeom         : geometry of cell
! In  modelDime        : global dimension of space (2 or 3)
! In  indxNode         : index of current node in cell
! Out normExte         : external norm to cell
!
! --------------------------------------------------------------------------------------------------
    subroutine calcNormExte(cellGeom, modelDime, indxNode, normExte)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Cell_Geom), intent(in) :: cellGeom
        integer, intent(in) :: modelDime, indxNode
        real(kind=8), intent(out) :: normExte(3)
! ----- Local
        real(kind=8) :: dff(3, 9), tau1(3), tau2(3), norm(3), coorPara(2)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNode .le. cellGeom%nbNode)
        coorPara = 0.d0
        coorPara(1) = cellGeom%coorNodePara(1, indxNode)
        if (modelDime .eq. 3) then
            coorPara(2) = cellGeom%coorNodePara(2, indxNode)
        end if
        normExte = 0.d0
        call elrfdf(cellGeom%cellCode, coorPara, dff)
        call compLocalBase(cellGeom, modelDime, dff, tau1, tau2)
        call calcNorm(modelDime, tau1, tau2, norm)
        normExte = -norm
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! calcNorm
!
! Compute norm
!
! In  modelDime        : global dimension of space (2 or 3)
! In  tau1             : first tangent
! In  tau2             : second tangent
! Out norm             : norm
!
! --------------------------------------------------------------------------------------------------
    subroutine calcNorm(modelDime, tau1, tau2, norm)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: modelDime
        real(kind=8), intent(in) :: tau1(3), tau2(3)
        real(kind=8), intent(out) :: norm(3)
! ----- Local
        real(kind=8) :: noor
!   ------------------------------------------------------------------------------------------------
!
        norm = 0.d0
        noor = r8prem()
        if (modelDime .eq. 2) then
            norm(1) = -tau1(2)
            norm(2) = tau1(1)
            norm(3) = 0.d0
        else if (modelDime .eq. 3) then
            call provec(tau2, tau1, norm)
        else
            ASSERT(ASTER_FALSE)
        end if
        call normev(norm, noor)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isPointInCell
!
! Is point in cell ?
!
! In  cellCode         : type of cell
! In  toleInside       : tolerance to detect
! In  ksi              : parametric coordinates of point
! Out pointIsInside
!
! --------------------------------------------------------------------------------------------------
    subroutine isPointInCell(cellCode, ksi, toleInside, pointIsInside)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: cellCode
        real(kind=8), intent(in) :: ksi(2), toleInside
        aster_logical, intent(out) :: pointIsInside
!   ------------------------------------------------------------------------------------------------
!
        pointIsInside = ASTER_FALSE
        if (cellCode .eq. 'SE2' .or. cellCode .eq. 'SE3') then
            if (ksi(1) .ge. (-1.d0-toleInside) .and. &
                ksi(1) .le. (1.d0+toleInside)) then
                pointIsInside = ASTER_TRUE
            end if
        elseif (cellCode .eq. 'TR3' .or. cellCode .eq. 'TR6') then
            if (ksi(1) .ge. -toleInside .and. &
                ksi(2) .ge. -toleInside .and. &
                (ksi(2)+ksi(1)) .le. (1.d0+toleInside)) then
                pointIsInside = ASTER_TRUE
            end if
        elseif (cellCode .eq. 'QU4' .or. cellCode .eq. 'QU8' .or. cellCode .eq. 'QU9') then
            if (ksi(1) .ge. (-1.d0-toleInside) .and. &
                ksi(1) .le. (1.d0+toleInside) .and. &
                ksi(2) .ge. (-1.d0-toleInside) .and. &
                ksi(2) .le. (1.d0+toleInside)) then
                pointIsInside = ASTER_TRUE
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! projPointVector
!
! Projection of a point on surface by given vector
!
! In  projAlgoPara     : parameters of projection algorithm
! IO  projPara         : parameters of projection
!
! --------------------------------------------------------------------------------------------------
    subroutine projPointVector(projAlgoPara, projPara, pointIsInside)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Contact_ProjAlgoPara), intent(in) :: projAlgoPara
        type(Contact_ProjPara), intent(inout) :: projPara
        aster_logical, intent(out) :: pointIsInside
! ----- Local
        type(Cell_Geom) :: cellTarget
        !aster_logical :: pointIsInside
        integer :: modelDime
!   ------------------------------------------------------------------------------------------------
!
        modelDime = projPara%modelDime
        pointIsInside = ASTER_FALSE
        if (projAlgoPara%withPrepLine) then
            cellTarget = projPara%geomTargetLine
            call projPointVector_(projAlgoPara, modelDime, cellTarget, projPara)
            ASSERT(projPara%errorCode == 0)
            call isPointInCell(cellTarget%cellCode, projPara%ksi, &
                               PROJ_TOLE, pointIsInside)
        end if
        if (pointIsInside .or. (.not. projAlgoPara%withPrepLine)) then
            cellTarget = projPara%geomTarget
            call projPointVector_(projAlgoPara, modelDime, cellTarget, projPara)
            ASSERT(projPara%errorCode == 0)
            call isPointInCell(cellTarget%cellCode, projPara%ksi, &
                               PROJ_TOLE, pointIsInside)
            ASSERT(pointIsInside)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! fromParaToGLob
!
! Projection of point from parametric space to global space
!
! In  cellGeom         : geometry of cell
! In  modelDime        : global dimension of space (2 or 3)
! In  coorNodePara     : coordinates of point in parametric space
! Out coorNodeGlob     : coordinates of point in global space
!
! --------------------------------------------------------------------------------------------------
    subroutine fromParaToGLob(cellGeom, modelDime, coorNodePara, coorNodeGlob)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: modelDime
        type(Cell_Geom), intent(in) :: cellGeom
        real(kind=8), intent(in) :: coorNodePara(2)
        real(kind=8), intent(out) :: coorNodeGlob(3)
! ----- Local
        real(kind=8) :: ff(9)
        integer :: iNode, iDime
!   ------------------------------------------------------------------------------------------------
!
        coorNodeGlob = 0.d0
        ff = 0.d0
        call elrfvf(cellGeom%cellCode, coorNodePara, ff)
        do iDime = 1, modelDime
            do iNode = 1, cellGeom%nbNode
                coorNodeGlob(iDime) = coorNodeGlob(iDime)+ &
                                      cellGeom%coorNodeGlob(iDime, iNode)*ff(iNode)
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module
