! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
subroutine getKinematicForHHO(valeType, model, numeddl, keywordFact, cnsForCharci)
!
    use HHO_precalc_module, only: hhoAddInputField
    use HHO_Dirichlet_module, only: hasHHODoFFromNodes
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterc/getmjm.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/cescel.h"
#include "asterfort/cescns.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/cncinv.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getnode.h"
#include "asterfort/getvr8.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/megeom.h"
!
    character(len=1), intent(in) :: valeType
    character(len=8), intent(in) :: model
    character(len=14), intent(in) :: numeddl
    character(len=16), intent(in) :: keywordFact
    character(len=19), intent(in) :: cnsForCharci
!
! --------------------------------------------------------------------------------------------------
!
! AFFE_CHAR_CINE
!
! Get kinematic values for HHO
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: elnoSVale = "&HHOELNS.VALE"
    character(len=24), parameter :: elnoVale = "&HHOELNO.VALE"
    character(len=24), parameter :: listNode = "&HHO    .LIST_NODE"
    character(len=4), parameter :: prolZero = "OUI"
    character(len=8) :: mesh
    character(len=16) :: phenom
    character(len=24) :: modelLigrel
    integer(kind=8) :: iocc, nocc, nbret, iad, nncp, iret, modelDime, ibid, cmpNume
    integer(kind=8) :: iCell, iNode, iCellNode, iKeywordLieu, iUserDOF
    integer(kind=8) :: userNodeNb, userDOFNb, nbin
    integer(kind=8) :: cellTypeNume, cellDime, cellNume, cellNbNode
    integer(kind=8) :: nodeNume, nodeNbCell, nodeNumeLoca, nbCmpVale
    aster_logical :: l_hho_dofs
    integer(kind=8), pointer :: userNodeNume(:) => null()
    integer(kind=8), pointer :: meshTypmail(:) => null()
    real(kind=8) :: userDOFVale
    character(len=16) :: currentDOF
    integer(kind=8), parameter :: mxcmp = 100
    character(len=16) :: userDOFName(mxcmp)
    character(len=8) :: chcity(mxcmp), physQuanVale
    integer(kind=8), parameter :: npg = -1, nspt = -1
    character(len=8), parameter ::  paraNameVale = "PCMPVALE"
    character(len=16) :: cmpNameVale(3)
    integer(kind=8) :: jvElnoSValeL, jvElnoSValeD, jvElnoSValeV
    integer(kind=8), parameter :: nbKeywordLieu = 5
    character(len=16), parameter :: keywordLieu(nbKeywordLieu) = (/'GROUP_MA', 'MAILLE  ', &
                                                                   'GROUP_NO', 'NOEUD   ', &
                                                                   'TOUT    '/)
    character(len=16) :: option
    integer(kind=8), parameter :: nbxin = 5
    integer(kind=8), parameter :: nbxout = 1
    character(len=8) :: lpaout(nbxout), lpain(nbxin)
    character(len=24) :: lchout(nbxout), lchin(nbxin)
    character(len=24) :: chgeom
    character(len=24), parameter :: elnoHHO = "&HHOELNO.CINE"
    character(len=24), parameter :: elnoHHOS = "&HHOELNS.CINE"
    character(len=19), parameter :: connexInvName = '&&HHOMEC.CONINV'
    integer(kind=8), pointer :: connexInv(:) => null(), connexInvLongCum(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: cataTmTmdim(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

! - Informations sur model
    call dismoi('PHENOMENE', model, 'MODELE', repk=phenom)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('DIM_GEOM', model, 'MODELE', repi=modelDime)

    if (phenom == "MECANIQUE") then
        option = "HHO_CINE_R_MECA"
        physQuanVale = "DEPL_R"
        nbCmpVale = 3
        cmpNameVale(1:3) = (/"DX", "DY", "DZ"/)
    else
        option = "HHO_CINE_R_THER"
        physQuanVale = "TEMP_R"
        nbCmpVale = 1
        cmpNameVale(1) = "TEMP"
    end if

! - Inverse connectivity and mesh paramters
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call cncinv(mesh, [ibid], 0, 'V', connexInvName)
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=meshTypmail)

! - Préparation CHAM_ELEM / ELNO (pour valeurs utilisateurs)
    call cescre("V", elnoSVale, "ELNO", mesh, physQuanVale, &
                nbCmpVale, cmpNameVale, [npg], [nspt], [-nbCmpVale])
    call jeveuo(elnoSVale(1:19)//'.CESD', 'L', jvElnoSValeD)
    call jeveuo(elnoSVale(1:19)//'.CESV', 'E', jvElnoSValeV)
    call jeveuo(elnoSVale(1:19)//'.CESL', 'E', jvElnoSValeL)

! - Lecture des données utilisateurs
    call getfac(keywordFact, nocc)
    do iocc = 1, nocc
! ----- NOEUDS A CONTRAINDRE
        call getnode(mesh, keywordfact, iocc, 'F', listNode, userNodeNb)
        if (userNodeNb .eq. 0) cycle
        call jeveuo(listNode, 'L', vi=userNodeNume)
!
! --------- Test only nodes with HHO dofs
        l_hho_dofs = hasHHODoFFromNodes(numeddl, userNodeNb, listNode)
        if (.not. l_hho_dofs) cycle

        do iNode = 1, userNodeNb
            nodeNume = userNodeNume(iNode)

! --------- Récupération de la maille attachée
            call jeveuo(jexatr(connexInvName, 'LONCUM'), 'L', vi=connexInvLongCum)
            nodeNbCell = connexInvLongCum(nodeNume+1)-connexInvLongCum(nodeNume)
            call jeveuo(jexnum(connexInvName, nodeNume), 'L', vi=connexInv)
!
            ! Optimization - too many cell to be a middle node
            if (nodeNbCell > 3) cycle
!
!           Search element with same dimenison
            nodeNumeLoca = 0
            do iCell = 1, nodeNbCell
! --------- Get element
                cellNume = connexInv(iCell)

! ----------Get type of element and dimension
                cellTypeNume = meshTypmail(cellNume)
                call jeveuo(jexnum('&CATA.TM.TMDIM', cellTypeNume), 'L', vi=cataTmTmdim)
                cellDime = cataTmTmdim(1)

! --------- Get local index of node in element
                if (cellDime .eq. modelDime) then
                    call jeveuo(jexnum(mesh//'.CONNEX', cellNume), 'L', vi=connex)
                    call jelira(jexnum(mesh//'.CONNEX', cellNume), 'LONMAX', cellNbNode)
                    do iCellNode = 1, cellNbNode
                        if (connex(iCellNode) .eq. nodeNume) then
                            nodeNumeLoca = iCellNode
                            exit
                        end if
                    end do

                    exit
                end if
            end do
            ASSERT(nodeNumeLoca .gt. 0)
! --------- DDL A CONTRAINDRE
            call getmjm(keywordFact, iocc, mxcmp, userDOFName, chcity, userDOFNb)
            do iUserDOF = 1, userDOFNb
! ------------- On ne garde que DX, DY, DZ et pas GROUP_MA dans la liste
                currentDOF = userDOFName(iUserDOF)
                do iKeywordLieu = 1, nbKeywordLieu
                    if (currentDOF .eq. keywordLieu(iKeywordLieu)) goto 110
                end do

! ------------- Calcul de l'index de la composante
                if (currentDOF .eq. "DX") then
                    cmpNume = 1
                elseif (currentDOF .eq. "DY") then
                    cmpNume = 2
                elseif (currentDOF .eq. "DZ") then
                    cmpNume = 3
                elseif (currentDOF .eq. "TEMP") then
                    cmpNume = 1
                else
                    ASSERT(ASTER_FALSE)
                end if

! ------------- On récupère la valeur du DDL donné par l'utilisateur
                if (valeType .eq. 'R') then
                    call getvr8(keywordFact, currentDOF, iocc=iocc, &
                                scal=userDOFVale, nbret=nbret)
                else
                    ASSERT(ASTER_FALSE)
                end if

! ------------- On sauve la valeur du DDL dans le CHAM_ELNO_S
                call cesexi('S', jvElnoSValeD, jvElnoSValeL, cellNume, nodeNumeLoca, &
                            1, cmpNume, iad)
                if (abs(iad) > 0) then
                    iad = abs(iad)
                    zl(jvElnoSValeL-1+iad) = ASTER_TRUE
                    zr(jvElnoSValeV-1+iad) = userDOFVale
                end if

110             continue
            end do
        end do
    end do
!
! - Conversion des CHAM_ELNO_S en CHAM_ELNO pour calcul
    call cescel(elnoSVale, modelLigrel, option, paraNameVale, prolZero, &
                nncp, 'V', elnoVale, 'A', iret)
!
    ASSERT(nncp .eq. 0)
!
! - Appel calcul élémentaire (HHO_CINE_MECA)
    call megeom(model, chgeom)
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "
    lpain(1) = "PGEOMER"
    lchin(1) = chgeom
    lpain(2) = paraNameVale
    lchin(2) = elnoVale
    nbin = 2
    call hhoAddInputField(model, nbxin, lchin, lpain, nbin)
!
    lpaout(1) = "PCINE"
    lchout(1) = elnoHHO
    call calcul('S', option, modelLigrel, nbin, lchin, &
                lpain, nbxout, lchout, lpaout, "V", &
                'OUI')
!
!
! - Conversion output in CHAM_NO_S
    call celces(elnoHHO, "V", elnoHHOS)
    call cescns(elnoHHOS, ' ', 'V', cnsForCharci, ' ', ibid)
!
!
! - Cleaning
    call detrsd("CHAM_ELEM_S", elnoSVale)
    call detrsd("CHAM_ELEM", elnoHHO)
    call detrsd("CHAM_ELEM_S", elnoHHOS)
!
end subroutine
