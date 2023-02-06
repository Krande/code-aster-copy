! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
! Module for damping
!
! ==================================================================================================
!
module Damping_module
! ==================================================================================================
    use Damping_type
    use NonLinearElem_module, only: elemDamp, asseDamp
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: dampModalGetParameters, dampModalPreparation
    public :: dampComputeMatrix, dampComputeVector
! ==================================================================================================
    private
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8pi.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mginfo.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "blas/dcopy.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! dampModalGetParameters
!
! Get parameters for modal damping
!
! In  factorKeyword    : factor keyword
! In  jvDataDamp       : name of datastructure to save parameters
! Out modalDamping     : parameters for modal damping
!
! --------------------------------------------------------------------------------------------------
    subroutine dampModalGetParameters(factorKeyword, jvDataDamp, modalDamping)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=16), intent(in) :: factorKeyword
        character(len=24), intent(in) :: jvDataDamp
        type(MODAL_DAMPING), intent(out) :: modalDamping
! - Local
        character(len=24), parameter :: jvListDamp = "&&NMMOAM.AMORTISSEMENT"
        character(len=16) :: answer
        character(len=8) :: listDamp
        integer :: iret, iMode
        integer :: nbMode, nbModeDS, nbModeMax, nbDampVale
        aster_logical :: lReacVite
        character(len=8) :: dampMode
        aster_logical :: lReducedDampPy, lReducedDampList
        real(kind=8), pointer :: dampVale(:) => null()
        real(kind=8), pointer :: dampValeTemp(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        lReacVite = ASTER_FALSE
        call getvtx(factorKeyword, 'REAC_VITE', iocc=1, scal=answer, nbret=iret)
        if (answer .eq. 'OUI') then
            lReacVite = ASTER_TRUE
        end if

! - Get modes
        call getvid(factorKeyword, 'MODE_MECA', iocc=1, scal=dampMode, nbret=iret)
        if (iret .eq. 0) then
            call utmess('F', 'DAMPING1_20')
        end if
        call jelira(dampMode(1:8)//'           .ORDR', 'LONMAX', nbModeDS)
        nbMode = nbModeDS
        call getvis(factorKeyword, 'NB_MODE', iocc=1, scal=nbModeMax, nbret=iret)
        if (iret .ne. 0) then
            nbMode = nbModeMax
        end if
        if (nbModeMax .ne. nbModeDS) then
            nbMode = min(nbModeDS, nbModeMax)
            call utmess('I', 'DAMPING1_30', ni=3, vali=[nbModeDS, nbModeMax, nbMode])
        end if

! - Get list of reduced damping values: by vector from Python or by list_r8 datastructure
        call getvr8(factorKeyword, 'AMOR_REDUIT', iocc=1, nbval=0, nbret=nbDampVale)
        nbDampVale = -nbDampVale
        lReducedDampPy = nbDampVale .ne. 0
        call getvid(factorKeyword, 'LIST_AMOR', iocc=1, nbval=0, nbret=iret)
        lReducedDampList = iret .ne. 0
        if (.not. lReducedDampPy .and. .not. lReducedDampList) then
            call utmess('F', 'DAMPING1_21')
        end if
        ASSERT(.not. (lReducedDampPy .and. lReducedDampList))

        if (lReducedDampPy) then
            call wkvect(jvListDamp, 'V V R', nbDampVale, vr=dampVale)
            call getvr8(factorKeyword, 'AMOR_REDUIT', iocc=1, nbval=nbDampVale, &
                        vect=dampVale, nbret=iret)
        else
            call getvid(factorKeyword, 'LIST_AMOR', iocc=1, scal=listDamp, nbret=iret)
            call jelira(listDamp//'           .VALE', 'LONMAX', ival=nbDampVale)
            call jeveuo(listDamp//'           .VALE', 'L', vr=dampValeTemp)
            call wkvect(jvListDamp, 'V V R', nbDampVale, vr=dampVale)
            dampVale(1:nbDampVale) = dampValeTemp(1:nbDampVale)
        end if

        if (nbDampVale .gt. nbMode) then
            call utmess('A', 'DAMPING1_19')
        end if
        if (nbDampVale .lt. nbMode) then
            AS_ALLOCATE(vr=dampValeTemp, size=nbMode)
            dampValeTemp(1:nbMode) = dampVale(1:nbMode)
            do iMode = nbDampVale+1, nbMode
                dampValeTemp(iMode) = dampVale(nbDampVale)
            end do
            call jedetr(jvListDamp)
            nbDampVale = nbMode
            call wkvect(jvListDamp, 'V V R', nbDampVale, vr=dampVale)
            dampVale(1:nbDampVale) = dampValeTemp(1:nbDampVale)
            AS_DEALLOCATE(vr=dampValeTemp)
        end if

! - Save parameters
        modalDamping%lReacVite = lReacVite
        modalDamping%dampMode = dampMode
        modalDamping%nbMode = nbMode
        modalDamping%jvListDamp = jvListDamp
        modalDamping%nbDampVale = nbDampVale
        modalDamping%jvDataDamp = jvDataDamp

! - Debug
        if (modalDamping%debug) then
            WRITE (6, *) "Présence d'amortissement modal"
            WRITE (6, *) " Modes mécaniques: ", dampMode
            WRITE (6, *) " Nombre de modes retenus: ", nbMode
            WRITE (6, *) " Nombre d'amortissement modaux: ", nbDampVale
            call jeveuo(jvListDamp, 'L', vr=dampVale)
            WRITE (6, *) "    => ", dampVale(1:nbDampVale)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dampModalPreparation
!
! Preparation for modal damping
!
! In  modalDamping     : parameters for modal damping
!
! --------------------------------------------------------------------------------------------------
    subroutine dampModalPreparation(modalDamping)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(MODAL_DAMPING), intent(in) :: modalDamping
! - Local
        character(len=24) :: matrRigi, modeVect
        character(len=24) :: jvListDamp, jvDataDamp
        character(len=14) :: numeDof
        character(len=8) :: dampMode
        integer :: nbMode, nbEqua
        integer :: iret, lmat, jvDeeq, iMode
        integer :: jvPara, jvDataDampBase
        real(kind=8), pointer :: dampVale(:) => null()
        real(kind=8), pointer :: valeTemp(:) => null()
        real(kind=8), pointer :: modeVale(:) => null()
        real(kind=8), pointer :: dataDampVale(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        dampMode = modalDamping%dampMode
        nbMode = modalDamping%nbMode
        jvListDamp = modalDamping%jvListDamp
        jvDataDamp = modalDamping%jvDataDamp

! - List of damping values
        call jeveuo(jvListDamp, 'L', vr=dampVale)

! - Get parameters from modes
        call mginfo(dampMode, numeDof_=numeDof, nbEqua_=nbEqua)
        call jeveuo(numeDof(1:14)//'.NUME.DEEQ', 'L', jvDeeq)

! - Get rigidity matrix
        call dismoi('REF_RIGI_PREM', dampMode, 'RESU_DYNA', repk=matrRigi)
        call mtdscr(matrRigi(1:8))
        call jeveuo(matrRigi(1:19)//'.&INT', 'E', lmat)

! - Create modal values
        call wkvect(jvDataDamp(1:19)//'.VALM', 'V V R', 3*nbMode, vr=dataDampVale)
        do iMode = 1, nbMode
            call rsadpa(dampMode, 'L', 1, 'MASS_GENE', iMode, 0, sjv=jvPara)
            dataDampVale(3*(iMode-1)+1) = zr(jvPara)
            call rsadpa(dampMode, 'L', 1, 'FREQ', iMode, 0, sjv=jvPara)
            dataDampVale(3*(iMode-1)+2) = zr(jvPara)*2.d0*r8pi()
            dataDampVale(3*(iMode-1)+3) = dampVale(iMode)
        end do

! - Create base
        call wkvect(jvDataDamp(1:19)//'.BASM', 'V V R', nbMode*nbEqua, jvDataDampBase)
        AS_ALLOCATE(vr=valeTemp, size=nbEqua)
        do iMode = 1, nbMode
            call rsexch('F', dampMode, 'DEPL', iMode, modeVect, iret)
            call jeveuo(modeVect(1:19)//'.VALE', 'L', vr=modeVale)
            call dcopy(nbEqua, modeVale, 1, valeTemp, 1)
            call zerlag(nbEqua, zi(jvDeeq), vectr=valeTemp)
            call mrmult('ZERO', lmat, valeTemp, zr(jvDataDampBase+(iMode-1)*nbEqua), 1, ASTER_TRUE)
            call zerlag(nbEqua, zi(jvDeeq), vectr=zr(jvDataDampBase+(iMode-1)*nbEqua))
        end do
        AS_DEALLOCATE(vr=valeTemp)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dampComputeMatrix
!
! Compute damping matrix
!
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  materialField    : field for material parameters
! In  materialCoding   : field for coding material parameters
! In  behaviourField   : field for behaviour parameters
! In  vari             : internal state variables
! In  time             : value of time
! In  listLoad         : name of datastructure for list of loads
! In  numeDof          : name of numbering (NUME_DDL)
! In  rigiElem         : name of elementary matrices for rigidity
! In  massElem         : name of elementary matrices for mass
! In  dampElem         : name of elementary matrices for damping
! In  dampAsse         : name of assembled matrice for damp
!
! --------------------------------------------------------------------------------------------------
    subroutine dampComputeMatrix(model, caraElem, &
                                 materialField, materialCoding, &
                                 behaviourField, &
                                 vari, time, listLoad, numeDof, &
                                 rigiElem, massElem, dampElem, &
                                 dampAsse)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=24), intent(in) :: model, caraElem
        character(len=24), intent(in) :: materialField, materialCoding
        character(len=24), intent(in) :: behaviourField
        character(len=24), intent(in) :: vari
        real(kind=8), intent(in) :: time
        character(len=14), intent(in) :: numeDof
        character(len=19), intent(in) :: listLoad
        character(len=24), intent(in) :: dampElem, rigiElem, massElem, dampAsse
!   ------------------------------------------------------------------------------------------------
!

! - Compute elementary matrices for damp
        call elemDamp(model, caraElem, &
                      materialField, materialCoding, &
                      behaviourField, &
                      vari, time, &
                      rigiElem, massElem, &
                      dampElem)

! - Assemble elementary matrices for damp
        call asseDamp(numeDof, listLoad, dampElem, dampAsse)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dampComputeVector
!
! Compute damping vector
!
! In  modalDamping     : parameters for modal damping
!
! --------------------------------------------------------------------------------------------------
    subroutine dampComputeVector()
!   ------------------------------------------------------------------------------------------------
! - Parameters

! - Local

!   ------------------------------------------------------------------------------------------------
!

!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module Damping_module
