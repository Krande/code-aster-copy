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
subroutine numer2(nbLigr, listLigr, base, numeDofZ, &
                  numeDofOldZ, modeLocZ, modelZ, idenRelaZ_)
!
    implicit none
!
#include "asterc/cheksd.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/idensd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jemarq.h"
#include "asterfort/matdis.h"
#include "asterfort/nueffe.h"
#include "asterfort/nugllo.h"
#include "asterfort/promor.h"
!
    integer(kind=8), intent(in) :: nbLigr
    character(len=24), pointer :: listLigr(:)
    character(len=2), intent(in) :: base
    character(len=*), intent(inout) :: numeDofZ
    character(len=*), intent(in) :: numeDofOldZ, modeLocZ, modelZ
    character(len=*), optional, intent(in) :: idenRelaZ_
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering - Create objects
!
! --------------------------------------------------------------------------------------------------
!
! In  nbLigr         : number of LIGREL in list
! Ptr listLigr       : pointer to list of LIGREL
! In  base           : JEVEUX base to create objects
!                      base(2:2) => NUME_EQUA objects
!                      base(1:1) => NUME_DDL objects
! IO  numeDof        : name of numbering object (NUME_DDL)
! In  modeloc        : local mode for GRANDEUR numbering
! In  numeDofOld     : name of previous numeDof object
! In  idenRela       : name of object for identity relations between dof
!
! If numeDofOld is present
!   -> try to know if NUME_EQUA in numeDofOld can be reuse
!      In this case numeDof = numeDofOld
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: debug = ASTER_FALSE
    character(len=24), parameter :: renumSans = "SANS"
    character(len=19) :: numeEqua, numeEquaOld
    character(len=14) :: numeDof, numeDofOld, modeLoc
    character(len=24) :: idenRela
    character(len=3) :: matd
    character(len=8) :: model
    aster_logical :: l_matr_dist, printt
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    model = modelZ
    numeDof = numeDofZ
    modeLoc = modeLocZ
    numeDofOld = numeDofOldZ
    numeEqua = numeDof//'.NUME'
    numeEquaOld = numeDofOld//'.NUME'
    idenRela = ' '
    if (present(idenRelaZ_)) then
        idenRela = idenRelaZ_
    end if

! - Delete previous numbering
    call detrsd('NUME_DDL', numeDof)

! - For MATR_DISTRIBUE
    call matdis(matd, debug)
    ASSERT(matd .eq. 'OUI' .or. matd .eq. 'NON')
    l_matr_dist = matd .eq. 'OUI'

! - Create NUME_EQUA objects
    call nueffe(nbLigr, listLigr, base, numeDof, renumSans, model, &
                modeLocZ_=modeLoc, idenRelaZ_=idenRela)

! - Create NUML_EQUA objects
    if (l_matr_dist) then
        call nugllo(numeDofZ, base)
    end if

! - Trying to reuse old numeDof
    if (numeDofOld .ne. ' ') then
        if (idensd('NUME_EQUA', numeEqua, numeEquaOld)) then
            call detrsd('NUME_DDL', numeDof)
            call jedupo(numeDof//'     .ADNE', 'V', numeDofOld//'     .ADNE', .false._1)
            call jedupo(numeDof//'     .ADLI', 'V', numeDofOld//'     .ADLI', .false._1)
            call jedetr(numeDof//'     .ADLI')
            call jedetr(numeDof//'     .ADNE')
            numeDof = numeDofOld
        end if
    end if

! - Create matrix topology
    printt = modeLoc .eq. ' '
    call promor(numeDof, base(1:1), printt)

! - Cleaning
    call jedetr(numeDof//'     .ADLI')
    call jedetr(numeDof//'     .ADNE')
!
    numeDofZ = numeDof
!
    if (debug) then
        call cheksd(numeDofZ, 'SD_NUME_DDL', iret)
    end if
!
    call jedema()
end subroutine
