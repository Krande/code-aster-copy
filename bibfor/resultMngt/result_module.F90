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
! ==================================================================================================
!
! Module for management of result datastructures
!
! ==================================================================================================
!
module result_module
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: rsCopyPara
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsnopa.h"
#include "jeveux.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! rsCopyPara
!
! Copy parameters from one result to another one
!
! In  resultIn          : name of datastructure for input results
! In  resultOut         : name of datastructure for output results
! In  nbStore           : number of storing indexes
! Ptr listStore         : pointer to list of storing indexes
!
! --------------------------------------------------------------------------------------------------
    subroutine rsCopyPara(resultInZ, resultOutZ, nbStore, listStore)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: resultInZ, resultOutZ
        integer(kind=8), intent(in) :: nbStore
        integer(kind=8), pointer :: listStore(:)
! ----- Local
        character(len=24), parameter :: paraJvName = '&&CCBCOP.NOMS_PARA '
        integer(kind=8) :: nbParaAccess, nbPara, nbParaTotal, jvPara, iParaTotal
        integer(kind=8) :: numeStore, iStore
        character(len=19) :: resultIn, resultOut
        character(len=8) :: paraType
        integer(kind=8) :: jvResultIn, jvResultOut
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

! ----- Initializations
        resultIn = resultInZ
        resultOut = resultOutZ

! ----- Acces to parameters
        call rsnopa(resultIn, 2, paraJvName, nbParaAccess, nbPara)
        nbParaTotal = nbParaAccess+nbPara
        call jeveuo(paraJvName, 'L', jvPara)

! ----- Copy parameters
        do iStore = 1, nbStore
            numeStore = listStore(iStore)
            do iParaTotal = 1, nbParaTotal
                call rsadpa(resultIn, 'L', 1, zk16(jvPara+iParaTotal-1), numeStore, &
                            1, sjv=jvResultIn, styp=paraType, istop=0)
                call rsadpa(resultOut, 'E', 1, zk16(jvPara+iParaTotal-1), numeStore, &
                            1, sjv=jvResultOut, styp=paraType)
                if (paraType(1:1) .eq. 'I') then
                    zi(jvResultOut) = zi(jvResultIn)
                else if (paraType(1:1) .eq. 'R') then
                    zr(jvResultOut) = zr(jvResultIn)
                else if (paraType(1:1) .eq. 'C') then
                    zc(jvResultOut) = zc(jvResultIn)
                else if (paraType(1:3) .eq. 'K80') then
                    zk80(jvResultOut) = zk80(jvResultIn)
                else if (paraType(1:3) .eq. 'K32') then
                    zk32(jvResultOut) = zk32(jvResultIn)
                else if (paraType(1:3) .eq. 'K24') then
                    zk24(jvResultOut) = zk24(jvResultIn)
                else if (paraType(1:3) .eq. 'K16') then
                    zk16(jvResultOut) = zk16(jvResultIn)
                else if (paraType(1:2) .eq. 'K8') then
                    zk8(jvResultOut) = zk8(jvResultIn)
                end if
            end do
        end do
        call jedetr(paraJvName)
!
        call jedema()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module result_module
