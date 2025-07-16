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
subroutine cbsint(load, mesh)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
!
    character(len=8), intent(in) :: load, mesh
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load PRE_SIGM
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  mesh             : mesh
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'PRE_SIGM'
    character(len=5), parameter :: param = 'SIINT'
    integer(kind=8) :: ibid, nbOcc, ncmp
    character(len=19) :: carte
    character(len=24) :: chsig
    character(len=8), pointer :: valv(:) => null()
    character(len=8), pointer :: vncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getfac(keywordFact, nbOcc)
    if (nbOcc .ne. 0) then
        carte = load//'.CHME.'//param
        call alcart('G', carte, mesh, 'NEUT_K8')
        call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
        call jeveuo(carte//'.VALV', 'E', vk8=valv)
        ncmp = 1
        vncmp(1) = 'Z1'
        call getvid(keywordFact, 'SIGM', iocc=1, scal=chsig, nbret=ibid)
        valv(1) = chsig(1:8)
        call nocart(carte, 1, ncmp)
    end if

    call jedema()
end subroutine
