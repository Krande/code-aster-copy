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
subroutine ve0124(typres)
    implicit none
    character(len=16), intent(inout) :: typres
! ----------------------------------------------------------------------
!     COMMANDE: CREA_RESU
!     VERIFICATION DE NIVEAU 1
! ----------------------------------------------------------------------
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/utmess.h"
    character(len=8) :: k8bid, resu
    character(len=16) :: type, oper
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid, iocc, k, n0, n1
    real(kind=8) :: r8bid
!-----------------------------------------------------------------------
    call getres(resu, type, oper)
    call getvtx(' ', 'TYPE_RESU', scal=typres, nbret=n1)
    call getfac('AFFE', iocc)
!
    if (typres .eq. 'EVOL_THER') then
        do k = 1, iocc
            call getvtx('AFFE', 'NOM_CAS', iocc=k, scal=k8bid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_7')
            end if
            call getvis('AFFE', 'NUME_MODE', iocc=k, scal=ibid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_8')
            end if
        end do
!
    else if (typres .eq. 'MULT_ELAS') then
        do k = 1, iocc
            call getvis('AFFE', 'NUME_MODE', iocc=k, scal=ibid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_9')
            end if
            call getvr8('AFFE', 'INST', iocc=k, scal=r8bid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_10')
            end if
            call getvid('AFFE', 'LIST_INST', iocc=k, scal=k8bid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_10')
            end if
        end do
!
    else if (typres .eq. 'FOURIER_ELAS') then
        do k = 1, iocc
            call getvtx('AFFE', 'NOM_CAS', iocc=k, scal=k8bid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_11')
            end if
            call getvr8('AFFE', 'INST', iocc=k, scal=r8bid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_12')
            end if
            call getvid('AFFE', 'LIST_INST', iocc=k, scal=k8bid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_12')
            end if
        end do
!
    else if (typres .eq. 'FOURIER_THER') then
        do k = 1, iocc
            call getvtx('AFFE', 'NOM_CAS', iocc=k, scal=k8bid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_13')
            end if
            call getvr8('AFFE', 'INST', iocc=k, scal=r8bid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_14')
            end if
            call getvid('AFFE', 'LIST_INST', iocc=k, scal=k8bid, nbret=n0)
            if (n0 .ne. 0) then
                call utmess('E', 'ALGORITH11_14')
            end if
        end do
    end if
!
    call getfac('PERM_CHAM', iocc)
    if (iocc .gt. 0) then
        call getvid(' ', 'RESU_INIT', nbval=0, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('E', 'ALGORITH11_15')
        end if
        call getvid(' ', 'MAILLAGE_INIT', nbval=0, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('E', 'ALGORITH11_16')
        end if
        call getvid(' ', 'RESU_FINAL', nbval=0, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('E', 'ALGORITH11_17')
        end if
        call getvid(' ', 'MAILLAGE_FINAL', nbval=0, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('E', 'ALGORITH11_18')
        end if
    end if
!
end subroutine
