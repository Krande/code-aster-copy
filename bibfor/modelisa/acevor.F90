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

subroutine acevor(nbocc, nlg, ier)
!
!
    implicit none
    integer :: nbocc, nlg, ier
!
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM
!     VERIFICATION DES MOTS CLES POUR LES ORIENTATIONS
!
! --------------------------------------------------------------------------------------------------
!
! IN  : NBOCC  : NOMBRE D'OCCURENCE
! OUT : NLG    : NOMBRE TOTAL DE GROUPE DE MAILLE
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
#include "asterc/getres.h"
#include "asterfort/codent.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/utmess.h"
! --------------------------------------------------------------------------------------------------
    integer :: ioc, jj, kk, nbcar, nbval, nc, ncar
    integer :: nco, ng, nv
    integer :: nval
! --------------------------------------------------------------------------------------------------
    parameter(nbcar=100, nbval=1000, nco=4)
    real(kind=8) :: val(nbval)
    character(len=6) :: kioc
    character(len=8) :: car(nbcar), nomu, carori(nco)
    character(len=16) :: cmd, concep
    character(len=24) :: valk(2)
    data carori/'VECT_Y  ', 'VECT_X_Y', 'ANGL_NAU', 'ANGL_VRI'/
! --------------------------------------------------------------------------------------------------
!
    call getres(nomu, concep, cmd)
    nlg = 0
!
    do ioc = 1, nbocc
        call codent(ioc, 'G', kioc)
        call getvtx('ORIENTATION', 'GROUP_MA', iocc=ioc, nbval=0, nbret=ng)
        call getvtx('ORIENTATION', 'CARA', iocc=ioc, nbval=0, nbret=nc)
        call getvtx('ORIENTATION', 'CARA', iocc=ioc, nbval=nbcar, vect=car, nbret=ncar)
        call getvr8('ORIENTATION', 'VALE', iocc=ioc, nbval=0, nbret=nv)
        call getvr8('ORIENTATION', 'VALE', iocc=ioc, nbval=nbval, vect=val, nbret=nval)
!
        if (ioc .eq. 1) then
            if (nv .eq. 0) then
                call utmess('E', 'MODELISA_57')
                ier = ier+1
            end if
            if (nc .eq. 0) then
                call utmess('E', 'MODELISA_58')
                ier = ier+1
            end if
        end if
!       CARA
        kk = 0
        if (ncar .gt. 0) then
            if (nval .eq. 0) then
                call utmess('E', 'MODELISA_59', sk=kioc)
                ier = ier+1
            end if
            do jj = 1, nco
                if (car(1) .eq. carori(jj)) kk = jj
            end do
        end if
!       VALE
        if (nval .gt. 0) then
            if ((kk .eq. 1 .and. nval .ne. 3) .or. (kk .eq. 2 .and. nval .ne. 6) .or. &
                (kk .eq. 3 .and. nval .ne. 3) .or. (kk .eq. 4 .and. nval .ne. 1)) then
                valk(1) = kioc
                valk(2) = carori(kk)
                call utmess('E', 'MODELISA_60', nk=2, valk=valk)
                ier = ier+1
            end if
        end if
!
        nlg = max(nlg, -ng)
    end do
!
end subroutine
