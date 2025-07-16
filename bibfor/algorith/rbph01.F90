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
subroutine rbph01(trange, nbcham, typea, itresu, nfonct, &
                  basemo, typref, typbas, tousno, multap, i_cham)

    use DynaGene_module

    implicit none
#include "asterf_types.h"
#include "asterfort/getvtx.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbcham, itresu(*), nfonct
    character(len=8) :: basemo
    character(len=16) :: typea(*), typbas(*)
    character(len=19) :: trange, typref(*)
    aster_logical :: tousno, multap
    integer(kind=8), dimension(8), intent(out), optional :: i_cham
!     OPERATEUR REST_BASE_PHYS
!               TRAITEMENT DES MOTS CLES "TOUT_CHAM" ET "NOM_CHAM"
!     ------------------------------------------------------------------
    integer(kind=8) :: n1, i, iret
    character(len=8) :: blanc, mode
    character(len=16) :: champ(8)
    character(len=19) :: nomcha
    type(DynaGene) :: dyna_gene
!     ------------------------------------------------------------------
    data blanc/'        '/
!     ------------------------------------------------------------------
!
    mode = basemo
    call dyna_gene%init(trange(1:8))
!
    champ(1) = ' '
    call getvtx(' ', 'TOUT_CHAM', scal=champ(1), nbret=n1)
!
    if (champ(1) .eq. 'OUI') then
        nbcham = 3
        typea(1) = 'DEPL            '
        typea(2) = 'VITE            '
        typea(3) = 'ACCE            '
        call dyna_gene%has_field(dyna_gene%depl, iret)
        if (iret .eq. 0) then
            call utmess('F', 'ALGORITH10_11')
        else
            if (present(i_cham)) then
                i_cham(1) = dyna_gene%depl
            else
                call jeveuo(trange//'.DEPL', 'L', itresu(1))
            end if
        end if
!
        call dyna_gene%has_field(dyna_gene%vite, iret)
        if (iret .eq. 0) then
            call utmess('F', 'ALGORITH10_12')
        else
            if (present(i_cham)) then
                i_cham(2) = dyna_gene%vite
            else
                call jeveuo(trange//'.VITE', 'L', itresu(2))
            end if
        end if
        call dyna_gene%has_field(dyna_gene%acce, iret)
        if (iret .eq. 0) then
            call utmess('F', 'ALGORITH10_13')
        else
            if (present(i_cham)) then
                i_cham(3) = dyna_gene%acce
            else
                call jeveuo(trange//'.ACCE', 'L', itresu(3))
            end if
        end if
        if (nfonct .ne. 0) then
            nbcham = 4
            if (present(i_cham)) then
                i_cham(4) = dyna_gene%acce
            else
                itresu(4) = itresu(3)
            end if
            typea(4) = 'ACCE_ABSOLU     '
        end if
        if (mode .eq. blanc) then
            typref(1) = ' '
            typref(2) = ' '
            typref(3) = ' '
            typref(4) = ' '
        else
            call rsexch(' ', basemo, 'DEPL', 1, nomcha, &
                        iret)
            typref(1) = nomcha
            typref(2) = nomcha
            typref(3) = nomcha
            typref(4) = nomcha
        end if
        typbas(1) = 'DEPL'
        typbas(2) = 'DEPL'
        typbas(3) = 'DEPL'
        typbas(4) = 'DEPL'
!
    else
!
        call getvtx(' ', 'NOM_CHAM', nbval=0, nbret=n1)
        nbcham = -n1
        call getvtx(' ', 'NOM_CHAM', nbval=nbcham, vect=champ, nbret=n1)
!
        do i = 1, nbcham
            if (champ(i) .eq. 'DEPL') then
                typea(i) = 'DEPL'
                call dyna_gene%has_field(dyna_gene%depl, iret)
                if (iret .eq. 0) then
                    call utmess('F', 'ALGORITH10_11')
                else
                    if (present(i_cham)) then
                        i_cham(i) = dyna_gene%depl
                    else
                        call jeveuo(trange//'.DEPL', 'L', itresu(i))
                    end if
                end if
                if (mode .eq. blanc) then
                    typref(i) = ' '
                else
                    call rsexch(' ', basemo, typea(i), 1, nomcha, &
                                iret)
                    typref(i) = nomcha
                end if
                typbas(i) = 'DEPL'
!
            else if (champ(i) .eq. 'VITE') then
                typea(i) = 'VITE'
                call dyna_gene%has_field(dyna_gene%vite, iret)
                if (iret .eq. 0) then
                    call utmess('F', 'ALGORITH10_12')
                else
                    if (present(i_cham)) then
                        i_cham(i) = dyna_gene%vite
                    else
                        call jeveuo(trange//'.VITE', 'L', itresu(i))
                    end if
                end if
                if (mode .eq. blanc) then
                    typref(i) = ' '
                else
                    call rsexch(' ', basemo, 'DEPL', 1, nomcha, &
                                iret)
                    typref(i) = nomcha
                end if
                typbas(i) = 'DEPL'
!
            else if (champ(i) .eq. 'ACCE') then
                typea(i) = 'ACCE'
                call dyna_gene%has_field(dyna_gene%acce, iret)
                if (iret .eq. 0) then
                    call utmess('F', 'ALGORITH10_13')
                else
                    if (present(i_cham)) then
                        i_cham(i) = dyna_gene%acce
                    else
                        call jeveuo(trange//'.ACCE', 'L', itresu(i))
                    end if
                end if
                if (mode .eq. blanc) then
                    typref(i) = ' '
                else
                    call rsexch(' ', basemo, 'DEPL', 1, nomcha, &
                                iret)
                    typref(i) = nomcha
                end if
                typbas(i) = 'DEPL'
!
            else if (champ(i) .eq. 'ACCE_ABSOLU') then
                typea(i) = 'ACCE_ABSOLU'
                call dyna_gene%has_field(dyna_gene%acce, iret)
                if (iret .eq. 0) then
                    call utmess('F', 'ALGORITH10_13')
                else
                    if (present(i_cham)) then
                        i_cham(i) = dyna_gene%acce
                    else
                        call jeveuo(trange//'.ACCE', 'L', itresu(i))
                    end if
                end if
                if (mode .eq. blanc) then
                    typref(i) = ' '
                else
                    call rsexch(' ', basemo, 'DEPL', 1, nomcha, &
                                iret)
                    typref(i) = nomcha
                end if
                typbas(i) = 'DEPL'
!
            elseif (champ(i) .eq. 'FORC_NODA' .or. champ(i) .eq. &
                    'REAC_NODA') then
                typea(i) = champ(i)
                call dyna_gene%has_field(dyna_gene%depl, iret)
                if (iret .eq. 0) then
                    call utmess('F', 'ALGORITH10_11')
                else
                    if (present(i_cham)) then
                        i_cham(i) = dyna_gene%depl
                    else
                        call jeveuo(trange//'.DEPL', 'L', itresu(i))
                    end if
                end if
                if (multap) then
                    call utmess('F', 'ALGORITH10_14')
                end if
                if (mode .eq. blanc) then
                    call utmess('F', 'ALGORITH10_15')
                else
                    call rsexch('F', basemo, typea(i), 1, nomcha, &
                                iret)
                    typref(i) = nomcha
                end if
                typbas(i) = typea(i)
!
            else
                typea(i) = champ(i)
                call dyna_gene%has_field(dyna_gene%depl, iret)
                if (iret .eq. 0) then
                    call utmess('F', 'ALGORITH10_11')
                else
                    if (present(i_cham)) then
                        i_cham(i) = dyna_gene%depl
                    else
                        call jeveuo(trange//'.DEPL', 'L', itresu(i))
                    end if
                end if
                if (.not. tousno) then
                    call utmess('F', 'ALGORITH10_17', sk=typea(i))
                end if
                if (multap) then
                    call utmess('F', 'ALGORITH10_14')
                end if
                if (mode .eq. blanc) then
                    call utmess('F', 'ALGORITH10_15')
                else
                    call rsexch('F', basemo, typea(i), 1, nomcha, &
                                iret)
                    typref(i) = nomcha
                end if
                typbas(i) = typea(i)
!
            end if
        end do
    end if

    call dyna_gene%free()
!
end subroutine
