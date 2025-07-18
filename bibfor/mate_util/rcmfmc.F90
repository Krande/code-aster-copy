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
subroutine rcmfmc(chmatz, chmacz, l_thm_, l_ther_, basename_, base)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rcmaco.h"
#include "asterfort/wkvect.h"
#include "asterfort/varc_prep.h"
!
    character(len=*), intent(in) :: chmatz
    character(len=*), intent(out) :: chmacz
    aster_logical, intent(in), optional :: l_thm_, l_ther_
    character(len=*), intent(in), optional :: basename_
    character(len=1), intent(in), optional :: base
!
! --------------------------------------------------------------------------------------------------
!
! Material
!
! Creation de la carte du materiau code a partir du champ_mater
!
! --------------------------------------------------------------------------------------------------
!
! In  chmate           : name of material field (CHAM_MATER)
! Out chmace           : name of CODED material field (CHAM_MATER)
! In  l_thm            : .true. if THM
! In  l_ther           : .true. if thermics
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbval, iret, jvale, igd, kk
    integer(kind=8) :: nbgrp, i, icompt, igrp, ingrp, nbcmp, j, k, nbmat
    integer(kind=8) :: inbmat
    character(len=1) :: bas
    character(len=4) :: knumat
    character(len=8) :: chmat, nomgd, basename
    character(len=19) :: codi
    character(len=19) :: chemat, chmace
    character(len=24) :: chmacegrp, chmacengrp
    character(len=8), pointer :: v_vale(:) => null()
    integer(kind=8), pointer :: v_desc(:) => null()
    aster_logical :: l_thm, l_ther
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    chmat = chmatz
!
    if (present(basename_)) then
        basename = basename_
    else
        basename = "&&MATECO"
    end if
!
    chemat = chmat//'.CHAMP_MAT'
    chmace = basename//'.MATE_CODE'
    chmacegrp = chmace//'.GRP'
    chmacengrp = chmace//'.NGRP'

!
    if (present(base)) then
        bas = base
    else
        bas = 'V'
    end if
!
    l_thm = ASTER_FALSE
    l_ther = ASTER_FALSE
    if (present(l_thm_)) then
        l_thm = l_thm_
    end if
    if (present(l_ther_)) then
        l_ther = l_ther_
    end if
!
    call exisd('CARTE', chmace, iret)
    if (iret .eq. 0) then
! ----- Preparation for external state variables (VARC)
        call varc_prep(chmat, l_thm)

! ----- Traitement du materiau par elements
        call jelira(chemat//'.VALE', 'LONMAX', nbval)
        call jeveuo(chemat//'.VALE', 'L', vk8=v_vale)
        call jeveuo(chemat//'.DESC', 'L', vi=v_desc)
        call jenuno(jexnum('&CATA.GD.NOMCMP', v_desc(1)), nomgd)
        call dismoi('NB_CMP_MAX', nomgd, 'GRANDEUR', repi=nbcmp)
        ASSERT(nbcmp .ge. 30)
        ASSERT((nbval/nbcmp)*nbcmp .eq. nbval)
        call copisd('CHAMP_GD', bas, chemat, chmace)
        call jedetr(chmace//'.VALE')
        nbgrp = nbval/nbcmp
        call wkvect(chmace//'.VALE', bas//' V I', nbgrp, jvale)
        call jenonu(jexnom('&CATA.GD.NOMGD', 'ADRSJEVE'), igd)
        call jeveuo(chmace//'.DESC', 'E', vi=v_desc)
        v_desc(1) = igd

! ----- Codage du materiau
        icompt = 0
        do i = 1, nbval
            if (v_vale(i) .ne. ' ') then
                icompt = icompt+1
            end if
        end do
        ASSERT(icompt .gt. 0)

        call jedetr(chmacegrp)
        call jedetr(chmacengrp)
        call wkvect(chmacegrp, bas//' V K8', icompt, igrp)
        call wkvect(chmacengrp, bas//' V I', nbgrp, ingrp)

        icompt = 0
        inbmat = 0
        do i = 1, nbgrp
            do j = 1, 26
                k = (i-1)*nbcmp+j
                if (v_vale(k) .eq. 'TREF=>') exit
                if (v_vale(k) .ne. ' ') then
                    zk8(igrp+icompt) = v_vale(k)
                    icompt = icompt+1
                    inbmat = inbmat+1
                end if
            end do
            zi(ingrp-1+i) = inbmat
            inbmat = 0
        end do

        codi = ' '
        call jeveuo(chmacegrp, 'L', igrp)
        call jeveuo(chmacengrp, 'L', ingrp)
        icompt = 0
        do kk = 1, nbgrp
            nbmat = zi(ingrp-1+kk)
            if (nbmat .ne. 0) then
                call rcmaco(chmat(1:8), chmacegrp, icompt, nbmat, kk, l_ther, basename, bas)
                call codent(kk, 'D0', knumat)
!       -- le nom du codi est celui du premier materiau du groupe kk
                codi(1:8) = basename
                codi(9:13) = '.'//knumat
!
                call jeveuo(codi//'.CODI', 'L', zi(jvale+kk-1))
                icompt = icompt+nbmat
            end if
        end do

    end if
    chmacz = chmace
!
    call jedema()
end subroutine
