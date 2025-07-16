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
subroutine rcadme(nommaz, phenom, nomres, valres, icodre, &
                  iarret)
    implicit none
#include "jeveux.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rccome.h"
#include "asterfort/rcvals.h"
#include "asterfort/tbexlr.h"
    character(len=*) :: nommaz, phenom, nomres
    integer(kind=8) :: icodre, iarret
    integer(kind=8) :: valres(*)
!
!     OBTENTION DES ADRESSES DES COMPOSANTES D'UN MATERIAU METALLURGIQUE
!
!     ARGUMENTS D'ENTREE:
!        NOMMAT : NOM UTILISATEUR DU MATERIAU
!     ARGUMENTS DE SORTIE:
!       VALRES : ADRESSE DU TRC IADTRC(1)=NBHIST IADTRC(2)=NBTRC
!       ICODRE : 0 SI ON A TROUVE, 1 SINON
! ----------------------------------------------------------------------
!
    integer(kind=8) :: iret, nbr, nbc, nbk, nbco, ik, nbcb1, nbcb2, nblb2
    integer(kind=8) :: nbhist, nbtrc
    character(len=11) :: k11
    character(len=8) :: nommat
    character(len=32) :: nomphe
    character(len=19) :: ch19, listr, noobrc
    real(kind=8), pointer :: vale(:) => null()
    character(len=16), pointer :: valk(:) => null()
! DEB ------------------------------------------------------------------
!
    call jemarq()
    nommat = nommaz
    nomphe = phenom
    k11 = ' '
!
    call rccome(nommat, nomphe, iret, k11_ind_nomrc=k11)
    if (iret .eq. 1) then
        icodre = 1
        goto 999
    end if
    noobrc = nommat//k11
    call jeexin(noobrc//'.VALR', iret)
    if (iret .eq. 0) then
        icodre = 1
        goto 999
    else
        call jelira(noobrc//'.VALR', 'LONUTI', nbr)
    end if
!
    call jeexin(noobrc//'.VALC', iret)
    if (iret .eq. 0) then
        icodre = 1
        goto 999
    else
        call jelira(noobrc//'.VALC', 'LONUTI', nbc)
    end if
!
    call jeexin(noobrc//'.VALK', iret)
    if (iret .eq. 0) then
        icodre = 1
        goto 999
    else
        call jeveuo(noobrc//'.VALK', 'L', vk16=valk)
        call jelira(noobrc//'.VALK', 'LONUTI', nbk)
    end if
!
    nbco = (nbk-nbr-nbc)/2
    do ik = 1, nbk
        if (nomres .eq. valk(ik)) then
            icodre = 0
            ch19 = valk(1+nbco+ik-1)
            listr = '&&RCADME.LR8'
            call tbexlr(ch19, listr, 'V')
            call jeveuo(listr//'.VALE', 'L', vr=vale)
            nbcb1 = nint(vale(2))
            nbhist = nint(vale(3))
            nbcb2 = nint(vale(1+1+2+nbcb1*nbhist))
            nblb2 = nint(vale(1+1+2+nbcb1*nbhist+1))
            nbtrc = nint(vale(1+1+2+nbcb1*nbhist+2+nbcb2*nblb2+1))
! --- NBHIST
            valres(1) = nbhist
! --- NBTRC
            valres(2) = nbtrc
! ---       MENAGE
            call detrsd('LISTR8', listr)
!
            goto 999
        end if
    end do
!
999 continue
!
    call rcvals(iarret, [icodre], 1, nomres)
!
    call jedema()
!
end subroutine
