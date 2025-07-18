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
subroutine focrr3(nomfon, resu, nopara, base, ier)
    implicit none
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsnopa.h"
#include "asterfort/rsutn1.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ier
    character(len=1) :: base
    character(len=16) :: nopara
    character(len=19) :: nomfon, resu
!     RECUPERATION D'UNE FONCTION DANS UNE STRUCTURE "RESULTAT"
!                                 PARAMETRE = F(VARIABLE D'ACCES)
!     ------------------------------------------------------------------
! VAR : NOMFON : NOM DE LA FONCTION
! IN  : RESU   : NOM DE LA STRUCTURE RESULTAT
! IN  : NOPARA : NOM DU PARAMETRE
! IN  : BASE   : BASE OU L'ON CREE LA FONCTION
! OUT : IER    : CODE RETOUR, = 0 : OK
!     ------------------------------------------------------------------
    integer(kind=8) :: nbordr, iret, kordr, lpro, lfon, lvar, iordr, nbacc, nbpar
    integer(kind=8) :: iad1, iad2, nbpt
    real(kind=8) :: rundf
    character(len=8) :: type
    character(len=16) :: nomacc
    character(len=19) :: knume
    character(len=16), pointer :: acces(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
    ier = 0
    knume = '&&FOCRR3.NUME_ORDR'
    rundf = r8vide()
!
!     --- RECUPERATION DES NUME_ORDRE FOURNIS PAR L'UTILISATEUR ---
!
    call rsutn1(resu, nopara, ' ', 1, knume, &
                nbordr)
    call jeveuo(knume, 'L', kordr)
!
!     --- RECUPERATION DE LA VARIABLE D'ACCES ---
!
    call rsnopa(resu, 0, '&&FOCRR3.VAR.ACCES', nbacc, nbpar)
    call jeexin('&&FOCRR3.VAR.ACCES', iret)
    if (iret .gt. 0) then
        call jeveuo('&&FOCRR3.VAR.ACCES', 'L', vk16=acces)
        nomacc = acces(1)
    else
        call utmess('F', 'UTILITAI2_4')
    end if
    call jedetr('&&FOCRR3.VAR.ACCES')
!
!     --- CREATION DE LA FONCTION SORTIE ---
!
!     --- REMPLISSAGE DU .PROL ---
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', base//' V K24', 6, lpro)
    zk24(lpro) = 'FONCTION'
    zk24(lpro+1) = 'LIN LIN '
    zk24(lpro+2) = nomacc(1:8)
    zk24(lpro+3) = nopara(1:8)
    zk24(lpro+4) = 'EE      '
    zk24(lpro+5) = nomfon
!
!
!   -- calcul du nombre de points de la fonction :
    nbpt = 0
    do iordr = 1, nbordr
        call rsadpa(resu, 'L', 1, nopara, zi(kordr+iordr-1), &
                    1, sjv=iad2, styp=type, istop=0)
        if (type(1:1) .ne. 'R') call utmess('F', 'UTILITAI2_6')
!
        if (zr(iad2) .eq. rundf) cycle
        nbpt = nbpt+1
    end do
!
!
!     --- REMPLISSAGE DU .VALE ---
    call wkvect(nomfon//'.VALE', base//' V R', 2*nbpt, lvar)
    lfon = lvar+nbpt
!
    nbpt = 0
    do iordr = 1, nbordr
        call rsadpa(resu, 'L', 1, nomacc, zi(kordr+iordr-1), &
                    1, sjv=iad1, styp=type)
        if (type(1:1) .ne. 'R') then
            call utmess('F', 'UTILITAI2_5')
        end if
!
        call rsadpa(resu, 'L', 1, nopara, zi(kordr+iordr-1), &
                    1, sjv=iad2, styp=type, istop=0)
        if (type(1:1) .ne. 'R') then
            call utmess('F', 'UTILITAI2_6')
        end if
!
        if (zr(iad2) .eq. rundf) cycle
!
        nbpt = nbpt+1
        zr(lvar+nbpt-1) = zr(iad1)
        zr(lfon+nbpt-1) = zr(iad2)
!
    end do
!
    call jedetr(knume)
!
    call jedema()
end subroutine
