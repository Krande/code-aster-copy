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

subroutine focrrs(nomfon, resu, base, nomcha, maille, &
                  noeud, cmp, npoint, nusp, ivari, nomvari, ier)
    implicit none
#include "jeveux.h"
#include "asterfort/focrr0.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsutn2.h"
    integer(kind=8), intent(in) :: npoint, nusp, ivari
    character(len=1), intent(in) :: base
    character(len=8), intent(in) :: maille, noeud, cmp
    character(len=16), intent(in) :: nomcha
    character(len=16), intent(in) :: nomvari
    character(len=19), intent(in) :: nomfon, resu
    integer(kind=8), intent(out) :: ier
!     RECUPERATION D'UNE FONCTION DANS UNE STRUCTURE "RESULTAT"
!     ------------------------------------------------------------------
! VAR : NOMFON : NOM DE LA FONCTION
! IN  : RESU   : NOM DE LA STRUCTURE RESULTAT
! IN  : BASE   : BASE OU L'ON CREE LA FONCTION
! IN  : NOMCHA : NOM DU CHAMP
! IN  : NOEUD  : NOEUD
! IN  : MAILLE : MAILLE
! IN  : CMP    : COMPOSANTE
! IN  : NPOINT : NUMERO DU POINT ( CAS DES CHAM_ELEMS )
! IN  : NUSP   : NUMERO DU SOUS-POINT ( CAS DES CHAM_ELEMS )
! IN  : IVARI  : NUMERO DE LA CMP (POUR VARI_R)
! IN  : NOMVARI: NOM DE LA CMP (POUR VARI_R)
! OUT : IER    : CODE RETOUR, = 0 : OK
!     ------------------------------------------------------------------
    integer(kind=8) :: nbordr, lordr
    character(len=8) :: interp
    character(len=19) :: knume
!     ------------------------------------------------------------------
!
    call jemarq()
    ier = 0
    knume = '&&FOCRRS.NUME_ORDR'
!
!     --- RECUPERATION DES NUME_ORDRE FOURNIS PAR L'UTILISATEUR ---
!
    call rsutn2(resu, nomcha, ' ', 1, knume, &
                nbordr)
    call jeveuo(knume, 'L', lordr)
!
    interp = 'LIN LIN '
!
    call focrr0(nomfon, interp, base, resu, nomcha, &
                maille, noeud, cmp, npoint, nusp, &
                ivari, nomvari, nbordr, zi(lordr))
!
    call jedetr(knume)
!
    call jedema()
end subroutine
