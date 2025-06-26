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
subroutine csmbr8(nommat, ccll, ccii, neq, vcine, &
                  vsmb)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    character(len=*) :: nommat
    real(kind=8) :: vsmb(*), vcine(*)
    integer(kind=8) :: ccll(*), ccii(*), neq
! BUT : CALCUL DE LA CONTRIBUTION AU SECOND MEMBRE DES DDLS IMPOSES
!       LORSQU'ILS SONT TRAITEES PAR ELIMINATION (CAS REEL)
! C.F. EXPLICATIONS DANS LA ROUTINE CSMBGG
!-----------------------------------------------------------------------
! IN  NOMMAT K19 : NOM DE LA MATR_ASSE
! IN  CCLL   I(*): TABLEAU .CCLL DE LA MATRICE
! IN  CCII   I(*): TABLEAU .CCII DE LA MATRICE
! IN  NEQ    I   : NOMBRE D'EQUATIONS
! VAR VSMB   R(*): VECTEUR SECOND MEMBRE
! IN  VCINE  R(*): VECTEUR DE CHARGEMENT CINEMATIQUE ( LE U0 DE U = U0
!                 SUR G AVEC VCINE = 0 EN DEHORS DE G )
!-----------------------------------------------------------------------
!     FONCTIONS JEVEUX
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     VARIABLES LOCALES
!-----------------------------------------------------------------------
    integer(kind=8) :: nelim, ielim, ieq, j, ieqg
    integer(kind=8) :: deciel, kterm, nterm, imatd
    real(kind=8) :: coef
    character(len=14) :: nu
    character(len=19) :: mat
    integer(kind=8), pointer :: ccid(:) => null()
    integer(kind=8), pointer :: nugl(:) => null()
    real(kind=8), pointer :: ccva(:) => null()
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer :: nulg(:) => null()
!-----------------------------------------------------------------------
!     DEBUT
    call jemarq()
!-----------------------------------------------------------------------
    mat = nommat
!
    call jeveuo(mat//'.CCVA', 'L', vr=ccva)
    call jelira(mat//'.CCLL', 'LONMAX', nelim)
    nelim = nelim/3
!
    call jeveuo(mat//'.REFA', 'L', vk24=refa)
    if (refa(11) .eq. 'MATR_DISTR') then
        imatd = 1
        nu = refa(2) (1:14)
        call jeveuo(nu//'.NUML.NULG', 'L', vi=nulg)
        call jeveuo(nu//'.NUML.NUGL', 'L', vi=nugl)
    else
        imatd = 0
    end if
!
    do ielim = 1, nelim
        ieq = ccll(3*(ielim-1)+1)
        nterm = ccll(3*(ielim-1)+2)
        deciel = ccll(3*(ielim-1)+3)
!
        if (imatd .eq. 0) then
            ieqg = ieq
        else
            ieqg = nulg(ieq)
        end if
        coef = vcine(ieqg)
!
        if (coef .ne. 0.d0) then
            do kterm = 1, nterm
                if (imatd .eq. 0) then
                    j = ccii(deciel+kterm)
                else
                    j = nulg(ccii(deciel+kterm))
                end if
                vsmb(j) = vsmb(j)-coef*ccva(deciel+kterm)
            end do
        end if
!
    end do
    call jelibe(mat//'.CCVA')
!
    if (imatd .ne. 0) then
        do ieq = 1, neq
            if (nugl(ieq) .eq. 0) vcine(ieq) = 0.d0
        end do
    end if
!
!
    call jeveuo(mat//'.CCID', 'L', vi=ccid)
    do ieq = 1, neq
        if (ccid(ieq) .eq. 1) then
            vsmb(ieq) = vcine(ieq)
        else
            if (vcine(ieq) .ne. 0.d0) then
                call utmess('F', 'ALGELINE_32')
            end if
        end if
!
    end do
!
    call jedema()
end subroutine
