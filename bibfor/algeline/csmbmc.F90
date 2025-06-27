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
subroutine csmbmc(nommat, neq, vsmb)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=*) :: nommat
    complex(kind=8) :: vsmb(*)
    integer(kind=8) :: neq
! BUT :
!-----------------------------------------------------------------------
!     FONCTIONS JEVEUX
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     VARIABLES LOCALES
!-----------------------------------------------------------------------
    integer(kind=8) :: ieq, iccid
    character(len=14) :: nu
    character(len=19) :: mat
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer :: ccid(:) => null()
    integer(kind=8), pointer :: nugl(:) => null()
!-----------------------------------------------------------------------
!     DEBUT
    call jemarq()
!-----------------------------------------------------------------------
    mat = nommat
!
    call jeveuo(mat//'.REFA', 'L', vk24=refa)
    if (refa(11) .eq. 'MATR_DISTR') then
        nu = refa(2) (1:14)
        call jeveuo(nu//'.NUML.NUGL', 'L', vi=nugl)
!
        call jeexin(mat//'.CCID', iccid)
!
        if (iccid .ne. 0) then
            call jeveuo(mat//'.CCID', 'L', vi=ccid)
            do ieq = 1, neq
!         SI LE DDL NE N'APPARTIENT PAS AU PROC COURANT ET QU'IL Y A
!         UNE CHARGE CINEMATIQUE DESSUS, ON MET LE SECOND MEMBRE A ZERO
!         SUR LE PROC COURANT POUR EVITER DES INTERFERENCES AVEC
!         LE PROC QUI POSSEDE EFFECTIVEMENT LE DDL BLOQUE
                if ((nugl(ieq) .eq. 0) .and. (ccid(ieq) .eq. 1)) then
                    vsmb(ieq) = 0.d0
                end if
            end do
        end if
    end if
!
    call jedema()
end subroutine
