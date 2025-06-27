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
subroutine ddllag(nume, iddl, neq, lagr1, lagr2, ideb)

    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8), intent(in)           :: iddl, neq
    integer(kind=8), intent(out)          :: lagr1, lagr2
    character(len=*), intent(in)  :: nume
    integer(kind=8), intent(in), optional :: ideb

!
!     RECHERCHE LES DEUX LAGRANGES ASSOCIES AU DDL IDDL.
!     CE IDDL DDL EST BLOQUE ET ON NE LE VERIFIE PAS.
!     DANS LE CAS OU IDDL N'EST PAS BLOQUE, LAGR1=LAGR2=0
!
!     SI LE PARAMETRE OPTIONNEL IDEB EST RENSEIGNE, ON PEUT PARCOURIR
!     LA LISTE DES DDLS POTENTIELLEMENT PLUS EFFICACEMMENT EN PARTANT
!     D'UN POINT DE DEPART IDEB > 1. SI ON N'A RIEN TROUVE ON FAIT, EN
!     DERNIER RESSORT LE PARCOURS DE 1 A IDEB-1.
!
! IN  : NUME   : NOM D'UN NUME_DDL
! IN  : IDDL   : NUMERO D'UN DDL BLOQUE
! IN  : NEQ    : NOMBRE D'EQUATIONS
! OUT : LAGR1  : PREMIER LAGRANGE ASSOCIE
! OUT : LAGR2  : DEUXIEME LAGRANGE ASSOCIE
! ----------------------------------------------------------------------
    character(len=24) :: nomnu
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icas, icmp, inoe, nc, nn, idebu
    integer(kind=8), pointer :: deeq(:) => null()
!
!-----------------------------------------------------------------------
    call jemarq()

    idebu = 1
    if (present(ideb)) then
! VALEUR LICITE
        if ((ideb .ge. 1) .and. (ideb .le. neq)) then
            idebu = ideb
        end if
    end if
    lagr1 = 0
    lagr2 = 0
    nomnu(1:14) = nume
    nomnu(15:19) = '.NUME'
    call jeveuo(nomnu(1:19)//'.DEEQ', 'L', vi=deeq)
!
    inoe = deeq(1+(2*(iddl-1))+1-1)
    icmp = -deeq(1+(2*(iddl-1))+2-1)
    icas = 1
    do i = idebu, neq
        nn = deeq(1+(2*(i-1))+1-1)
        nc = deeq(1+(2*(i-1))+2-1)
        if (nn .eq. inoe .and. nc .eq. icmp) then
            if (icas .eq. 1) then
                lagr1 = i
                icas = 2
            else
                lagr2 = i
!               write(6,*)'<ddlag> standard',neq,inoe,icmp,lagr1,lagr2
                goto 999
            end if
        end if
    end do
! SOLUTION DE RATTRAPAGE AU CAS OU (SI IDEB >1)
    do i = 1, idebu-1
        nn = deeq(1+(2*(i-1))+1-1)
        nc = deeq(1+(2*(i-1))+2-1)
        if (nn .eq. inoe .and. nc .eq. icmp) then
            if (icas .eq. 1) then
                lagr1 = i
                icas = 2
            else
                lagr2 = i
!               write(6,*)'<ddlag> secours',neq,inoe,icmp,lagr1,lagr2
                goto 999
            end if
        end if
    end do
!
999 continue
    call jedema()
end subroutine
