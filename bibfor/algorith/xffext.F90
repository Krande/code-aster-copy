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
subroutine xffext(jinfo, nfon, nmafon, listpt, ptextr, &
                  nbptex)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: jinfo, nfon, nmafon, nbptex
    character(len=19) :: listpt, ptextr
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM
!
!              RECHERCHE DES POINTS EXTREMITES DU FOND DE FISSURE
!
! ----------------------------------------------------------------------
!
!
! IN  JINFO      : ADRESSE DU VECTEUR INFO DE LA SD
! IN  NFON       : NOMBRE DE POINTS AU FOND DE FISSURE
! IN  NMAFON     : NOMBRE DE MAILLES CONNECTEES AU FOND
! IN  LISTPT     : VECTEUR CONTENANT LES INDICES DES POINTS DU FOND PAR
!                  MAILLE
!
! IN/OUT PTEXTR  : VECTEUR CONTENANT LES INDICES DES POINTS EXTREMITES
!                  DU OU DES FONDS DE FISSURE
! OUT    NBPTEX  : NOMBRE D'EXTREMITES DU OU DES FONDS
!
!
    integer(kind=8) :: ima, ipt, jlistp, jptext, nocc, ptasso
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    call wkvect(ptextr, 'V V I', nfon, jptext)
!
    call jeveuo(listpt, 'L', jlistp)
!
    nbptex = 0
!
    do ipt = 1, nfon
!       COMPTEUR DU NOMBRE D'OCCURENCES DU POINT IPT DANS LE TABLEAU
!       LISTPT
        nocc = 0
!       INDICE DU POINT ASSOCIE AU POINT IPT
        ptasso = 0
!
        do ima = 1, nmafon
            if ((zi(jlistp-1+2*(ima-1)+1) .eq. ipt) .and. (zi(jlistp-1+2*(ima-1)+2) .ne. 0) &
                .and. (ptasso .ne. zi(jlistp-1+2*(ima-1)+2))) then
!
                ptasso = zi(jlistp-1+2*(ima-1)+2)
                nocc = nocc+1
!
            elseif ((zi(jlistp-1+2*(ima-1)+2) .eq. ipt) .and. ( &
                    ptasso .ne. zi(jlistp-1+2*(ima-1)+1))) then
!
                ptasso = zi(jlistp-1+2*(ima-1)+1)
                nocc = nocc+1
            end if
!
            if (nocc .eq. 2) goto 10
!
        end do
!
!       UNE SEULE OCCURENCE DU POINT IPT:
!       C'EST UN POINT EXTREMITE DU FOND DE FISSURE
        if (nocc .eq. 1) then
            nbptex = nbptex+1
            zi(jptext-1+nbptex) = ipt
        end if
10      continue
    end do
!
!     ON DOIT AVOIR UN NOMBRE PAIR D'EXTREMITES
    ASSERT(mod(nbptex, 2) .eq. 0)
!
    zk16(jinfo-1+3) = 'OUVERT'
!     CAS D'UN FOND FERME SI ABSCENCE D'EXTREMITES
    if (nbptex .eq. 0) then
        zk16(jinfo-1+3) = 'FERME'
        zi(jptext-1+1) = 1
    end if
!
    call jedema()
end subroutine
