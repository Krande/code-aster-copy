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
subroutine rcvalc(jmat, phenom, nbres, nomres, valres, &
                  icodre, iarret)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/rcvals.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: imat, nbres, jmat, nbmat
    character(len=*) :: phenom, nomres(nbres)
    integer(kind=8) :: iarret
    integer(kind=8) :: icodre(nbres)
    complex(kind=8) :: valres(nbres)
! ----------------------------------------------------------------------
!     OBTENTION DE LA VALEUR VALRES C D'UN "ELEMENT" D'UNE RELATION DE
!     COMPORTEMENT D'UN MATERIAU DONNE (NOUVELLE FORMULE RAPIDE)
!
!     ARGUMENTS D'ENTREE:
!        IMAT   : ADRESSE DU MATERIAU CODE
!        PHENOM : NOM DU PHENOMENE
!        NBRES  : NOMBRE DE RESULTATS
!        NOMRES : NOM DES RESULTATS (EX: E,NU,... )
!                 TELS QU'IL FIGURENT DANS LA COMMANDE MATERIAU
!     ARGUMENTS DE SORTIE:
!     VALRES : VALEURS DES RESULTATS APRES RECUPERATION ET INTERPOLATION
!     ICODRE : POUR CHAQUE RESULTAT, 0 SI ON A TROUVE, 1 SINON
!
!
!
!
!
!
    character(len=32) :: nomphe
! ----------------------------------------------------------------------
! PARAMETER ASSOCIE AU MATERIAU CODE
! DEB ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: icomp, idf, ik, ipi, ir, ires, ivalc
    integer(kind=8) :: ivalk, nbc, nbf, nbobj, nbr, nbt
!-----------------------------------------------------------------------
    nbmat = zi(jmat)
    ASSERT(nbmat .eq. 1)
    imat = jmat+zi(jmat+nbmat+1)
!
!
    do ires = 1, nbres
        icodre(ires) = 1
    end do
    nomphe = phenom
    do icomp = 1, zi(imat+1)
        if (nomphe .eq. zk32(zi(imat)+icomp-1)) then
            ipi = zi(imat+2+icomp-1)
            goto 11
        end if
    end do
    call utmess('A', 'ELEMENTS2_63')
    goto 999
11  continue
!
    nbobj = 0
    nbr = zi(ipi)
    nbc = zi(ipi+1)
    ivalk = zi(ipi+3)
    ivalc = zi(ipi+5)
    nbt = nbr+nbc
    do ir = 1, nbt
        do ires = 1, nbres
            if (nomres(ires) .eq. zk16(ivalk+ir-1)) then
                valres(ires) = zc(ivalc-1+ir)
                icodre(ires) = 0
                nbobj = nbobj+1
            end if
        end do
    end do
    if (nbobj .ne. nbres) then
        idf = zi(ipi)+zi(ipi+1)
        nbf = zi(ipi+2)
        do ires = 1, nbres
            do ik = 1, nbf
                if (nomres(ires) .eq. zk16(ivalk+idf+ik-1)) then
                    call utmess('F', 'MODELISA6_93')
!              CALL FOINTA (IFON,NBPAR,NOMPAR,VALPAR,VALRES(IRES))
                    icodre(ires) = 0
                end if
            end do
        end do
    end if
999 continue
!
    call rcvals(iarret, icodre, nbres, nomres)
!
end subroutine
