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
subroutine rrssm2(neq, smhcr, smhci, smdir, smdii, &
                  idlexc, coef, valmi, valmr)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
    integer(kind=4) :: smhci(*), smhcr(*)
    integer(kind=8) :: idlexc(*)
    integer(kind=8) :: smdir(*), smdii(*)
    integer(kind=8) :: neq, idebl1, idebl2, kin1, kin2, iequa, ifinl1
    integer(kind=8) :: ifinl2, j2, k, i2, i1, j1, ind1
    real(kind=8) :: coef, valmi(*), valmr(*)
    integer(kind=8), pointer :: ind_lig(:) => null()
!--------------------------------------------------------
!
!
! --- CREATION DU TABLEAU D'INDIRECTION DES INDICES DES
! --- TERMES  D'UNE LIGNE DE LA MATRICE A COMBINER
! --- VERS LES INDICES DES TERMES DE LA MEME LIGNE
! --- DE LA MATRICE RESULTANTE :
!     ------------------------
    AS_ALLOCATE(vi=ind_lig, size=neq)
!
    idebl1 = 1
    idebl2 = 1
    kin1 = 0
    kin2 = 0
!
!
    do iequa = 1, neq
        ifinl1 = smdir(iequa)
        ifinl2 = smdii(iequa)
!
!       -- CALCUL DE .IND_LIG :
!       ------------------------
        do j2 = idebl2, ifinl2
            k = 0
            i2 = smhci(j2)
            do j1 = idebl1, ifinl1
                i1 = smhcr(j1)
                k = k+1
                if (i1 .eq. i2) goto 30
            end do
30          continue
            ind_lig(i2) = k
        end do
!
!
!       -- CUMUL DANS LA MATRICE RESULTAT :
!       ------------------------------------
        kin1 = idebl1-1
        do j2 = idebl2, ifinl2
            kin2 = kin2+1
            i2 = smhci(j2)
            ind1 = ind_lig(i2)
            valmr(kin1+ind1) = valmr(kin1+ind1)+coef*valmi(kin2)*(1-idlexc(i2))*(1-idlexc(ie&
                               &qua))
            ind_lig(i2) = 0
        end do
!
        idebl1 = smdir(iequa)+1
        idebl2 = smdii(iequa)+1
    end do
    AS_DEALLOCATE(vi=ind_lig)
!
end subroutine
