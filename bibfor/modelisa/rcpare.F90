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

subroutine rcpare(nommat, pheno, para, icodre)
    implicit none
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
    character(len=*) :: nommat, pheno, para
    integer(kind=8) :: icodre
! ----------------------------------------------------------------------
!     VERIFICATION DE LA PRESENCE D'UNE CARACTERISTIQUE DANS UN
!     COMPORTEMENT DONNE
!
!     ARGUMENTS D'ENTREE:
!        NOMMAT  : NOM DU MATERIAU
!        PHENO   : NOM DE LA LOI DE COMPORTEMENT
!        PARA    : NOM DU PARAMETRE
!     ARGUMENTS DE SORTIE:
!     ICODRE : POUR CHAQUE RESULTAT, 0 SI ON A TROUVE, 1 SINON
!
!
!
! ----------------------------------------------------------------------
! DEB ------------------------------------------------------------------
    character(len=6) ::  k6
    character(len=8) ::  nomma2
    character(len=32) :: pheno2
    character(len=32) :: ncomr, ncomc, ncomk, ncomp
    integer(kind=8) :: nbpar, nbr, nbc, nbk, nbcomp, icomp
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, ipar
!-----------------------------------------------------------------------
    icodre = 1
    pheno2 = pheno
    nomma2 = nommat
!
    ncomp = nommat//'.MATERIAU.NOMRC         '
    call jelira(ncomp, 'LONUTI', nbcomp)
    call jeveuo(ncomp, 'L', icomp)
    do i = 1, nbcomp
        if (pheno2 .eq. zk32(icomp+i-1)) then
            call codent(i, 'D0', k6)
            ncomr = nomma2//'.CPT.'//k6//'.VALR        '
            ncomc = nomma2//'.CPT.'//k6//'.VALC        '
            ncomk = nomma2//'.CPT.'//k6//'.VALK        '
            call jelira(ncomr, 'LONUTI', nbr)
            call jelira(ncomc, 'LONUTI', nbc)
            call jelira(ncomk, 'LONUTI', nbk)
            call jeveuo(ncomk, 'L', ipar)
            nbpar = nbr+nbc+nbk/2
            do j = 1, nbpar
                if (para .eq. zk16(ipar+j-1)) then
                    icodre = 0
                end if
            end do
        end if
    end do
! FIN ------------------------------------------------------------------
end subroutine
