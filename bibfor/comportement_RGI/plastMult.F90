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
subroutine plastMult(na, nf0, ngf, x, irr, nc, indic2, ipla2, imax, ig, &
                     a, b, dgfa_ds, goto20)
! person_in_charge: etienne.grimal@edf.fr
!-----------------------------------------------------------------------
!       verif positivite des multiplicateurs plastiques
!-----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"

!   declaration des arguments
    integer(kind=8), intent(in) :: nf0, irr, nc, imax, ngf
    integer(kind=8), intent(inout) :: na, ipla2, ig(nc)
    real(kind=8), intent(inout) :: x(ngf), a(ngf, ngf+1), b(ngf), dgfa_ds(nc, 6)
    aster_logical, intent(inout) :: indic2
    aster_logical, intent(out) :: goto20
! ----------------------------------------------------------------------
!   variables internes
    integer(kind=8) :: nc0
    parameter(nc0=10)
    integer(kind=8) :: nsupr, i, j, k, supr(nc0), nared
! ----------------------------------------------------------------------
    ASSERT(nc0 .eq. nc)
    supr(:) = 0
    goto20 = ASTER_FALSE
    if (na .ne. 0) then
!       initialisation du compteur de nbre de lignes a supprimer
        nsupr = 0
        do i = (nf0+1), (nf0+na)
!   test positivite des multiplicateurs
            if ((x(i) .lt. 0.d0) .and. (irr .eq. 0)) then
!   actualisation des numero de ligne a supprimer
!   pour la mise a zero des multiplicateurs negatifs
                nsupr = nsupr+1
                supr(nsupr) = i
            end if
        end do
        if (nsupr .gt. 0) then
!    indicateur de reduction matrice de couplage
            indic2 = .true.
!    compteur reduction
            ipla2 = ipla2+1
            if (ipla2 .le. imax) then
                nared = na
!     reduction du systeme lineaire des couplages
!     on remonte toute les lignes au dessous de celles supprimees
                do i = 1, nsupr
!     decallage vers le haut des lignes du dessous
                    do j = supr(i), nf0+(nared-1)
                        ig(j-nf0) = ig(j-nf0+1)
!       on boucle sur les colonnes
                        do k = 1, nf0+nared
                            a(j, k) = a(j+1, k)
                        end do
!       pareil au second membre
                        b(j) = b(j+1)
                    end do
!     decalage des colonnes vers la gauche
                    do j = supr(i), nf0+(nared-1)
                        do k = 1, nf0+(nared-1)
                            a(k, j) = a(k, j+1)
                        end do
                    end do
!     decalage des derives des fonctions de charge
                    do j = supr(i)-nf0, nf0+(nared-1)-nf0
                        do k = 1, 6
                            dgfa_ds(j, k) = dgfa_ds(j+1, k)
                        end do
                    end do
!     mise a jour de la taille de la matrice
                    nared = nared-1
!     mise a jour des numeros de lignes a supprimer
                    do j = i, nsupr
!       comme on vient d eliminer une ligne
!       les lignes restantes sont remontees de 1
                        supr(j) = supr(j)-1
                    end do
                end do
!     resolution du systeme reduit
                na = nared
                goto20 = ASTER_TRUE
                goto 20
            else
                call utmess('F', 'COMPOR3_31', si=ipla2)
            end if
        end if
    end if
20  continue
end subroutine
