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
subroutine vpordc(type, iordre, nbpro, valpro, vecpro, &
                  neq)
    implicit none
    integer(kind=8) :: type, nbpro, neq
    complex(kind=8) :: valpro(*), vecpro(neq, nbpro)
!     TRIE DES VALEURS (ET DES VECTEURS) PROPRES COMPLEXES
!     PAR ORDRE CROISSANT
!     ------------------------------------------------------------------
! IN  TYPE   : IS : TYPE DU TRI SUR LES VALEURS.
!        * SI TYPE = 0  TRI EN VALEUR RELATIVE
!        * SI TYPE = 1  TRI EN VALEUR ABSOLUE
! IN  IORDRE : IS : ORDRE DU TRI SUR LES VALEURS.
!        * SI IORDRE = 0  TRI PAR ORDRE CROISSANT
!        * SI IORDRE = 1  TRI PAR ORDRE DECROISSANT
! IN  NBPRO  : IS : NOMBRE DE VALEUR PROPRE
!     VALPRO : R8 : TABLEAU DES VALEURS PROPRES
!     VECPRO : R8 : MATRICE DES VECTEURS PROPRES
!     NEQ    : IS : NOMBRE D'EQUATIONS
!                 SI NEQ < NBPRO ALORS ON NE TRIE PAS DE VECTEURS
!     ------------------------------------------------------------------
    integer(kind=8) :: iperm, i, iordre, j
    real(kind=8) :: rperm, eps
    complex(kind=8) :: cperm
!
!
!     --- TRI PAR ORDRE CROISSANT ---
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    eps = 1.d-7
    if (iordre .eq. 0) then
!
        do i = 1, nbpro
            iperm = i
            if (type .eq. 0) then
                rperm = dble(valpro(i))
                do j = i+1, nbpro
                    if (dble(valpro(j)) .le. rperm) then
                        iperm = j
                        rperm = dble(valpro(iperm))
                    end if
                end do
            else if (type .eq. 1) then
                rperm = abs(valpro(i))
                do j = i+1, nbpro
                    if (abs(valpro(j)) .lt. (rperm*(1.d0-eps))) then
                        iperm = j
                        rperm = abs(valpro(iperm))
                    end if
                    if ((abs(valpro(j))-rperm) .le. (eps*rperm)) then
                        if ((dble(valpro(j)*valpro(iperm)) .ge. 0.d0) .and. &
                            (abs(valpro(j)) .lt. rperm)) then
                            iperm = j
                            rperm = abs(valpro(iperm))
                        end if
                        if ((dble(valpro(j)*valpro(iperm)) .lt. 0.d0) .and. &
                            (dble(valpro(j)) .lt. 0.d0)) then
                            iperm = j
                            rperm = abs(valpro(iperm))
                        end if
                    end if
                end do
            end if
!
            if (iperm .ne. i) then
                cperm = valpro(iperm)
                valpro(iperm) = valpro(i)
                valpro(i) = cperm
                if (neq .ge. nbpro) then
                    do j = 1, neq
                        cperm = vecpro(j, i)
                        vecpro(j, i) = vecpro(j, iperm)
                        vecpro(j, iperm) = cperm
                    end do
                end if
            end if
        end do
!
    else if (iordre .eq. 1) then
!
        do i = 1, nbpro
            iperm = i
            if (type .eq. 0) then
                rperm = dble(valpro(i))
                do j = i+1, nbpro
                    if (dble(valpro(j)) .ge. rperm) then
                        iperm = j
                        rperm = dble(valpro(iperm))
                    end if
                end do
            else if (type .eq. 1) then
                rperm = abs(valpro(i))
                do j = i+1, nbpro
                    if (abs(valpro(j)) .gt. (rperm*(1.d0+eps))) then
                        iperm = j
                        rperm = abs(valpro(iperm))
                    end if
                end do
                if ((abs(valpro(j))-rperm) .le. (eps*rperm)) then
                    if ((dble(valpro(j)*valpro(iperm)) .ge. 0.d0) .and. &
                        (abs(valpro(j)) .gt. rperm)) then
                        iperm = j
                        rperm = abs(valpro(iperm))
                    end if
                    if ((dble(valpro(j)*valpro(iperm)) .lt. 0.d0) .and. &
                        (dble(valpro(j)) .lt. 0.d0)) then
                        iperm = j
                        rperm = abs(valpro(iperm))
                    end if
                end if
            end if
!
            if (iperm .ne. i) then
                cperm = valpro(iperm)
                valpro(iperm) = valpro(i)
                valpro(i) = cperm
                if (neq .ge. nbpro) then
                    do j = 1, neq
                        cperm = vecpro(j, i)
                        vecpro(j, i) = vecpro(j, iperm)
                        vecpro(j, iperm) = cperm
                    end do
                end if
            end if
        end do
!
    end if
!
end subroutine
