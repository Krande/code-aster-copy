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
subroutine dfmdf(dim, f, dsidep)
!
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/diago2.h"
#include "asterfort/diago3.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: dim
    real(kind=8) :: f(dim), dsidep(dim, dim)
!
! ----------------------------------------------------------------------
!    CALCUL DE LA DERIVEE DE LA PARTIE NEGATIVE D UN TENSEUR SYMETRIQUE
!    PAR RAPPORT A CE TENSEUR
!
!    DEFINITION:
!    F SYMETRIQUE DONC DIAGONALISABLE:
!                      L1 0  0
!             F=  MT.  0  L2 0   .M
!                      0  0  L3
!
!    LA PARTIE NEGATIVE S ECRIT (H FONCTION DE HEAVISIDE):
!                        H(-L1)  0       0
!             FM=  MT.   0       H(-L2)  0       .M
!                        0       0       H(-L3)
!
!    ON CHERCHE DFM/DF
!
!    TENSEURS ENTRES SOUS FORME VECTORIELLE:
!    MATRICE (3x3)---> VECTEUR(6) ORDRE (XX YY ZZ XY XZ YZ)
!    MATRICE (2x2)---> VECTEUR(3) ORDRE (XX YY XY)
!    MATRICE (1x1)---> VECTEUR(1)
!
!    IN DIM   : DIMENSION DU VECTEUR (6 POUR 3D, 3 POUR 2D, 1 POUR 1D)
!    IN F     : VECTEUR (=TENSEUR)
!    OUT      : DFM/DF
!
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, l, m, n, t(3, 3), t2(2, 2), ordre(2)
    real(kind=8) :: rtemp, rtemp2, dspdep(3, 3), epsp(3), vecp(3, 3), rigmin
    real(kind=8) :: sigp(3)
    real(kind=8) :: dspdeb(2, 2), epsp2(2), vecp2(2, 2), sigp2(2)
    real(kind=8) :: rtemp3, rtemp4, rac2
    aster_logical :: mtg, test
    parameter(rigmin=1.d-6)
!
    t(1, 1) = 1
    t(1, 2) = 4
    t(1, 3) = 5
    t(2, 1) = 4
    t(2, 2) = 2
    t(2, 3) = 6
    t(3, 1) = 5
    t(3, 2) = 6
    t(3, 3) = 3
!
    t2(1, 1) = 1
    t2(1, 2) = 3
    t2(2, 1) = 3
    t2(2, 2) = 2
!
    rac2 = sqrt(2.d0)
!
    call r8inir(dim*dim, 0.d0, dsidep, 1)
!
    if (dim .eq. 6) then
!
        call diago3(f, vecp, epsp)
        mtg = .true.
        rtemp = abs(epsp(1))
! TEST = .TRUE. PAR DEFAUT : SIGNIFIE QUE LA MATRICE EN ENTREE
! EST DEJA DIAGONALE. LA PROCEDURE UTILISEE EST DE COMPARER CHAQUE
! VALEUR PROPRE AU TERME DIAGONAL CORRESPONDANT DE LA MATRICE EN
! ENTREE. SI AU MOINS L'UNE DES VALEURS PROPRES DIFFERE DE LA VALEUR
! DE LA DIAGONALE, ALORS TEST=.FALSE. ET ON NE CORRIGE PAS LES VECTEURS
! PROPRES
        test = .true.
        do i = 1, 3
            if (epsp(i) .ne. 0.d0) then
                if (abs(f(i)-epsp(i))/abs(epsp(i)) .gt. 1.d-12) then
                    test = .false.
                end if
            else
                if (abs(f(i)-epsp(i)) .gt. 1.d-12) then
                    test = .false.
                end if
            end if
        end do
! SI TEST = .TRUE. : CORRECTION DES VECTEURS PROPRES
        if (test) then
            do i = 1, 3
                do j = 1, 3
                    vecp(i, j) = 0.d0
                end do
                vecp(i, i) = 1.d0
            end do
        end if
! FIN CORRECTION DES VALEURS PROPRES
        if (abs(epsp(2)) .gt. rtemp) rtemp = abs(epsp(2))
        if (abs(epsp(3)) .gt. rtemp) rtemp = abs(epsp(3))
        do i = 1, 2
            do j = (i+1), 3
                if (abs(epsp(i)-epsp(j)) .lt. 1d-12) then
                    epsp(i) = epsp(i)+rigmin*rtemp
                    epsp(j) = epsp(j)-rigmin*rtemp
                end if
            end do
        end do
        do i = 1, 3
            if (epsp(i) .lt. 0.d0) then
                sigp(i) = epsp(i)
            else
                sigp(i) = 0.d0
            end if
        end do
        mtg = .true.
        call r8inir(9, 0.d0, dspdep, 1)
        call r8inir(36, 0.d0, dsidep, 1)
        do k = 1, 3
            if (epsp(k) .lt. 0.d0) then
                dspdep(k, k) = 1.d0
            else
                dspdep(k, k) = 0.d0
            end if
        end do
        do i = 1, 3
            do j = i, 3
                do k = 1, 3
                    do l = 1, 3
                        do m = 1, 3
                            do n = 1, 3
                                if (i .eq. j) then
                                    rtemp3 = 1.d0
                                else
                                    rtemp3 = rac2
                                end if
                                if (k .eq. l) then
                                    rtemp4 = 1.d0
                                else
                                    rtemp4 = 1.d0/rac2
                                end if
                                dsidep(t(i, j), t(k, l)) = dsidep(t(i, j), t( &
                                                          k, l))+vecp(k, m)*vecp(i, n)*vecp(j, n)* &
                                                           vecp(l, m)*dspdep(n, m)*rtemp3*rtemp4
                                rtemp = abs(epsp(m)-epsp(n))
                                if ((m .ne. n)) then
                                    if ((rtemp .gt. 1.d-12)) then
                                        rtemp2 = (vecp(k, m)*vecp(l, n))/( &
                                                 epsp(n)-epsp(m))
                                        rtemp2 = rtemp2*vecp(i, m)*vecp( &
                                                 j, n)*sigp(n)*rtemp3*rtemp4
                                        dsidep(t(i, j), t(k, l)) = dsidep( &
                                                                   t(i, j), t(k, l))+rtemp2
                                        rtemp2 = (vecp(k, n)*vecp(l, m))/( &
                                                 epsp(n)-epsp(m))
                                        rtemp2 = rtemp2*vecp(j, m)*vecp( &
                                                 i, n)*sigp(n)*rtemp3*rtemp4
                                        dsidep(t(i, j), t(k, l)) = dsidep( &
                                                                   t(i, j), t(k, l))+rtemp2
                                    else
                                        mtg = .false.
                                    end if
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do
!
        if (.not. mtg) then
            do k = 1, 6
                do l = 1, 6
                    dsidep(k, l) = 0.d0
                end do
            end do
            dsidep(1, 1) = 1.d0
            dsidep(2, 2) = 1.d0
            dsidep(3, 3) = 1.d0
            dsidep(4, 4) = 1.d0
            dsidep(5, 5) = 1.d0
            dsidep(6, 6) = 1.d0
        end if
!
    else if (dim .eq. 3) then
!
        call diago2(f, vecp2, epsp2)
        mtg = .true.
        rtemp = abs(epsp2(1))
! VERIFICATION DE L'ORDRE DES VALEURS PROPRES, UTILE LORSQUE
! LA MATRICE EN ENTREE EST DEJA DIAGONALE AFIN D'IMPOSER LES
! VECTEURS PROPRES A (1 0 , 0 1) OU (0 -1 , 1 0) SELON L'ORDRE
! DES VALEURS PROPRES CALCULEES PAR DIAGO2
        ordre(1) = 1
        ordre(2) = 2
        if (f(2) .gt. f(1)) then
            ordre(1) = 2
            ordre(2) = 1
        end if
! TEST = .TRUE. PAR DEFAUT : SIGNIFIE QUE LA MATRICE EN ENTREE
! EST DEJA DIAGONALE. LA PROCEDURE UTILISEE EST DE COMPARER CHAQUE
! VALEUR PROPRE AU TERME DIAGONAL CORRESPONDANT DE LA MATRICE EN
! ENTREE. SI AU MOINS L'UNE DES VALEURS PROPRES DIFFERE DE LA VALEUR
! DE LA DIAGONALE, ALORS TEST=.FALSE. ET ON NE CORRIGE PAS LES VECTEURS
! PROPRES
        test = .true.
        do i = 1, 2
            if (epsp2(i) .ne. 0.d0) then
                if (abs(f(ordre(i))-epsp2(i))/abs(epsp2(i)) .gt. 1.d-12) then
                    test = .false.
                end if
            else
                if (abs(f(ordre(i))-epsp2(i)) .gt. 1.d-12) then
                    test = .false.
                end if
            end if
        end do
! SI TEST = .TRUE. : CORRECTION DES VECTEURS PROPRES
        if (test) then
            if (ordre(1) .eq. 1) then
                do i = 1, 2
                    do j = 1, 2
                        vecp2(i, j) = 0.d0
                    end do
                    vecp2(i, i) = 1.d0
                end do
            else
                vecp2(1, 1) = 0.d0
                vecp2(1, 2) = 1.d0
                vecp2(2, 1) = -1.d0
                vecp2(2, 2) = 0.d0
            end if
        end if
! FIN CORRECTION DES VALEURS PROPRES
        if (abs(epsp2(2)) .gt. rtemp) rtemp = abs(epsp2(2))
        if (abs(epsp2(1)-epsp2(2)) .lt. 1d-12) then
            epsp2(1) = epsp2(1)+rigmin*rtemp
            epsp2(2) = epsp2(2)-rigmin*rtemp
        end if
!
        do i = 1, 2
            if (epsp2(i) .lt. 0.d0) then
                sigp2(i) = epsp2(i)
            else
                sigp2(i) = 0.d0
            end if
        end do
        mtg = .true.
        call r8inir(4, 0.d0, dspdeb, 1)
!
        do k = 1, 2
            if (epsp2(k) .lt. 0.d0) then
                dspdeb(k, k) = 1.d0
            else
                dspdeb(k, k) = 0.d0
            end if
        end do
        do i = 1, 2
            do j = i, 2
                do k = 1, 2
                    do l = 1, 2
                        do m = 1, 2
                            do n = 1, 2
                                if (i .eq. j) then
                                    rtemp3 = 1.d0
                                else
                                    rtemp3 = rac2
                                end if
                                if (k .eq. l) then
                                    rtemp4 = 1.d0
                                else
                                    rtemp4 = 1.d0/rac2
                                end if
                                dsidep(t2(i, j), t2(k, l)) = dsidep(t2(i, j) &
                                                              , t2(k, l))+vecp2(k, m)*vecp2(i, n)* &
                                                             vecp2(j, n)*vecp2(l, m)*dspdeb(n, m)* &
                                                             rtemp3*rtemp4
                                rtemp = abs(epsp2(m)-epsp2(n))
                                if ((m .ne. n)) then
                                    if ((rtemp .gt. 1.d-12)) then
                                        rtemp2 = (vecp2(k, m)*vecp2(l, n)) &
                                                 /(epsp2(n)-epsp2(m))
                                        rtemp2 = rtemp2*vecp2(i, m)* &
                                                 vecp2(j, n)*sigp2(n)*rtemp3* &
                                                 rtemp4
                                        dsidep(t2(i, j), t2(k, l)) = &
                                            dsidep(t2(i, j), t2(k, l))+ &
                                            rtemp2
                                        rtemp2 = (vecp2(k, n)*vecp2(l, m)) &
                                                 /(epsp2(n)-epsp2(m))
                                        rtemp2 = rtemp2*vecp2(j, m)* &
                                                 vecp2(i, n)*sigp2(n)*rtemp3* &
                                                 rtemp4
                                        dsidep(t2(i, j), t2(k, l)) = &
                                            dsidep(t2(i, j), t2(k, l))+ &
                                            rtemp2
                                    else
                                        mtg = .false.
                                    end if
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do
        if (.not. mtg) then
            dsidep(1, 1) = 1.d0
            dsidep(2, 2) = 1.d0
            dsidep(3, 3) = 1.d0
        end if
!
    else if (dim .eq. 1) then
!
        if (f(1) .lt. 0.d0) then
            dsidep(1, 1) = 1.d0
        else
            dsidep(1, 1) = 0.d0
        end if
!
    end if
!
end subroutine
