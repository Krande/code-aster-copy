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
subroutine fpres(nomte, xi, nb1, vecl, vectpt)
    implicit none
#include "jeveux.h"
#include "asterfort/forsrg.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/r8inir.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vectci.h"
#include "asterfort/vexpan.h"
    integer(kind=8) :: nb1
    character(len=16) :: nomte
    real(kind=8) :: vecl(51), vectpt(9, 3, 3)
!
!
    real(kind=8) :: rnormc, f1, chg(6), kijkm1(40, 2), pgl(3, 3)
    real(kind=8) :: xi(3, *), vecl1(42), chgsrg(6, 8), chgsrl(6)
    integer(kind=8) :: lzi, nb2, npgsn, lzr, jpres, j, i, jp, ip, intsn, i1, k
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsn = zi(lzi-1+4)
!
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
    call r8inir(42, 0.d0, vecl1, 1)
!
! --- CAS DES CHARGEMENTS DE PRESSION : Z LOCAL
!
    call jevech('PPRESSR', 'L', jpres)
    do j = 1, nb1
        do i = 1, 6
            chgsrl(i) = 0.d0
        end do
!-----------------------------------------------------
!  LE SIGNE MOINS CORRESPOND A LA CONVENTION :
!      UNE PRESSION POSITIVE PROVOQUE UN GONFLEMENT
        chgsrl(3) = -zr(jpres-1+j)
!-----------------------------------------------------
        do jp = 1, 3
            do ip = 1, 3
                pgl(jp, ip) = vectpt(j, jp, ip)
            end do
        end do
        call utpvlg(1, 6, pgl, chgsrl, chg)
        do i = 1, 6
            chgsrg(i, j) = chg(i)
        end do
    end do
!
    do intsn = 1, npgsn
        call vectci(intsn, nb1, xi, zr(lzr), rnormc)
!
        call forsrg(intsn, nb1, nb2, zr(lzr), chgsrg, &
                    rnormc, vectpt, vecl1)
    end do
!
!     RESTITUTION DE KIJKM1 POUR CONDENSER LES FORCES
!     ATTENTION LA ROUTINE N'EST PAS UTILISEE DANS LE CAS DES
!     EFFORTS SUIVANTS (MOMENTS SURFACIQUES)
!
    i1 = 5*nb1
    do j = 1, 2
        do i = 1, i1
            k = (j-1)*i1+i
            kijkm1(i, j) = zr(lzr-1+1000+k)
        end do
    end do
!
    do i = 1, i1
        f1 = 0.d0
        do k = 1, 2
            f1 = f1+kijkm1(i, k)*vecl1(i1+k)
        end do
        vecl1(i) = vecl1(i)-f1
    end do
!
!     EXPANSION DU VECTEUR VECL1 : DUE A L'AJOUT DE LA ROTATION FICTIVE
!
    call vexpan(nb1, vecl1, vecl)
    do i = 1, 3
        vecl(6*nb1+i) = 0.d0
    end do
!
end subroutine
