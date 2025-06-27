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
subroutine fsurf(option, nomte, xi, nb1, vecl, &
                 vectpt)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/fointe.h"
#include "asterfort/forsrg.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/r8inir.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vectci.h"
#include "asterfort/vexpan.h"
    character(len=16) :: option, nomte
!
!
    integer(kind=8) :: nb1
    real(kind=8) :: rnormc, f1, pr
    real(kind=8) :: xi(3, *), vectpt(9, 3, 3)
!     REAL*8       VECTC(3),VECPTX(3,3)
    real(kind=8) :: vecl(51), vecl1(42)
    real(kind=8) :: chgsrg(6, 8), chgsrl(6), chg(6)
    real(kind=8) :: kijkm1(40, 2), pgl(3, 3)
    aster_logical :: global, locapr
    real(kind=8) :: valpar(4)
    character(len=8) :: nompar(4)
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, ier, intsn, ip, itemps, j
    integer(kind=8) :: jp, jpres, k, lzi, lzr, nb2, npgsn
!
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
!
! --- CAS DES CHARGEMENTS DE FORME REEL
    if (option .eq. 'CHAR_MECA_FRCO3D') then
        call jevech('PFRCO3D', 'L', jpres)
        global = abs(zr(jpres-1+7)) .lt. 1.d-3
        locapr = abs(zr(jpres-1+7)-3.d0) .lt. 1.d-3
        if (global) then
            do j = 1, nb1
                do i = 1, 6
                    chgsrg(i, j) = zr(jpres-1+7*(j-1)+i)
                end do
            end do
        else
            do j = 1, nb1
                chgsrl(:) = 0.d0
                if (locapr) then
                    chgsrl(3) = -zr(jpres-1+7*(j-1)+3)
                else
                    do i = 1, 5
                        chgsrl(i) = zr(jpres-1+7*(j-1)+i)
                    end do
                end if
!
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
        end if
!
!
! --- CAS DES CHARGEMENTS DE FORME FONCTION
    else if (option .eq. 'CHAR_MECA_FFCO3D') then
        call jevech('PFFCO3D', 'L', jpres)
        call jevech('PINSTR', 'L', itemps)
        valpar(4) = zr(itemps)
        nompar(4) = 'INST'
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        global = zk8(jpres+6) .eq. 'GLOBAL'
        locapr = zk8(jpres+6) .eq. 'LOCAL_PR'
!
        if (global) then
! --        LECTURE DES INTERPOLATIONS DE FX, FY, FZ, MX, MY, MZ
            do j = 1, nb1
                valpar(1) = xi(1, j)
                valpar(2) = xi(2, j)
                valpar(3) = xi(3, j)
                call fointe('FM', zk8(jpres), 4, nompar, valpar, &
                            chgsrg(1, j), ier)
                call fointe('FM', zk8(jpres+1), 4, nompar, valpar, &
                            chgsrg(2, j), ier)
                call fointe('FM', zk8(jpres+2), 4, nompar, valpar, &
                            chgsrg(3, j), ier)
                call fointe('FM', zk8(jpres+3), 4, nompar, valpar, &
                            chgsrg(4, j), ier)
                call fointe('FM', zk8(jpres+4), 4, nompar, valpar, &
                            chgsrg(5, j), ier)
                call fointe('FM', zk8(jpres+5), 4, nompar, valpar, &
                            chgsrg(6, j), ier)
            end do
!
        else if (locapr) then
! --        BASE LOCALE - CAS D UNE PRESSION
! --        LECTURE DES INTERPOLATIONS DE LA PRESSION PRES
            do j = 1, nb1
                valpar(1) = xi(1, j)
                valpar(2) = xi(2, j)
                valpar(3) = xi(3, j)
                call fointe('FM', zk8(jpres+2), 4, nompar, valpar, &
                            pr, ier)
                chgsrl(3) = -1*pr
                chgsrl(1) = 0.d0
                chgsrl(2) = 0.d0
                chgsrl(4) = 0.d0
                chgsrl(5) = 0.d0
                chgsrl(6) = 0.d0
! --           CHANGEMENT DE BASE LOCAL --> GLOBAL
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
        else
!
! --        BASE LOCALE - CAS DE F1, F2, F3, MF1, MF2
! --        LECTURE DES INTERPOLATIONS DE F1, F2, F3, MF1, MF2
            do j = 1, nb1
                valpar(1) = xi(1, j)
                valpar(2) = xi(2, j)
                valpar(3) = xi(3, j)
                call fointe('FM', zk8(jpres), 4, nompar, valpar, &
                            chgsrl(1), ier)
                call fointe('FM', zk8(jpres+1), 4, nompar, valpar, &
                            chgsrl(2), ier)
                call fointe('FM', zk8(jpres+2), 4, nompar, valpar, &
                            chgsrl(3), ier)
                call fointe('FM', zk8(jpres+3), 4, nompar, valpar, &
                            chgsrl(4), ier)
                call fointe('FM', zk8(jpres+4), 4, nompar, valpar, &
                            chgsrl(5), ier)
                chgsrl(6) = 0.d0
! --           CHANGEMENT DE BASE LOCAL --> GLOBAL
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
        end if
!
    end if
!
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
