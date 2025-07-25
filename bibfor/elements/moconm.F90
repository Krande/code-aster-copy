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
subroutine moconm(dir, sigb, siga, hh, nlit, &
                  om, rr, nufsup, nufinf, nufsd1, &
                  nufid1, nufsd2, nufid2, prec)
    implicit none
! person_in_charge: sebastien.fayolle at edf.fr
!
#include "jeveux.h"
#include "asterfort/lsqpol.h"
#include "asterfort/wkvect.h"
    character(len=8) :: nufsup, nufinf, nufsd1, nufid1, nufsd2, nufid2
    character(len=1) :: dir
    integer(kind=8) :: nlit
    real(kind=8) :: sigb, siga(nlit), hh, om(nlit), rr(nlit), prec, e1, sigma
    integer(kind=8) :: ptmax, ordlu
    parameter(ptmax=50)
    parameter(ordlu=2)
    real(kind=8) :: nn(nlit*ptmax), mm(nlit*ptmax), eta, rhol(nlit+2)
    real(kind=8) :: xi(nlit+2), omm(nlit+2), nn0, mm0, poly(ordlu+1), xx
    integer(kind=8) :: i, j, k, ii, ilit, deb, npt, tri(nlit+2)
    integer(kind=8) :: ordok, jvale, jfon, jprol
!
! --- TRI SELON LA POSITION DANS L'EPAISSEUR
!    POUR AVOIR: RR(TRI(I))<=RR(TRI(I+1))
!    --> DU PLUS PETIT AU PLUS GRAND
!
    do i = 1, nlit
        rhol(i+1) = rr(i)
        omm(i+1) = om(i)
    end do
    rhol(1) = -1.0d0
    rhol(nlit+2) = 1.0d0
    omm(1) = hh
    omm(nlit+2) = hh
!
    do i = 1, nlit+2
        tri(i) = i
    end do
!
    if (nlit .gt. 1) then
        do j = 1, nlit-1
            do i = 2, nlit+1-j
                if (rr(tri(i)-1) .gt. rr(tri(i+1)-1)) then
                    ii = tri(i)
                    tri(i) = tri(i+1)
                    tri(i+1) = ii
                end if
            end do
        end do
    end if
!
! --- POSITIVE BENDING
!
! LE CAS DES POINTS SUPERPOSES EST TRAITE: OM=0 OU RR(I)=RR(I-1)
!    (PAR EX LINER SI RR=-1)
!
    do i = 1, nlit+2
        xi(i) = -1.0d0
    end do
    ii = 0
    do ilit = 0, nlit
        i = nlit-ilit+2
        xi(tri(i)) = 1.0d0
        nn0 = 0.d0
        mm0 = 0.d0
        do j = 1, nlit
            nn0 = nn0+xi(j+1)*om(j)*siga(j)
            mm0 = mm0+xi(j+1)*om(j)*siga(j)*rhol(j+1)*hh/2.0d0
        end do
        if (omm(tri(i)) .lt. 1.d-8*omm(1)) then
            deb = 1
        else
            deb = 0
        end if
        npt = int(abs(rhol(tri(i-1))-rhol(tri(i)))/2.d0*ptmax)-1
        npt = max(npt, 0)
        do k = deb, npt
            if (npt .eq. 0) then
                eta = rhol(tri(i))
            else
                eta = rhol(tri(i))+k*(rhol(tri(i-1))-rhol(tri(i)))/npt
            end if
            ii = ii+1
            nn(ii) = nn0+sigb*hh*(1+eta)/2.0d0-prec
            mm(ii) = mm0-sigb*hh*hh*(1-eta*eta)/8.0d0
        end do
    end do
!
! --- AJOUT DE LA FONCTION
    e1 = 0.d0
    npt = ii
    call lsqpol(ordlu, e1, npt, nn, mm, &
                ordok, poly, sigma)
!
!     --- CREATION ET REMPLISSAGE DE L'OBJET NUFSUP.VALE ---
    call wkvect(nufsup//'           .VALE', 'G V R', 2*npt, jvale)
    jfon = jvale+npt
    do i = 0, npt-1
        xx = nn(1)+(nn(npt)-nn(1))*i/(npt-1)
        zr(jvale+i) = xx
        zr(jfon+i) = 0.0d0
        do j = 0, ordok
            zr(jfon+i) = zr(jfon+i)+poly(j+1)*(xx**j)
        end do
    end do
!
    call wkvect(nufsd1//'           .VALE', 'G V R', 2*npt, jvale)
    jfon = jvale+npt
    do i = 0, npt-1
        xx = nn(1)+(nn(npt)-nn(1))*i/(npt-1)
        zr(jvale+i) = xx
        zr(jfon+i) = 0.0d0
        do j = 1, ordok
            zr(jfon+i) = zr(jfon+i)+j*poly(j+1)*(xx**(j-1))
        end do
    end do
!
    call wkvect(nufsd2//'           .VALE', 'G V R', 2*npt, jvale)
    jfon = jvale+npt
    do i = 0, npt-1
        xx = nn(1)+(nn(npt)-nn(1))*i/(npt-1)
        zr(jvale+i) = xx
        zr(jfon+i) = 0.0d0
        do j = 2, ordok-1
            zr(jfon+i) = zr(jfon+i)+j*(j-1)*poly(j+1)*(xx**(j-2))
        end do
    end do
!
!     --- CREATION ET REMPLISSAGE DE L'OBJET NUFSUP.PROL ---
    call wkvect(nufsup//'           .PROL', 'G V K24', 6, jprol)
    zk24(jprol) = 'FONCTION                '
    zk24(jprol+1) = 'LIN LIN                 '
    zk24(jprol+2) = 'X                       '
    zk24(jprol+3) = 'FME'//dir//'1                   '
    zk24(jprol+4) = 'LL                      '
!
    call wkvect(nufsd1//'           .PROL', 'G V K24', 6, jprol)
    zk24(jprol) = 'FONCTION                '
    zk24(jprol+1) = 'LIN LIN                 '
    zk24(jprol+2) = 'X                       '
    zk24(jprol+3) = 'DFME'//dir//'1                  '
    zk24(jprol+4) = 'CC                      '
!
    call wkvect(nufsd2//'           .PROL', 'G V K24', 6, jprol)
    zk24(jprol) = 'FONCTION                '
    zk24(jprol+1) = 'LIN LIN                 '
    zk24(jprol+2) = 'X                       '
    zk24(jprol+3) = 'DDFME'//dir//'1                 '
    zk24(jprol+4) = 'CC                      '
!
!--- NEGATIVE BENDING
!
    ii = 0
    do i = 1, nlit+1
        nn0 = 0.d0
        mm0 = 0.d0
        do j = 1, nlit
            nn0 = nn0-xi(j+1)*om(j)*siga(j)
            mm0 = mm0-xi(j+1)*om(j)*siga(j)*rhol(j+1)*hh/2.0d0
        end do
        if (omm(tri(i)) .lt. 1.d-8*omm(1)) then
            deb = 1
        else
            deb = 0
        end if
        npt = int(abs(rhol(tri(i+1))-rhol(tri(i)))/2.d0*ptmax)-1
        npt = max(npt, 0)
        do k = deb, npt
            if (npt .eq. 0) then
                eta = rhol(tri(i))
            else
                eta = rhol(tri(i))+k*(rhol(tri(i+1))-rhol(tri(i)))/npt
            end if
            ii = ii+1
            nn(ii) = nn0+sigb*hh*(1-eta)/2.0d0-prec
            mm(ii) = mm0+sigb*hh*hh*(1-eta*eta)/8.0d0
        end do
        xi(tri(i+1)) = -1
    end do
!
!--- AJOUT DE LA FONCTION
    e1 = 0.d0
    npt = ii
    call lsqpol(ordlu, e1, npt, nn, mm, &
                ordok, poly, sigma)
!
!     --- CREATION ET REMPLISSAGE DE L'OBJET NUFINF.VALE ---
    call wkvect(nufinf//'           .VALE', 'G V R', 2*npt, jvale)
    jfon = jvale+npt
    do i = 0, npt-1
        xx = nn(1)+(nn(npt)-nn(1))*i/(npt-1)
        zr(jvale+i) = xx
        zr(jfon+i) = 0.0d0
        do j = 0, ordok
            zr(jfon+i) = zr(jfon+i)+poly(j+1)*(xx**j)
        end do
    end do
!
    call wkvect(nufid1//'           .VALE', 'G V R', 2*npt, jvale)
    jfon = jvale+npt
    do i = 0, npt-1
        xx = nn(1)+(nn(npt)-nn(1))*i/(npt-1)
        zr(jvale+i) = xx
        zr(jfon+i) = 0.0d0
        do j = 1, ordok
            zr(jfon+i) = zr(jfon+i)+j*poly(j+1)*(xx**(j-1))
        end do
    end do
!
    call wkvect(nufid2//'           .VALE', 'G V R', 2*npt, jvale)
    jfon = jvale+npt
    do i = 0, npt-1
        xx = nn(1)+(nn(npt)-nn(1))*i/(npt-1)
        zr(jvale+i) = xx
        zr(jfon+i) = 0.0d0
        do j = 2, ordok-1
            zr(jfon+i) = zr(jfon+i)+j*(j-1)*poly(j+1)*(xx**(j-2))
        end do
    end do
!
!     --- CREATION ET REMPLISSAGE DE L'OBJET NUFSUP.PROL ---
    call wkvect(nufinf//'           .PROL', 'G V K24', 6, jprol)
    zk24(jprol) = 'FONCTION                '
    zk24(jprol+1) = 'LIN LIN                 '
    zk24(jprol+2) = 'X                       '
    zk24(jprol+3) = 'FME'//dir//'2                   '
    zk24(jprol+4) = 'LL                      '
!
    call wkvect(nufid1//'           .PROL', 'G V K24', 6, jprol)
    zk24(jprol) = 'FONCTION                '
    zk24(jprol+1) = 'LIN LIN                 '
    zk24(jprol+2) = 'X                       '
    zk24(jprol+3) = 'DFME'//dir//'2                  '
    zk24(jprol+4) = 'CC                      '
!
    call wkvect(nufid2//'           .PROL', 'G V K24', 6, jprol)
    zk24(jprol) = 'FONCTION                '
    zk24(jprol+1) = 'LIN LIN                 '
    zk24(jprol+2) = 'X                       '
    zk24(jprol+3) = 'DDFME'//dir//'2                 '
    zk24(jprol+4) = 'CC                      '
!
end subroutine
