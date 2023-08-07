! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1501
!
subroutine elraga(elrefz, fapz, ndim, nbpg, coopg, poipg)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/elraca.h"
!
    character(len=*), intent(in) :: elrefz, fapz
    integer, intent(out) :: nbpg, ndim
    real(kind=8), intent(out) :: coopg(*), poipg(*)
!
! --------------------------------------------------------------------------------------------------
!
! Finite elements management
!
! Get parameters of integration scheme
!
! --------------------------------------------------------------------------------------------------
!
! In  elrefe           : name of geometric support for finite element
! In  fapg             : name of Gauss integration scheme
! Out ndim             : topological dimension (0/1/2/3)
! Out nbpg             : number of points of integration schemes
! Out coopg            : coordinataes of Gauss points
! Out poipg            : weight of Gauss points
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elrefa, fapg, nofpg(MT_NBFAMX)
    integer :: i, npar, npi, ix, iy, iz, npx, npyz
    integer :: nno, nnos, nbfpg, nbpg1(MT_NBFAMX), iNode, ifam
    real(kind=8) :: xpg(MT_NBPGMX), ypg(MT_NBPGMX), zpg(MT_NBPGMX), hpg(MT_NBPGMX)
    real(kind=8) :: h(4), a(4)
    real(kind=8) :: aty(7), ht(7), atz(7)
    real(kind=8) :: lobWeight(7), lobCoor(7)
    real(kind=8) :: a1, a2, b1, b2, c1, c2, d1, e1
    real(kind=8) :: h1, h2, h3, h5
    real(kind=8) :: p1, p2, p3, p4, p5
    real(kind=8) :: xno(3*MT_NNOMAX), vol
    real(kind=8), parameter :: zero = 0.d0, undemi = 0.5d0
    real(kind=8), parameter :: un = 1.d0, deux = 2.d0

    real(kind=8), parameter :: rac5 = sqrt(5.d0), rac15 = sqrt(15.d0), rac30 = sqrt(30.d0)
    real(kind=8), parameter :: rac_1div3 = sqrt(1.d0/3.d0), rac_3div5 = sqrt(3.d0/5.d0)
    real(kind=8), parameter :: rac_3div7 = sqrt(3.d0/7.d0)

    real(kind=8), parameter :: gauss4p12 = sqrt(3.d0/7.d0-2.d0/7.d0*sqrt(6.d0/5.d0))
    real(kind=8), parameter :: gauss4p34 = sqrt(3.d0/7.d0+2.d0/7.d0*sqrt(6.d0/5.d0))

    real(kind=8), parameter :: lobatto7p35 = sqrt(5.d0/11.d0-2.d0/11.d0*sqrt(5.d0/3.d0))
    real(kind=8), parameter :: lobatto7p26 = sqrt(5.d0/11.d0+2.d0/11.d0*sqrt(5.d0/3.d0))

!
#define t(u) 2.0d0*(u) - 1.0d0
!
! --------------------------------------------------------------------------------------------------
!
    elrefa = elrefz
    fapg = fapz

! - Get list of integration schemes of geometric support
    call elraca(elrefa, &
                nbfpg, nofpg, nbpg1, &
                ndim, nno, nnos, &
                xno, vol)
    ASSERT((ndim .ge. 0) .and. (ndim .le. 3))

! - Get index for integration scheme
    ifam = indik8(nofpg, fapg, 1, nbfpg)
    ASSERT(ifam .gt. 0)

! - Get number of Gauss points
    nbpg = nbpg1(ifam)

! - For 'NOEU' scheme
    if (fapg .eq. 'NOEU') then
        ASSERT(nbpg .eq. nno)
        do iNode = 1, nno
            hpg(iNode) = vol/nno
            if (ndim .ge. 1) xpg(iNode) = xno(ndim*(iNode-1)+1)
            if (ndim .ge. 2) ypg(iNode) = xno(ndim*(iNode-1)+2)
            if (ndim .eq. 3) zpg(iNode) = xno(ndim*(iNode-1)+3)
        end do
        goto 170
    end if

! - For 'NOEU_S' scheme
    if (fapg .eq. 'NOEU_S') then
        ASSERT(nbpg .eq. nnos)
        do iNode = 1, nnos
            hpg(iNode) = vol/nnos
            if (ndim .ge. 1) xpg(iNode) = xno(ndim*(iNode-1)+1)
            if (ndim .ge. 2) ypg(iNode) = xno(ndim*(iNode-1)+2)
            if (ndim .eq. 3) zpg(iNode) = xno(ndim*(iNode-1)+3)
        end do
        goto 170
    end if

! - For 'FPG1' scheme
    if (fapg .eq. 'FPG1') then
        ASSERT(nbpg .eq. 1)
        xpg(1) = zero
        if (ndim .ge. 1) xpg(1) = zero
        if (ndim .ge. 2) ypg(1) = zero
        if (ndim .eq. 3) zpg(1) = zero
        do iNode = 1, nno
            if (ndim .ge. 1) xpg(1) = xpg(1)+xno(ndim*(iNode-1)+1)
            if (ndim .ge. 2) ypg(1) = ypg(1)+xno(ndim*(iNode-1)+2)
            if (ndim .eq. 3) zpg(1) = zpg(1)+xno(ndim*(iNode-1)+3)
        end do
        if (ndim .ge. 1) xpg(1) = xpg(1)/nno
        if (ndim .ge. 2) ypg(1) = ypg(1)/nno
        if (ndim .eq. 3) zpg(1) = zpg(1)/nno
        hpg(1) = vol
        goto 170
    end if

! - For other schemes
    if (elrefa .eq. 'HE8' .or. elrefa .eq. 'H20' .or. elrefa .eq. 'H27') then
        npar = 0
        if (fapg .eq. 'FPG1') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 1 POINTS ( ORDRE 1 )
            xpg(1) = zero
            ypg(1) = zero
            zpg(1) = zero
            hpg(1) = 8.d0

        else if (fapg .eq. 'FPG8') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 2 POINTS DANS CHAQUE DIRECTION ( ORDRE 3 )
            npar = 2
            a(1) = -rac_1div3
            a(2) = -a(1)
            h(1) = un
            h(2) = un

        else if (fapg .eq. 'FPG27') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 3 POINTS DANS CHAQUE DIRECTION ( ORDRE 5 )
            npar = 3
            a(1) = -rac_3div5
            a(2) = zero
            a(3) = -a(1)
            h(1) = 5.d0/9.d0
            h(2) = 8.d0/9.d0
            h(3) = h(1)

        else if (fapg .eq. 'FPG64') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 4 POINTS DANS CHAQUE DIRECTION ( ORDRE 7 )
            npar = 4
            a(1) = -gauss4p12
            a(2) = -a(1)
            a(3) = -gauss4p34
            a(4) = -a(3)
            h(1) = (18.d0+rac30)/36.d0
            h(2) = h(1)
            h(3) = (18.d0-rac30)/36.d0
            h(4) = h(3)

        else if (fapg .eq. 'FPG8NOS') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 2 POINTS DANS CHAQUE DIRECTION ( ORDRE 3 )
            npar = 2
            a(1) = -rac_1div3
            a(2) = -a(1)
            h(1) = un
            h(2) = un
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+8) = vol/nnos
                xpg(iNode+8) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+8) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+8) = xno(ndim*(iNode-1)+3)
            end do

        else
            ASSERT(ASTER_FALSE)

        end if
        npi = 0
        do ix = 1, npar
            do iy = 1, npar
                do iz = 1, npar
                    npi = npi+1
                    xpg(npi) = a(ix)
                    ypg(npi) = a(iy)
                    zpg(npi) = a(iz)
                    hpg(npi) = h(ix)*h(iy)*h(iz)
                end do
            end do
        end do

    else if (elrefa .eq. 'HE9') then
        if (fapg .eq. 'LOB5') then
! --------- FORMULE DE QUADRATURE DE GAUSS-LOBATTO A 5 POINTS DANS
!           L EPAISSEUR AU CENTRE DE L'ELEMENT
            lobCoor(1) = -un
            lobCoor(2) = -rac_3div7
            lobCoor(3) = zero
            lobCoor(4) = -lobCoor(2)
            lobCoor(5) = -lobCoor(1)
            lobWeight(1) = 0.1d0
            lobWeight(2) = 49.d0/90.d0
            lobWeight(3) = 32.0/45.d0
            lobWeight(4) = lobWeight(2)
            lobWeight(5) = lobWeight(1)
            do iz = 1, 5
                xpg(iz) = zero
                ypg(iz) = zero
                zpg(iz) = lobCoor(iz)
                hpg(iz) = lobWeight(iz)*4.d0
            end do

        else if (fapg .eq. 'LOB7') then
! --------- FORMULE DE QUADRATURE DE GAUSS-LOBATTO A 7 POINTS DANS
!           L EPAISSEUR AU CENTRE DE L'ELEMENT
            lobCoor(1) = -un
            lobCoor(2) = -lobatto7p26
            lobCoor(3) = -lobatto7p35
            lobCoor(4) = zero
            lobCoor(5) = -lobCoor(3)
            lobCoor(6) = -lobCoor(2)
            lobCoor(7) = -lobCoor(1)
            lobWeight(1) = un/21.d0
            lobWeight(2) = (124.d0-7.d0*rac15)/350.d0
            lobWeight(3) = (124.d0+7.d0*rac15)/350.d0
            lobWeight(4) = 256.d0/525.d0
            lobWeight(5) = lobWeight(3)
            lobWeight(6) = lobWeight(2)
            lobWeight(7) = lobWeight(1)
            do iz = 1, 7
                xpg(iz) = zero
                ypg(iz) = zero
                zpg(iz) = lobCoor(iz)
                hpg(iz) = lobWeight(iz)*4.d0
            end do

        else if (fapg .eq. 'FPG8') then
! --------- FORMULE DE QUADRATURE DE GAUSS A 2 POINTS DANS CHAQUE DIRECTION ( ORDRE 3 )
            npar = 2
            a(1) = -rac_1div3
            a(2) = -a(1)
            h(1) = un
            h(2) = un
            npi = 0
            do ix = 1, npar
                do iy = 1, npar
                    do iz = 1, npar
                        npi = npi+1
                        xpg(npi) = a(ix)
                        ypg(npi) = a(iy)
                        zpg(npi) = a(iz)
                        hpg(npi) = h(ix)*h(iy)*h(iz)
                    end do
                end do
            end do

        else
            ASSERT(ASTER_FALSE)
        end if
!
    else if (elrefa .eq. 'PE6' .or. elrefa .eq. 'P15' .or. elrefa .eq. 'P18' .or. &
             elrefa .eq. 'P21') then
        if (fapg .eq. 'FPG6') then
! --------- FORMULE A 4 * 2 POINTS (CF TOUZOT PAGE 297) -> ORDRE 2
! --------- FORMULE DE GAUSS - 2 POINTS DE GAUSS  EN X (ORDRE 3)
            npx = 2
            a(1) = rac_1div3
            a(2) = -a(1)
            h(1) = un
            h(2) = un

! --------- FORMULE DE HAMMER - 3 POINTS DE HAMMER EN Y Z (ORDRE 2 EN Y Z)
            npyz = 3
            aty(1) = undemi
            aty(2) = zero
            aty(3) = undemi
            atz(1) = undemi
            atz(2) = undemi
            atz(3) = zero

            ht(1) = un/6.d0
            ht(2) = ht(1)
            ht(3) = ht(1)

        else if (fapg .eq. 'FPG6B') then
! --------- FORMULE A 4 * 2 POINTS (CF TOUZOT PAGE 297) -> ORDRE 2
! --------- FORMULE DE GAUSS - 2 POINTS DE GAUSS  EN X (ORDRE 3)
            npx = 2
            a(1) = rac_1div3
            a(2) = -a(1)

            h(1) = un
            h(2) = un

! --------- FORMULE DE HAMMER - 3 POINTS DE HAMMER EN Y Z (ORDRE 2 EN Y Z)
            npyz = 3
            aty(1) = un/6.d0
            aty(2) = deux/3.d0
            aty(3) = un/6.d0
            atz(1) = un/6.d0
            atz(2) = un/6.d0
            atz(3) = deux/3.d0

            ht(1) = un/6.d0
            ht(2) = ht(1)
            ht(3) = ht(1)

        else if (fapg .eq. 'FPG8') then
! --------- FORMULE A 4 * 2 POINTS (CF TOUZOT PAGE 297) -> ORDRE 3
! --------- FORMULE DE GAUSS - 2 POINTS DE GAUSS  EN X (ORDRE 3)
            npx = 2
            a(1) = rac_1div3
            a(2) = -a(1)
            h(1) = un
            h(2) = un

! --------- FORMULE DE HAMMER - 4 POINTS DE HAMMER EN Y Z (ORDRE 3 EN Y Z)
            npyz = 4
            aty(1) = un/3.d0
            aty(2) = 0.6d0
            aty(3) = 0.2d0
            aty(4) = 0.2d0
            atz(1) = un/3.d0
            atz(2) = 0.2d0
            atz(3) = 0.6d0
            atz(4) = 0.2d0
            ht(1) = -27.d0/96.d0
            ht(2) = 25.d0/96.d0
            ht(3) = ht(2)
            ht(4) = ht(2)

        else if (fapg .eq. 'FPG21') then
! --------- FORMULE A 7 * 3 POINTS :   (CF TOUZOT PAGE 298)
! --------- FORMULE DE GAUSS - 3 POINTS DE GAUSS EN X (ORDRE 5)
            npx = 3
            a(1) = -rac_3div5
            a(2) = zero
            a(3) = -a(1)
            h(1) = 5.d0/9.d0
            h(2) = 8.d0/9.d0
            h(3) = h(1)

! --------- FORMULE DE HAMMER - 7 POINTS DE HAMMER EN Y Z (ORDRE 5 EN Y Z)
            npyz = 7
            aty(1) = un/3.d0
            atz(1) = un/3.d0
            aty(2) = (6.d0+rac15)/21.d0
            atz(2) = aty(2)
            aty(3) = un-deux*aty(2)
            atz(3) = aty(2)
            aty(4) = aty(2)
            atz(4) = un-deux*aty(2)
            aty(5) = (6.d0-rac15)/21.d0
            atz(5) = aty(5)
            aty(6) = un-deux*aty(5)
            atz(6) = aty(5)
            aty(7) = aty(5)
            atz(7) = un-deux*aty(5)
            ht(1) = 9.d0/80.d0
            ht(2) = (155.d0+rac15)/2400.d0
            ht(3) = ht(2)
            ht(4) = ht(2)
            ht(5) = (155.d0-rac15)/2400.d0
            ht(6) = ht(5)
            ht(7) = ht(5)

        else if (fapg .eq. 'FPG29') then
            ! ORDRE 6
            ! x = zeta
            ! y = (xi + 1) / 2
            ! z = (eta + 1) / 2
            ! h = weight / 4

            hpg(1) = 0.027191062410231d0
            hpg(2) = 0.027191062410231d0
            hpg(3) = 0.027191062410231d0
            hpg(4) = 0.040636041641220d0
            hpg(5) = 0.040636041641220d0
            hpg(6) = 0.040636041641220d0
            hpg(7) = 0.040636041641220d0
            hpg(8) = 0.040636041641220d0
            hpg(9) = 0.040636041641220d0
            hpg(10) = 0.050275140937507d0
            hpg(11) = 0.011774414962347d0
            hpg(12) = 0.011774414962347d0
            hpg(13) = 0.011774414962347d0
            hpg(14) = 0.041951149272741d0
            hpg(15) = 0.041951149272741d0
            hpg(16) = 0.041951149272741d0
            hpg(17) = 0.041951149272741d0
            hpg(18) = 0.041951149272741d0
            hpg(19) = 0.041951149272741d0
            hpg(20) = 0.050275140937507d0
            hpg(21) = 0.011774414962347d0
            hpg(22) = 0.011774414962347d0
            hpg(23) = 0.011774414962347d0
            hpg(24) = 0.041951149272741d0
            hpg(25) = 0.041951149272741d0
            hpg(26) = 0.041951149272741d0
            hpg(27) = 0.041951149272741d0
            hpg(28) = 0.041951149272741d0
            hpg(29) = 0.041951149272741d0

            xpg(1) = 0.000000000000000d0
            xpg(2) = 0.000000000000000d0
            xpg(3) = 0.000000000000000d0
            xpg(4) = 0.000000000000000d0
            xpg(5) = 0.000000000000000d0
            xpg(6) = 0.000000000000000d0
            xpg(7) = 0.000000000000000d0
            xpg(8) = 0.000000000000000d0
            xpg(9) = 0.000000000000000d0
            xpg(10) = 0.936241512371697d0
            xpg(11) = 0.948681147283254d0
            xpg(12) = 0.948681147283254d0
            xpg(13) = 0.948681147283254d0
            xpg(14) = 0.600638052820557d0
            xpg(15) = 0.600638052820557d0
            xpg(16) = 0.600638052820557d0
            xpg(17) = 0.600638052820557d0
            xpg(18) = 0.600638052820557d0
            xpg(19) = 0.600638052820557d0
            xpg(20) = -0.936241512371697d0
            xpg(21) = -0.948681147283254d0
            xpg(22) = -0.948681147283254d0
            xpg(23) = -0.948681147283254d0
            xpg(24) = -0.600638052820557d0
            xpg(25) = -0.600638052820557d0
            xpg(26) = -0.600638052820557d0
            xpg(27) = -0.600638052820557d0
            xpg(28) = -0.600638052820557d0
            xpg(29) = -0.600638052820557d0

            ypg(1) = 0.895512822481133d0
            ypg(2) = 0.052243588759434d0
            ypg(3) = 0.052243588759434d0
            ypg(4) = 0.198304865473555d0
            ypg(5) = 0.198304865473555d0
            ypg(6) = 0.270635256143164d0
            ypg(7) = 0.531059878383280d0
            ypg(8) = 0.531059878383280d0
            ypg(9) = 0.270635256143164d0
            ypg(10) = 0.333333333333333d0
            ypg(11) = 0.841699897299232d0
            ypg(12) = 0.079150051350384d0
            ypg(13) = 0.079150051350384d0
            ypg(14) = 0.054831294873304d0
            ypg(15) = 0.054831294873304d0
            ypg(16) = 0.308513201856883d0
            ypg(17) = 0.636655503269814d0
            ypg(18) = 0.636655503269814d0
            ypg(19) = 0.308513201856883d0
            ypg(20) = 0.333333333333333d0
            ypg(21) = 0.841699897299232d0
            ypg(22) = 0.079150051350384d0
            ypg(23) = 0.079150051350384d0
            ypg(24) = 0.054831294873304d0
            ypg(25) = 0.054831294873304d0
            ypg(26) = 0.308513201856883d0
            ypg(27) = 0.636655503269814d0
            ypg(28) = 0.636655503269814d0
            ypg(29) = 0.308513201856883d0

            zpg(1) = 0.052243588759434d0
            zpg(2) = 0.895512822481133d0
            zpg(3) = 0.052243588759434d0
            zpg(4) = 0.270635256143164d0
            zpg(5) = 0.531059878383280d0
            zpg(6) = 0.531059878383280d0
            zpg(7) = 0.270635256143164d0
            zpg(8) = 0.198304865473555d0
            zpg(9) = 0.198304865473555d0
            zpg(10) = 0.333333333333333d0
            zpg(11) = 0.079150051350384d0
            zpg(12) = 0.841699897299232d0
            zpg(13) = 0.079150051350384d0
            zpg(14) = 0.308513201856883d0
            zpg(15) = 0.636655503269814d0
            zpg(16) = 0.636655503269814d0
            zpg(17) = 0.308513201856883d0
            zpg(18) = 0.054831294873304d0
            zpg(19) = 0.054831294873304d0
            zpg(20) = 0.333333333333333d0
            zpg(21) = 0.079150051350384d0
            zpg(22) = 0.841699897299232d0
            zpg(23) = 0.079150051350384d0
            zpg(24) = 0.308513201856883d0
            zpg(25) = 0.636655503269814d0
            zpg(26) = 0.636655503269814d0
            zpg(27) = 0.308513201856883d0
            zpg(28) = 0.054831294873304d0
            zpg(29) = 0.054831294873304d0

        else if (fapg .eq. 'FPG6NOS') then
! --------- POUR LES POINTS DE GAUSS
            npx = 2
            npyz = 3
            a(1) = -rac_1div3
            a(2) = rac_1div3
            aty(1) = undemi
            aty(2) = zero
            aty(3) = undemi
            atz(1) = undemi
            atz(2) = undemi
            atz(3) = zero
            h(1) = un
            h(2) = un
            ht(1) = un/6.d0
            ht(2) = ht(1)
            ht(3) = ht(1)
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+6) = vol/nnos
                xpg(iNode+6) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+6) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+6) = xno(ndim*(iNode-1)+3)
            end do
        else
            ASSERT(ASTER_FALSE)

        end if

        if (fapg .ne. 'FPG29') then
            npi = 0
            do ix = 1, npx
                do iy = 1, npyz
                    npi = npi+1
                    xpg(npi) = a(ix)
                    ypg(npi) = aty(iy)
                    zpg(npi) = atz(iy)
                    hpg(npi) = h(ix)*ht(iy)
                end do
            end do
        end if

    else if (elrefa .eq. 'PE7') then
        if (fapg .eq. 'LOB5') then
! --------- FORMULE DE QUADRATURE DE GAUSS-LOBATTO A 5 POINTS DANS
!           L EPAISSEUR AU CENTRE DE L'ELEMENT
            lobCoor(1) = -un
            lobCoor(2) = -rac_3div7
            lobCoor(3) = zero
            lobCoor(4) = -lobCoor(2)
            lobCoor(5) = un
            lobWeight(1) = 0.1d0
            lobWeight(2) = 49.d0/90.d0
            lobWeight(3) = 32.0/45.d0
            lobWeight(4) = lobWeight(2)
            lobWeight(5) = lobWeight(1)
            do iz = 1, 5
                xpg(iz) = un/3.d0
                ypg(iz) = un/3.d0
                zpg(iz) = lobCoor(iz)
                hpg(iz) = lobWeight(iz)*undemi
            end do

        else if (fapg .eq. 'LOB7') then
! --------- FORMULE DE QUADRATURE DE GAUSS-LOBATTO A 7 POINTS DANS
!           L EPAISSEUR AU CENTRE DE L'ELEMENT
            lobCoor(1) = -un
            lobCoor(2) = -lobatto7p26
            lobCoor(3) = -lobatto7p35
            lobCoor(4) = zero
            lobCoor(5) = -lobCoor(3)
            lobCoor(6) = -lobCoor(2)
            lobCoor(7) = -lobCoor(1)
            lobWeight(1) = un/21.d0
            lobWeight(2) = (124.d0-7.d0*rac15)/350.d0
            lobWeight(3) = (124.d0+7.d0*rac15)/350.d0
            lobWeight(4) = 256.d0/525.d0
            lobWeight(5) = lobWeight(3)
            lobWeight(6) = lobWeight(2)
            lobWeight(7) = lobWeight(1)
            do iz = 1, 7
                xpg(iz) = un/3.d0
                ypg(iz) = un/3.d0
                zpg(iz) = lobCoor(iz)
                hpg(iz) = lobWeight(iz)
            end do

        else
            ASSERT(ASTER_FALSE)
        end if

    else if (elrefa .eq. 'TE4' .or. elrefa .eq. 'T10' .or. elrefa .eq. 'T15') then
        if (fapg .eq. 'FPG4') then
! --------- FORMULE A 4 POINTS :  (CF TOUZOT PAGE 300) - ORDRE 2 EN X Y Z
            a1 = (5.d0-rac5)/20.d0
            b1 = (5.d0+3.d0*rac5)/20.d0
            h5 = un/24.d0
            npi = 0
            do i = 1, 4
                npi = npi+1
                xpg(npi) = a1
                ypg(npi) = a1
                zpg(npi) = a1
                hpg(npi) = h5
            end do
            zpg(2) = b1
            ypg(3) = b1
            xpg(4) = b1

        else if (fapg .eq. 'FPG5') then
! --------- FORMULE A 5 POINTS :  (CF TOUZOT PAGE 300) - ORDRE 3 EN X Y Z
            a1 = 0.25d0
            b1 = un/6.d0
            c1 = undemi
            h1 = -deux/15.d0
            h2 = 3.d0/40.d0
            xpg(1) = a1
            ypg(1) = a1
            zpg(1) = a1
            hpg(1) = h1
            do i = 2, 5
                xpg(i) = b1
                ypg(i) = b1
                zpg(i) = b1
                hpg(i) = h2
            end do
            zpg(3) = c1
            ypg(4) = c1
            xpg(5) = c1

        else if (fapg .eq. 'FPG11') then
! --------- FORMULE A 11 POINTS :  (CF DUNAVANT) - ORDRE 4 EN X Y Z
            xpg(1) = 0.097204644587583d0
            ypg(1) = 0.106604172561993d0
            zpg(1) = 0.684390415453040d0
            hpg(1) = 0.106468034155490d0/6.d0
            xpg(2) = 0.029569495206479d0
            ypg(2) = 0.329232959742646d0
            zpg(2) = 0.317903560213394d0
            hpg(2) = 0.110234232428497d0/6.d0
            xpg(3) = 0.432710239047768d0
            ypg(3) = 0.103844116410993d0
            zpg(3) = 0.353823239209297d0
            hpg(3) = 0.154976116016246d0/6.d0
            xpg(4) = 0.240276664928072d0
            ypg(4) = 0.304448402434496d0
            zpg(4) = 0.126801725915392d0
            hpg(4) = 0.193410812049634d0/6.d0
            xpg(5) = 0.129411373788910d0
            ypg(5) = 0.538007203916185d0
            zpg(5) = 0.330190414837464d0
            hpg(5) = 0.076162715245558d0/6.d0
            xpg(6) = 0.121541991333927d0
            ypg(6) = 0.008991260093335d0
            zpg(6) = 0.306493988429690d0
            hpg(6) = 0.079426680068025d0/6.d0
            xpg(7) = 0.450765876091276d0
            ypg(7) = 0.432953490481355d0
            zpg(7) = 0.059456616299433d0
            hpg(7) = 0.069469965937635d0/6.d0
            xpg(8) = 0.419266313879513d0
            ypg(8) = 0.053341239535745d0
            zpg(8) = 0.047781435559086d0
            hpg(8) = 0.059933185146559d0/6.d0
            xpg(9) = 0.067223294893383d0
            ypg(9) = 0.741228882093622d0
            zpg(9) = 0.035183929773598d0
            hpg(9) = 0.055393798871576d0/6.d0
            xpg(10) = 0.752508507009654d0
            ypg(10) = 0.081404918402859d0
            zpg(10) = 0.068099370938206d0
            hpg(10) = 0.055273369155936d0/6.d0
            xpg(11) = 0.040490506727590d0
            ypg(11) = 0.174694058697230d0
            zpg(11) = 0.013560701879802d0
            hpg(11) = 0.039251090924839d0/6.d0

        else if (fapg .eq. 'FPG15') then
! --------- FORMULE A 15 POINTS :  (CF TOUZOT PAGE 300) - ORDRE 5 EN X Y Z
            a1 = 0.25d0
            b1 = (7.0d0+rac15)/34.0d0
            b2 = (7.0d0-rac15)/34.0d0
            c1 = (13.0d0-3.0d0*rac15)/34.0d0
            c2 = (13.0d0+3.0d0*rac15)/34.0d0
            d1 = (5.0d0-rac15)/20.0d0
            e1 = (5.0d0+rac15)/20.0d0
            h5 = 8.0d0/405.0d0
            h1 = (2665.0d0-14.0d0*rac15)/226800.0d0
            h2 = (2665.0d0+14.0d0*rac15)/226800.0d0
            h3 = 5.0d0/567.0d0

            xpg(1) = a1
            ypg(1) = a1
            zpg(1) = a1
            hpg(1) = h5
            xpg(2) = b1
            ypg(2) = b1
            zpg(2) = b1
            hpg(2) = h1
            xpg(3) = b1
            ypg(3) = b1
            zpg(3) = c1
            hpg(3) = h1
            xpg(4) = b1
            ypg(4) = c1
            zpg(4) = b1
            hpg(4) = h1
            xpg(5) = c1
            ypg(5) = b1
            zpg(5) = b1
            hpg(5) = h1
            xpg(6) = b2
            ypg(6) = b2
            zpg(6) = b2
            hpg(6) = h2
            xpg(7) = b2
            ypg(7) = b2
            zpg(7) = c2
            hpg(7) = h2
            xpg(8) = b2
            ypg(8) = c2
            zpg(8) = b2
            hpg(8) = h2
            xpg(9) = c2
            ypg(9) = b2
            zpg(9) = b2
            hpg(9) = h2
            xpg(10) = d1
            ypg(10) = d1
            zpg(10) = e1
            hpg(10) = h3
            xpg(11) = d1
            ypg(11) = e1
            zpg(11) = d1
            hpg(11) = h3
            xpg(12) = e1
            ypg(12) = d1
            zpg(12) = d1
            hpg(12) = h3
            xpg(13) = d1
            ypg(13) = e1
            zpg(13) = e1
            hpg(13) = h3
            xpg(14) = e1
            ypg(14) = d1
            zpg(14) = e1
            hpg(14) = h3
            xpg(15) = e1
            ypg(15) = e1
            zpg(15) = d1
            hpg(15) = h3

        else if (fapg .eq. 'FPG23') then
! --------- FORMULE A 23 POINTS :  (CF DUNAVANT) - ORDRE 6 EN X Y Z
            xpg(1) = 0.038836084344884d0
            ypg(1) = 0.024318974248143d0
            zpg(1) = 0.902928799013611d0
            hpg(1) = 0.001182632475277d0
            xpg(2) = 0.064769436930053d0
            ypg(2) = 0.267844198183576d0
            zpg(2) = 0.636767508558514d0
            hpg(2) = 0.005251568313784d0
            xpg(3) = 0.064775160447105d0
            ypg(3) = 0.023467795573055d0
            zpg(3) = 0.390862050671012d0
            hpg(3) = 0.004038547812907d0
            xpg(4) = 0.277903669330078d0
            ypg(4) = 0.063732895294998d0
            zpg(4) = 0.594909689021796d0
            hpg(4) = 0.008148345983740d0
            xpg(5) = 0.066098662414680d0
            ypg(5) = 0.083678814060055d0
            zpg(5) = 0.630054555110990d0
            hpg(5) = 0.008838887318028d0
            xpg(6) = 0.325119658577025d0
            ypg(6) = 0.329379718549198d0
            zpg(6) = 0.326833504619046d0
            hpg(6) = 0.007206549449246d0
            xpg(7) = 0.319194280348931d0
            ypg(7) = 0.304169265349782d0
            zpg(7) = 0.044383344357208d0
            hpg(7) = 0.011189302702093d0
            xpg(8) = 0.328388171231222d0
            ypg(8) = 0.038288670738245d0
            zpg(8) = 0.320287433697693d0
            hpg(8) = 0.009970224610238d0
            xpg(9) = 0.055099022490726d0
            ypg(9) = 0.351939197334705d0
            zpg(9) = 0.381084308906310d0
            hpg(9) = 0.010435745880219d0
            xpg(10) = 0.124649963637486d0
            ypg(10) = 0.152103811309931d0
            zpg(10) = 0.201234567364421d0
            hpg(10) = 0.010722336995515d0
            xpg(11) = 0.065924923160010d0
            ypg(11) = 0.624321363553430d0
            zpg(11) = 0.253593674743200d0
            hpg(11) = 0.007265066343438d0
            xpg(12) = 0.007354523838069d0
            ypg(12) = 0.211297658581586d0
            zpg(12) = 0.251184495277530d0
            hpg(12) = 0.003760944546357d0
            xpg(13) = 0.617455720147269d0
            ypg(13) = 0.063199980942570d0
            zpg(13) = 0.258449148983926d0
            hpg(13) = 0.007768855687763d0
            xpg(14) = 0.279420052945988d0
            ypg(14) = 0.255820784264986d0
            zpg(14) = 0.269569929633272d0
            hpg(14) = 0.018766567415678d0
            xpg(15) = 0.287725094826464d0
            ypg(15) = 0.577345781389727d0
            zpg(15) = 0.064620638073369d0
            hpg(15) = 0.008989168438052d0
            xpg(16) = 0.594717301875796d0
            ypg(16) = 0.065177992763370d0
            zpg(16) = 0.066603298007603d0
            hpg(16) = 0.008294771681919d0
            xpg(17) = 0.066789599781738d0
            ypg(17) = 0.530063275481017d0
            zpg(17) = 0.076992717100967d0
            hpg(17) = 0.010511060314253d0
            xpg(18) = 0.626540201708882d0
            ypg(18) = 0.248449540118895d0
            zpg(18) = 0.062115533183599d0
            hpg(18) = 0.007858005078710d0
            xpg(19) = 0.060010583020269d0
            ypg(19) = 0.213041183236186d0
            zpg(19) = 0.025842686260703d0
            hpg(19) = 0.004250720711174d0
            xpg(20) = 0.275786300469851d0
            ypg(20) = 0.053996140835915d0
            zpg(20) = 0.060016149166169d0
            hpg(20) = 0.006619016274847d0
            xpg(21) = 0.051325206165203d0
            ypg(21) = 0.841138951662319d0
            zpg(21) = 0.037264752138356d0
            hpg(21) = 0.002654246530834d0
            xpg(22) = 0.040576051066818d0
            ypg(22) = 0.008781957777519d0
            zpg(22) = 0.088600350468910d0
            hpg(22) = 0.001737222620616d0
            xpg(23) = 0.903770001332182d0
            ypg(23) = 0.022865823814023d0
            zpg(23) = 0.029335721083179d0
            hpg(23) = 0.001206879481978d0

        else if (fapg .eq. 'FPG4NOS') then
! --------- POUR LES POINTS DE GAUSS
            a1 = (5.d0-rac5)/20.d0
            b1 = (5.d0+3.d0*rac5)/20.d0
            h5 = un/24.d0
            npi = 0
            do i = 1, 4
                npi = npi+1
                xpg(npi) = a1
                ypg(npi) = a1
                zpg(npi) = a1
                hpg(npi) = h5
            end do
            zpg(2) = b1
            ypg(3) = b1
            xpg(4) = b1
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+4) = vol/nnos
                xpg(iNode+4) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+4) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+4) = xno(ndim*(iNode-1)+3)
            end do
        end if

    else if (elrefa .eq. 'PY5' .or. elrefa .eq. 'P13' .or. elrefa .eq. 'P19') then
        if (fapg .eq. 'FPG5') then
! --------- ORDRE 2
            p1 = deux/15.d0
            h1 = 0.1531754163448146d0
            h2 = 0.6372983346207416d0
            xpg(1) = undemi
            xpg(2) = zero
            xpg(3) = -undemi
            xpg(4) = zero
            xpg(5) = zero
            ypg(1) = zero
            ypg(2) = undemi
            ypg(3) = zero
            ypg(4) = -undemi
            ypg(5) = zero
            zpg(1) = h1
            zpg(2) = h1
            zpg(3) = h1
            zpg(4) = h1
            zpg(5) = h2
            hpg(1) = p1
            hpg(2) = p1
            hpg(3) = p1
            hpg(4) = p1
            hpg(5) = p1

        else if (fapg .eq. 'FPG6') then
!    Order 3 (REX 33086)
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.

            xpg(1) = 0.d0
            ypg(1) = 0.d0
            zpg(1) = 0.5714285703860683d0
            hpg(1) = 0.1681372559485071d0

            xpg(2) = 0.d0
            ypg(2) = 0.d0
            zpg(2) = 5.67585d-09
            hpg(2) = 0.07500000404404333d0

            xpg(3) = -0.5610836110587396d0
            ypg(3) = 0.d0
            zpg(3) = 0.1666666666666667d0
            hpg(3) = 0.1058823516685291d0

            xpg(4) = 0.d0
            ypg(4) = -0.5610836110587396d0
            zpg(4) = 0.1666666666666667d0
            hpg(4) = 0.1058823516685291d0

            xpg(5) = 0.d0
            ypg(5) = 0.5610836110587396d0
            zpg(5) = 0.1666666666666667d0
            hpg(5) = 0.1058823516685291d0

            xpg(6) = 0.5610836110587396d0
            ypg(6) = 0.d0
            zpg(6) = 0.1666666666666667d0
            hpg(6) = 0.1058823516685291d0

        else if (fapg .eq. 'FPG10') then
! Order 4 (REX 20813)
!    Freddie Witherden, Peter Vincent,
!    On the identification of symmetric quadrature rules for finite element methods,
!    Computers and Mathematics with Applications,
!    Volume 69, pages 1232-1241, 2015.

            xpg(1) = 0.d0
            ypg(1) = 0.d0
            zpg(1) = 0.6772327888861374d0
            hpg(1) = 0.07582792211376127d0

            xpg(2) = 0.d0
            ypg(2) = 0.d0
            zpg(2) = 0.1251369531087465d0
            hpg(2) = 0.1379222683930349d0

            xpg(3) = -0.3252907781991163d0
            ypg(3) = -0.3252907781991163d0
            zpg(3) = 0.3223841495782137d0
            hpg(3) = 0.07088305859288367d0

        else if (fapg .eq. 'FPG23') then
            ! ORDRE 6
            ! x = (xi - eta) / 2
            ! y = (xi + eta) / 2
            ! z = (zeta + 1) / 2
            ! h = weight / 4

            xpg(5) = 0.3252907781991163d0
            ypg(5) = 0.3252907781991163d0
            zpg(5) = 0.3223841495782137d0
            hpg(5) = 0.07088305859288367d0

            xpg(6) = 0.3252907781991163d0
            ypg(6) = -0.3252907781991163d0
            zpg(6) = 0.3223841495782137d0
            hpg(6) = 0.07088305859288367d0

            xpg(7) = -0.65796699712169d0
            ypg(7) = 0.d0
            zpg(7) = 0.0392482838988154d0
            hpg(7) = 0.04234606044708394d0

            xpg(8) = 0.d0
            ypg(8) = -0.65796699712169d0
            zpg(8) = 0.0392482838988154d0
            hpg(8) = 0.04234606044708394d0

            xpg(9) = 0.d0
            ypg(9) = 0.65796699712169d0
            zpg(9) = 0.0392482838988154d0
            hpg(9) = 0.04234606044708394d0

            xpg(10) = 0.65796699712169d0
            ypg(10) = 0.d0
            zpg(10) = 0.0392482838988154d0
            hpg(10) = 0.04234606044708394d0

        else if (fapg .eq. 'FPG5NOS') then
! --------- POUR LES POINTS DE GAUSS
            p1 = deux/15.d0
            h1 = 0.1531754163448146d0
            h2 = 0.6372983346207416d0
            xpg(1) = undemi
            xpg(2) = zero
            xpg(3) = -undemi
            xpg(4) = zero
            xpg(5) = zero
            ypg(1) = zero
            ypg(2) = undemi
            ypg(3) = zero
            ypg(4) = -undemi
            ypg(5) = zero
            zpg(1) = h1
            zpg(2) = h1
            zpg(3) = h1
            zpg(4) = h1
            zpg(5) = h2
            hpg(1) = p1
            hpg(2) = p1
            hpg(3) = p1
            hpg(4) = p1
            hpg(5) = p1
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+5) = vol/nnos
                xpg(iNode+5) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+5) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+5) = xno(ndim*(iNode-1)+3)
            end do
        else
            ASSERT(ASTER_FALSE)

        end if

    else if (elrefa .eq. 'TR3' .or. elrefa .eq. 'TR6' .or. elrefa .eq. 'TR7') then
        if (fapg .eq. 'FPG1') then
            xpg(1) = un/3.d0
            ypg(1) = un/3.d0
            hpg(1) = un/deux

        else if (fapg .eq. 'FPG3') then
            xpg(1) = un/6.d0
            ypg(1) = un/6.d0
            xpg(2) = deux/3.d0
            ypg(2) = un/6.d0
            xpg(3) = un/6.d0
            ypg(3) = deux/3.d0
            hpg(1) = un/6.d0
            hpg(2) = un/6.d0
            hpg(3) = un/6.d0

        else if (fapg .eq. 'FPG4') then
            xpg(1) = 0.2d0
            ypg(1) = 0.2d0
            xpg(2) = 0.6d0
            ypg(2) = 0.2d0
            xpg(3) = 0.2d0
            ypg(3) = 0.6d0
            xpg(4) = un/3.d0
            ypg(4) = un/3.d0
            hpg(1) = 25.d0/96.d0
            hpg(2) = 25.d0/96.d0
            hpg(3) = 25.d0/96.d0
            hpg(4) = -27.d0/96.d0

        else if (fapg .eq. 'FPG6') then
            h1 = 0.111690794839005d0
            h2 = 0.054975871827661d0
            a1 = 0.445948490915965d0
            b1 = 0.091576213509771d0
            xpg(3) = (t(b1)+un)/deux
            ypg(3) = (t(un-deux*b1)+un)/deux
            xpg(1) = (t(b1)+un)/deux
            ypg(1) = (t(b1)+un)/deux
            xpg(2) = (t(un-deux*b1)+un)/deux
            ypg(2) = (t(b1)+un)/deux
            xpg(6) = (t(un-deux*a1)+un)/deux
            ypg(6) = (t(a1)+un)/deux
            xpg(4) = (t(a1)+un)/deux
            ypg(4) = (t(un-deux*a1)+un)/deux
            xpg(5) = (t(a1)+un)/deux
            ypg(5) = (t(a1)+un)/deux
            hpg(1) = h2
            hpg(2) = h2
            hpg(3) = h2
            hpg(4) = h1
            hpg(5) = h1
            hpg(6) = h1

        else if (fapg .eq. 'FPG7') then
            p1 = (155.d0+rac15)/2400.d0
            p2 = (155.d0-rac15)/2400.d0
            a2 = (6.d0+rac15)/21.d0
            b2 = (6.d0-rac15)/21.d0
            xpg(1) = un/3.d0
            ypg(1) = un/3.d0
            xpg(2) = a2
            ypg(2) = a2
            xpg(3) = un-deux*a2
            ypg(3) = a2
            xpg(4) = a2
            ypg(4) = un-deux*a2
            xpg(5) = b2
            ypg(5) = b2
            xpg(6) = un-deux*b2
            ypg(6) = b2
            xpg(7) = b2
            ypg(7) = un-deux*b2
            hpg(1) = 9.d0/80.d0
            hpg(2) = p1
            hpg(3) = p1
            hpg(4) = p1
            hpg(5) = p2
            hpg(6) = p2
            hpg(7) = p2

        else if (fapg .eq. 'FPG12') then
            a1 = 0.063089014491502d0
            b1 = 0.249286745170910d0
            c1 = 0.310352451033785d0
            d1 = 0.053145049844816d0
            xpg(1) = a1
            ypg(1) = a1
            xpg(2) = un-deux*a1
            ypg(2) = a1
            xpg(3) = a1
            ypg(3) = un-deux*a1
            xpg(4) = b1
            ypg(4) = b1
            xpg(5) = un-deux*b1
            ypg(5) = b1
            xpg(6) = b1
            ypg(6) = un-deux*b1
            xpg(7) = c1
            ypg(7) = d1
            xpg(8) = d1
            ypg(8) = c1
            xpg(9) = un-c1-d1
            ypg(9) = c1
            xpg(10) = un-c1-d1
            ypg(10) = d1
            xpg(11) = c1
            ypg(11) = un-c1-d1
            xpg(12) = d1
            ypg(12) = un-c1-d1
            p1 = 0.025422453185103d0
            p2 = 0.058393137863189d0
            p3 = 0.041425537809187d0
            hpg(1) = p1
            hpg(2) = p1
            hpg(3) = p1
            hpg(4) = p2
            hpg(5) = p2
            hpg(6) = p2
            hpg(7) = p3
            hpg(8) = p3
            hpg(9) = p3
            hpg(10) = p3
            hpg(11) = p3
            hpg(12) = p3

        else if (fapg .eq. 'FPG13') then
!         FORMULE A 13 POINTS : ORDRE 7  (CF BATHE, PAGE 280)
            xpg(1) = 0.0651301029022d0
            ypg(1) = 0.0651301029022d0
            xpg(2) = 0.8697397941956d0
            ypg(2) = 0.0651301029022d0
            xpg(3) = 0.0651301029022d0
            ypg(3) = 0.8697397941956d0
            xpg(4) = 0.3128654960049d0
            ypg(4) = 0.0486903154253d0
            xpg(5) = 0.6384441885698d0
            ypg(5) = 0.3128654960049d0
            xpg(6) = 0.0486903154253d0
            ypg(6) = 0.6384441885698d0
            xpg(7) = 0.6384441885698d0
            ypg(7) = 0.0486903154253d0
            xpg(8) = 0.3128654960049d0
            ypg(8) = 0.6384441885698d0
            xpg(9) = 0.0486903154253d0
            ypg(9) = 0.3128654960049d0
            xpg(10) = 0.2603459660790d0
            ypg(10) = 0.2603459660790d0
            xpg(11) = 0.4793080678419d0
            ypg(11) = 0.2603459660790d0
            xpg(12) = 0.2603459660790d0
            ypg(12) = 0.4793080678419d0
            xpg(13) = 0.3333333333333d0
            ypg(13) = 0.3333333333333d0
            p1 = 0.0533472356088d0/deux
            p2 = 0.0771137608903d0/deux
            p3 = 0.1756152574332d0/deux
            p4 = -0.1495700444677d0/deux
            hpg(1) = p1
            hpg(2) = p1
            hpg(3) = p1
            hpg(4) = p2
            hpg(5) = p2
            hpg(6) = p2
            hpg(7) = p2
            hpg(8) = p2
            hpg(9) = p2
            hpg(10) = p3
            hpg(11) = p3
            hpg(12) = p3
            hpg(13) = p4

        else if (fapg .eq. 'FPG16') then
            xpg(1) = un/3.d0
            ypg(1) = un/3.d0
            xpg(2) = 0.081414823414554d0
            ypg(2) = 0.459292588292723d0
            xpg(3) = 0.459292588292723d0
            ypg(3) = 0.081414823414554d0
            xpg(4) = 0.459292588292723d0
            ypg(4) = 0.459292588292723d0
            xpg(5) = 0.658861384496480d0
            ypg(5) = 0.170569307751760d0
            xpg(6) = 0.170569307751760d0
            ypg(6) = 0.658861384496480d0
            xpg(7) = 0.170569307751760d0
            ypg(7) = 0.170569307751760d0
            xpg(8) = 0.898905543365938d0
            ypg(8) = 0.050547228317031d0
            xpg(9) = 0.050547228317031d0
            ypg(9) = 0.898905543365938d0
            xpg(10) = 0.050547228317031d0
            ypg(10) = 0.050547228317031d0
            xpg(11) = 0.008394777409958d0
            ypg(11) = 0.728492392955404d0
            xpg(12) = 0.728492392955404d0
            ypg(12) = 0.008394777409958d0
            xpg(13) = 0.263112829634638d0
            ypg(13) = 0.008394777409958d0
            xpg(14) = 0.008394777409958d0
            ypg(14) = 0.263112829634638d0
            xpg(15) = 0.263112829634638d0
            ypg(15) = 0.728492392955404d0
            xpg(16) = 0.728492392955404d0
            ypg(16) = 0.263112829634638d0
            p1 = 0.144315607677787d0/deux
            p2 = 0.095091634267285d0/deux
            p3 = 0.103217370534718d0/deux
            p4 = 0.032458497623198d0/deux
            p5 = 0.027230314174435d0/deux
            hpg(1) = p1
            hpg(2) = p2
            hpg(3) = p2
            hpg(4) = p2
            hpg(5) = p3
            hpg(6) = p3
            hpg(7) = p3
            hpg(8) = p4
            hpg(9) = p4
            hpg(10) = p4
            hpg(11) = p5
            hpg(12) = p5
            hpg(13) = p5
            hpg(14) = p5
            hpg(15) = p5
            hpg(16) = p5

        else if (fapg .eq. 'COT3') then
            xpg(1) = undemi
            ypg(1) = undemi
            xpg(2) = zero
            ypg(2) = undemi
            xpg(3) = undemi
            ypg(3) = zero
            hpg(1) = un/6.d0
            hpg(2) = un/6.d0
            hpg(3) = un/6.d0

        else if (fapg .eq. 'SIMP') then
            xpg(1) = zero
            ypg(1) = zero
            xpg(2) = un
            ypg(2) = zero
            xpg(3) = zero
            ypg(3) = un
            xpg(4) = undemi
            ypg(4) = zero
            xpg(5) = undemi
            ypg(5) = undemi
            xpg(6) = zero
            ypg(6) = undemi
            hpg(1) = un/30.d0
            hpg(2) = un/30.d0
            hpg(3) = un/30.d0
            hpg(4) = 4.d0/30.d0
            hpg(5) = 4.d0/30.d0
            hpg(6) = 4.d0/30.d0

        else if (fapg .eq. 'FPG3NOS') then
! --------- POUR LES POINTS DE GAUSS
            xpg(1) = un/6.d0
            ypg(1) = un/6.d0
            xpg(2) = deux/3.d0
            ypg(2) = un/6.d0
            xpg(3) = un/6.d0
            ypg(3) = deux/3.d0
            hpg(1) = un/6.d0
            hpg(2) = un/6.d0
            hpg(3) = un/6.d0
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+3) = vol/nnos
                xpg(iNode+3) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+3) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+3) = xno(ndim*(iNode-1)+3)
            end do

        else
            ASSERT(ASTER_FALSE)

        end if

    else if (elrefa .eq. 'QU4' .or. elrefa .eq. 'QU8' .or. elrefa .eq. 'QU9') then
        if (fapg .eq. 'FPG1') then
            xpg(1) = zero
            ypg(1) = zero
            hpg(1) = 4.d0

        else if (fapg .eq. 'FIS2') then
! ------- ELEMENT PARTICULIER DE FISSURE, S'APPUIE SUR UN SEG2
            xpg(1) = -rac_1div3
            ypg(1) = zero
            xpg(2) = rac_1div3
            ypg(2) = zero
            hpg(1) = deux
            hpg(2) = deux

        else if (fapg .eq. 'FPG4') then
            xpg(1) = -rac_1div3
            ypg(1) = -rac_1div3
            xpg(2) = rac_1div3
            ypg(2) = -rac_1div3
            xpg(3) = rac_1div3
            ypg(3) = rac_1div3
            xpg(4) = -rac_1div3
            ypg(4) = rac_1div3
            hpg(1) = un
            hpg(2) = un
            hpg(3) = un
            hpg(4) = un

        else if (fapg .eq. 'FPG9') then
            hpg(1) = 25.d0/81.0d0
            hpg(2) = 25.d0/81.0d0
            hpg(3) = 25.d0/81.0d0
            hpg(4) = 25.d0/81.0d0
            hpg(5) = 40.d0/81.0d0
            hpg(6) = 40.d0/81.0d0
            hpg(7) = 40.d0/81.0d0
            hpg(8) = 40.d0/81.0d0
            hpg(9) = 64.d0/81.0d0
            xpg(1) = -rac_3div5
            ypg(1) = -rac_3div5
            xpg(2) = rac_3div5
            ypg(2) = -rac_3div5
            xpg(3) = rac_3div5
            ypg(3) = rac_3div5
            xpg(4) = -rac_3div5
            ypg(4) = rac_3div5
            xpg(5) = zero
            ypg(5) = -rac_3div5
            xpg(6) = rac_3div5
            ypg(6) = zero
            xpg(7) = zero
            ypg(7) = rac_3div5
            xpg(8) = -rac_3div5
            ypg(8) = zero
            xpg(9) = zero
            ypg(9) = zero

        else if (fapg .eq. 'FPG9COQ') then
            hpg(7) = 25.d0/81.0d0
            hpg(1) = 25.d0/81.0d0
            hpg(3) = 25.d0/81.0d0
            hpg(5) = 25.d0/81.0d0
            hpg(8) = 40.d0/81.0d0
            hpg(2) = 40.d0/81.0d0
            hpg(4) = 40.d0/81.0d0
            hpg(6) = 40.d0/81.0d0
            hpg(9) = 64.d0/81.0d0
            xpg(1) = -rac_3div5
            ypg(1) = -rac_3div5
            xpg(3) = rac_3div5
            ypg(3) = -rac_3div5
            xpg(5) = rac_3div5
            ypg(5) = rac_3div5
            xpg(7) = -rac_3div5
            ypg(7) = rac_3div5
            xpg(2) = zero
            ypg(2) = -rac_3div5
            xpg(4) = rac_3div5
            ypg(4) = zero
            xpg(6) = zero
            ypg(6) = rac_3div5
            xpg(8) = -rac_3div5
            ypg(8) = zero
            xpg(9) = zero
            ypg(9) = zero

        else if (fapg .eq. 'FPG16') then
            h(1) = (18.d0+rac30)/36.d0
            h(2) = h(1)
            h(3) = (18.d0-rac30)/36.d0
            h(4) = h(3)
            a(1) = -gauss4p12
            a(2) = -a(1)
            a(3) = -gauss4p34
            a(4) = -a(3)
            npar = 4
            npi = 0
            do ix = 1, npar
                do iy = 1, npar
                    npi = npi+1
                    xpg(npi) = a(ix)
                    ypg(npi) = a(iy)
                    hpg(npi) = h(ix)*h(iy)
                end do
            end do

        else if (fapg .eq. 'FPG4NOS') then
! --------- POUR LES POINTS DE GAUSS
            xpg(1) = -rac_1div3
            ypg(1) = -rac_1div3
            xpg(2) = rac_1div3
            ypg(2) = -rac_1div3
            xpg(3) = rac_1div3
            ypg(3) = rac_1div3
            xpg(4) = -rac_1div3
            ypg(4) = rac_1div3
            hpg(1) = un
            hpg(2) = un
            hpg(3) = un
            hpg(4) = un
! --------- POUR LES SOMMETS
            do iNode = 1, nnos
                hpg(iNode+4) = vol/nnos
                xpg(iNode+4) = xno(ndim*(iNode-1)+1)
                if (ndim .ge. 2) ypg(iNode+4) = xno(ndim*(iNode-1)+2)
                if (ndim .eq. 3) zpg(iNode+4) = xno(ndim*(iNode-1)+3)
            end do

        else
            ASSERT(ASTER_FALSE)

        end if

    else if (elrefa .eq. 'SE2' .or. elrefa .eq. 'SE3' .or. elrefa .eq. 'SE4') then
        if (fapg .eq. 'FPG1') then
            xpg(1) = zero
            hpg(1) = deux

        else if (fapg .eq. 'FPG2') then
            xpg(1) = rac_1div3
            xpg(2) = -xpg(1)
            hpg(1) = un
            hpg(2) = hpg(1)

        else if (fapg .eq. 'FPG3') then
            xpg(1) = -rac_3div5
            xpg(2) = zero
            xpg(3) = rac_3div5
            hpg(1) = 5.d0/9.d0
            hpg(2) = 8.d0/9.d0
            hpg(3) = 5.d0/9.d0

        else if (fapg .eq. 'FPG4') then
            xpg(1) = gauss4p12
            xpg(2) = -xpg(1)
            xpg(3) = gauss4p34
            xpg(4) = -xpg(3)
            hpg(1) = (18.d0+rac30)/36.d0
            hpg(2) = hpg(1)
            hpg(3) = (18.d0-rac30)/36.d0
            hpg(4) = hpg(3)

        else if (fapg .eq. 'FPG2NOS') then
            xpg(1) = rac_1div3
            xpg(2) = -xpg(1)
            xpg(3) = xno(1)
            xpg(4) = xno(2)
            hpg(1) = un
            hpg(2) = hpg(1)
            hpg(3) = vol/nnos
            hpg(4) = hpg(3)

        else if (fapg .eq. 'FPG3NOS') then
            xpg(1) = -rac_3div5
            xpg(2) = zero
            xpg(3) = rac_3div5
            xpg(4) = xno(1)
            xpg(5) = xno(nnos)
            hpg(1) = 5.d0/9.d0
            hpg(2) = 8.d0/9.d0
            hpg(3) = 5.d0/9.d0
            hpg(4) = vol/nnos
            hpg(5) = hpg(4)

        else if (fapg .eq. 'SIMP') then
            xpg(1) = -un
            xpg(2) = zero
            xpg(3) = un
            hpg(1) = un/3.d0
            hpg(2) = 4.d0/3.d0
            hpg(3) = un/3.d0

        else if (fapg .eq. 'SIMP1') then
            xpg(1) = -un
            xpg(2) = -undemi
            xpg(3) = zero
            xpg(4) = undemi
            xpg(5) = un
            hpg(1) = un/6.d0
            hpg(2) = deux/3.d0
            hpg(3) = un/3.d0
            hpg(4) = deux/3.d0
            hpg(5) = un/6.d0

        else if (fapg .eq. 'COTES') then
            xpg(1) = -un
            xpg(2) = -un/3.d0
            xpg(3) = un/3.d0
            xpg(4) = un
            hpg(1) = un/4.d0
            hpg(2) = 3.d0/4.d0
            hpg(3) = 3.d0/4.d0
            hpg(4) = un/4.d0

        else if (fapg .eq. 'COTES1') then
            xpg(1) = -un
            xpg(2) = -un/deux
            xpg(3) = zero
            xpg(4) = un/deux
            xpg(5) = un
            hpg(1) = 7.d0/45.d0
            hpg(2) = 32.d0/45.d0
            hpg(3) = 12.d0/45.d0
            hpg(4) = 32.d0/45.d0
            hpg(5) = 7.d0/45.d0

        else if (fapg .eq. 'COTES2') then
            xpg(1) = -un
            xpg(2) = -7.d0/9.d0
            xpg(3) = -5.d0/9.d0
            xpg(4) = -un/3.d0
            xpg(5) = -un/9.d0
            xpg(6) = un/9.d0
            xpg(7) = un/3.d0
            xpg(8) = 5.d0/9.d0
            xpg(9) = 7.d0/9.d0
            xpg(10) = un
            hpg(1) = un/12.d0
            hpg(2) = un/4.d0
            hpg(3) = un/4.d0
            hpg(4) = un/6.d0
            hpg(5) = un/4.d0
            hpg(6) = un/4.d0
            hpg(7) = un/6.d0
            hpg(8) = un/4.d0
            hpg(9) = un/4.d0
            hpg(10) = un/12.d0

        else
            ASSERT(ASTER_FALSE)

        end if

    else if (elrefa .eq. 'PO1') then
        hpg(1) = un

    else
        ASSERT(ASTER_FALSE)
    end if
!
170 continue
!
    do i = 1, nbpg
        poipg(i) = hpg(i)
        if (ndim .ge. 1) coopg(ndim*(i-1)+1) = xpg(i)
        if (ndim .ge. 2) coopg(ndim*(i-1)+2) = ypg(i)
        if (ndim .eq. 3) coopg(ndim*(i-1)+3) = zpg(i)
    end do

end subroutine
