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
subroutine efcoq3d(nomte, nb1, nb2, cara, geom, &
                   lzr, chg, matr, effg, nbcou, &
                   npgsn, npgsr, npge, nso, npgt)
!     CALCUL DE EFGE_ELNO
!     ------------------------------------------------------------------
    implicit none
!
#include "asterfort/vdefgn.h"
#include "asterfort/vdefro.h"
#include "asterfort/vdrepe.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
!
!
    character(len=16) :: nomte
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ic, icomp, ii
    integer(kind=8) :: inte, intsn, intsr, isom
    integer(kind=8) :: j, jj
    integer(kind=8) :: k, k1, kpgs, l
    integer(kind=8) :: nbcou, ncmp
    integer(kind=8) :: npge, npgt
    integer(kind=8) :: nso
    real(kind=8) :: hic, s, zero, zic, zmin
!-----------------------------------------------------------------------
    integer(kind=8) :: icou
    integer(kind=8) :: nb1, nb2, npgsr, npgsn
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: epais
    real(kind=8) :: matevn(2, 2, npgt), matevg(2, 2, npgt)
    real(kind=8) :: geom(*), cara(*), chg(*), matr(*), effg(*), lzr(*)
    real(kind=8) :: sigm(6, 270), sigma(6, 120), effgc(8, 9), effgt(8, 9)
!
!
    zero = 0.0d0
!
    epais = cara(1)
    zmin = -epais/2.d0
    hic = epais/nbcou
!
    call vectan(nb1, nb2, geom, lzr, vecta, &
                vectn, vectpt)
!
    kpgs = 0
    do icou = 1, nbcou
        do inte = 1, npge
            if (inte .eq. 1) then
                zic = zmin+(icou-1)*hic
            else if (inte .eq. 2) then
                zic = zmin+hic/2.d0+(icou-1)*hic
            else
                zic = zmin+hic+(icou-1)*hic
            end if
!
            do intsn = 1, npgsn
                kpgs = kpgs+1
                k1 = 6*((intsn-1)*npge*nbcou+(icou-1)*npge+inte-1)
                do i = 1, 6
                    sigm(i, kpgs) = chg(k1+i)
                end do
            end do
        end do
    end do
    ncmp = 6
!
!
! --- DETERMINATION DES REPERES  LOCAUX DE L'ELEMENT AUX POINTS
! --- D'INTEGRATION ET STOCKAGE DE CES REPERES DANS LE VECTEUR .DESR
!     --------------------------------------------------------------
    k = 0
    do intsr = 1, npgsr
        call vectgt(0, nb1, geom, zero, intsr, &
                    lzr, epais, vectn, vectg, vectt)
!
        do j = 1, 3
            do i = 1, 3
                k = k+1
                lzr(2000+k) = vectt(i, j)
            end do
        end do
    end do
! !
! !--- EXTRAPOLATION VERS LES NOEUDS SOMMETS
! !
!
!
    do icou = 1, nbcou
        do ic = 1, ncmp
            do i = 1, npge*nso
                l = npge*npgsn*(i-1)
                s = 0.d0
                do j = 1, npge*npgsn
                    jj = (icou-1)*npge*npgsn+j
                    s = s+matr(l+j)*sigm(ic, jj)
                end do
                ii = (icou-1)*npge*nso+i
                sigma(ic, ii) = s
            end do
        end do
    end do
! !
! ! --- DETERMINATION DES MATRICE DE PASSAGE DES REPERES INTRINSEQUES
! ! --- AUX NOEUDS ET AUX POINTS D'INTEGRATION DE L'ELEMENT
! ! --- AU REPERE UTILISATEUR :
! !     ---------------------
    call vdrepe(nomte, matevn, matevg)
! !
    do i = 1, nb2
        do j = 1, 8
            effgt(j, i) = 0.d0
        end do
    end do
    do ic = 1, nbcou
        j = (ic-1)*npge*nso+1
        zic = zmin+(ic-1)*hic
        call vdefgn(nomte, nb2, hic, zic, sigma(1, j), &
                    effgc)
        do isom = 1, nb2
            do icomp = 1, 8
                effgt(icomp, isom) = effgt(icomp, isom)+effgc(icomp, isom)
            end do
        end do
    end do
! !
! ! --- PASSAGE DU VECTEUR DES EFFORTS GENERALISES DEFINI AUX NOEUDS
! ! --- DE L'ELEMENT DU REPERE INTRINSEQUE AU REPERE UTILISATEUR :
! !     --------------------------------------------------------
!
    call vdefro(nb2, matevn, effgt, effg)
!
!
!
end subroutine
