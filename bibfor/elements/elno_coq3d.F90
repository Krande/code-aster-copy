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
subroutine elno_coq3d(option, nomte, nb1, nb2, npgsr, &
                      npgsn, nso, nbcou, geom, cara, &
                      valpg, outno, lzr, matr, lgreen)
!
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/caurtg.h"
#include "asterfort/pk2cau.h"
#include "asterfort/utmess.h"
#include "asterfort/vdrepe.h"
#include "asterfort/vdsiro.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
    character(len=16) :: option, nomte
!     CALCUL DES OPTIONS DES ELEMENTS DE COQUE 3D
!     OPTIONS : EPSI_ELNO
!               SIEF_ELNO
!               SIGM_ELNO
!          -----------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ic, icmp, ii, ino
    integer(kind=8) :: inte, intsn, intsr, isp, j
    integer(kind=8) :: jj, k, k1, kpgs, l
    integer(kind=8) :: nbcou, ncmp, npge, npgt, nso
!
    real(kind=8) :: s, zero
!-----------------------------------------------------------------------
    parameter(npge=3)
    parameter(npgt=10)
!
    integer(kind=8) :: icou, nordo
    integer(kind=8) :: nb1, nb2, npgsr, npgsn
!
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: epais
    real(kind=8) :: matevn(2, 2, npgt), matevg(2, 2, npgt)
    real(kind=8) :: geom(*), cara(*), valpg(*), outno(*), lzr(*), matr(*)
    real(kind=8) :: matpg(6, 27*nbcou), matno(6, 12*nbcou), matgn(6, 12*nbcou)
    real(kind=8) :: pk2(6, 27*nbcou), matgnu(6, 12*nbcou), signo(6, 12*nbcou)
!
    aster_logical :: lgreen
!
! ----------------------------------------------------------------------
!
    zero = 0.0d0
!
    if (nbcou .le. 0) then
        call utmess('F', 'ELEMENTS_12')
    end if
!
    epais = cara(1)
!
    call vectan(nb1, nb2, geom, lzr, vecta, &
                vectn, vectpt)
!
    kpgs = 0
    do icou = 1, nbcou
        do inte = 1, npge
            do intsn = 1, npgsn
                kpgs = kpgs+1
                k1 = 6*((intsn-1)*npge*nbcou+(icou-1)*npge+inte-1)
                do i = 1, 6
                    matpg(i, kpgs) = valpg(k1+i)
                end do
            end do
        end do
    end do
!
    ncmp = 6
!
    if (lgreen) then
!
! --- AFFECTATION DES CONTRAINTES DE PIOLA-KIRCHHOFF DE
! --- SECONDE ESPECE :
!     --------------
        do i = 1, 6
            do j = 1, kpgs
                pk2(i, j) = matpg(i, j)
            end do
        end do
!
! --- TRANSFORMATION DES CONTRAINTES DE PIOLA-KIRCHHOFF DE
! --- SECONDE ESPECE PK2 EN CONTRAINTES DE CAUCHY :
!     -------------------------------------------
        call pk2cau(nomte, ncmp, pk2, matpg)
    end if
!
! ---  DETERMINATION DES REPERES  LOCAUX DE L'ELEMENT AUX POINTS
! ---  D'INTEGRATION ET STOCKAGE DE CES REPERES DANS LE VECTEUR .DESR
!      --------------------------------------------------------------
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
!
!
    do icou = 1, nbcou
        do ic = 1, ncmp
            do i = 1, npge*nso
                l = npge*npgsn*(i-1)
                s = 0.d0
                do j = 1, npge*npgsn
                    jj = (icou-1)*npge*npgsn+j
                    s = s+matr(l+j)*matpg(ic, jj)
                end do
                ii = (icou-1)*npge*nso+i
                matno(ic, ii) = s
            end do
        end do
    end do
!
!
! --- DETERMINATION DES MATRICE DE PASSAGE DES REPERES INTRINSEQUES
! --- AUX NOEUDS ET AUX POINTS D'INTEGRATION DE L'ELEMENT
! --- AU REPERE UTILISATEUR :
!     ---------------------
    call vdrepe(nomte, matevn, matevg)
!
! --- PASSAGE DU VECTEUR DES CONTRAINTES DEFINI AUX NOEUDS
! --- DE L'ELEMENT DU REPERE INTRINSEQUE AU REPERE UTILISATEUR :
!     --------------------------------------------------------
!
    do icou = 1, nbcou
        do nordo = -1, 1
!
            isp = npge*(icou-1)+nordo+2
!
            do i = 1, ncmp
                do j = 1, nso
                    jj = nso*(nordo+1)+nso*npge*(icou-1)+j
                    matgn(i, j) = matno(i, jj)
                end do
                if (nomte .eq. 'MEC3QU9H') then
                    matgn(i, 5) = (matgn(i, 1)+matgn(i, 2))/2.d0
                    matgn(i, 6) = (matgn(i, 2)+matgn(i, 3))/2.d0
                    matgn(i, 7) = (matgn(i, 3)+matgn(i, 4))/2.d0
                    matgn(i, 8) = (matgn(i, 4)+matgn(i, 1))/2.d0
                    matgn(i, 9) = (matgn(i, 1)+matgn(i, 2)+matgn(i, 3)+ &
                                   matgn(i, 4))/4.d0
                else if (nomte .eq. 'MEC3TR7H') then
                    matgn(i, 4) = (matgn(i, 1)+matgn(i, 2))/2.d0
                    matgn(i, 5) = (matgn(i, 2)+matgn(i, 3))/2.d0
                    matgn(i, 6) = (matgn(i, 3)+matgn(i, 1))/2.d0
                    matgn(i, 7) = (matgn(i, 1)+matgn(i, 2)+matgn(i, 3))/ &
                                  3.d0
                end if
            end do
!
            if (lgreen) then
                call vdsiro(nb2, 1, matevn, 'IU', 'N', &
                            matgn, matgnu)
                call caurtg(nomte, ncmp, matgnu, signo)
            else
                call vdsiro(nb2, 1, matevn, 'IU', 'N', &
                            matgn, signo)
            end if
!
            if (option .eq. 'EPSI_ELNO') then
                do icmp = 1, ncmp
                    do ino = 1, nb2
                        outno((ino-1)*ncmp*nbcou*npge+(isp-1)*ncmp+icmp) = matgn(icmp, ino)
                    end do
                end do
            else if ((option .eq. 'SIEF_ELNO') .or. ( &
                     option .eq. 'SIGM_ELNO')) then
                do icmp = 1, ncmp
                    do ino = 1, nb2
                        outno((ino-1)*ncmp*nbcou*npge+(isp-1)*ncmp+icmp) = signo(icmp, ino)
                    end do
                end do
            else
                ASSERT(.false.)
            end if
!
        end do
    end do
end subroutine
