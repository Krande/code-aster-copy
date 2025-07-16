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
subroutine mltpas(nbnd, nbsn, supnd, xadj, adjncy, &
                  anc, nouv, seq, global, adress, &
                  nblign, lgsn, nbloc, ncbloc, lgbloc, &
                  diag, col, lmat, place)
! person_in_charge: olivier.boiteau at edf.fr
!
    implicit none
    integer(kind=8) :: nbnd, nbsn, nbloc, ncbloc(*), lgbloc(*)
    integer(kind=8) :: supnd(nbsn+1), diag(0:nbnd), seq(nbsn)
    integer(kind=8) :: col(*)
    integer(kind=8) :: xadj(nbnd+1), adjncy(*), lmat
    integer(kind=8) :: anc(nbnd), nouv(nbnd)
    integer(kind=4) :: global(*)
    integer(kind=8) :: adress(nbsn+1)
    integer(kind=8) :: nblign(nbsn), lgsn(nbsn)
!
!=========================================================
!     CALCUL DES ADRESSES DANS LA FACTORISEE DES TERMES INITIAUX
!     VERSION ASTER AVEC MATRICE INITIALE COMPACTE PAR LIGNES
!     ET DDL DE LAGRANGE
!     DANS CETTE VERSION LES ADRESSES DES TERMES INITIAUX
!     SONT RANGEES DANS COL, QUI NE SERT PLUS.
!     AUPARAVANT ON UTILISAIT UN TABLEAU ADINIT
!==========================================================
    integer(kind=8) :: place(nbnd)
    integer(kind=8) :: i, j, ndj, sni, andi, andj, code, haut
    integer(kind=8) :: ndi, lfac, depart, ad, isn, longb, ib, ic
    isn = 0
    longb = 0
    lmat = diag(nbnd)
    do ib = 1, nbloc
        lfac = longb
        do ic = 1, ncbloc(ib)
            isn = isn+1
            sni = seq(isn)
            do i = adress(sni), adress(sni+1)-1
                place(global(i)) = i-adress(sni)+1
            end do
            haut = nblign(sni)
            do i = 0, lgsn(sni)-1
                ndi = supnd(sni)+i
                andi = anc(ndi)
                depart = diag(andi-1)+1
                col(diag(andi)) = lfac+i*haut+i+1+nbnd
                do j = xadj(andi), xadj(andi+1)-1
                    andj = adjncy(j)
                    ndj = nouv(andj)
                    if (ndj .ge. ndi) then
                        if (andj .le. andi) then
                            do ad = depart, diag(andi)
                                if (col(ad) .eq. andj) goto 140
                            end do
                            goto 170
140                         continue
                            code = -1
                            depart = ad
                        else
                            do ad = diag(andj-1)+1, diag(andj)
                                if (col(ad) .eq. andi) goto 160
                            end do
                            goto 170
160                         continue
                            code = 1
                        end if
!     CHANGTDGEMV                  ADINIT(AD) = LFAC + PLACE(NDJ) - I
!     ADINIT(AD) = LFAC + I*HAUT + PLACE(NDJ)
                        col(ad) = lfac+i*haut+place(ndj)+nbnd
!
                        if (code .lt. 0) then
                            col(ad) = -col(ad)
                        end if
                    end if
170                 continue
                end do
!
            end do
            lfac = lfac+haut*lgsn(sni)
        end do
        longb = lgbloc(ib)+longb
    end do
    do i = 1, lmat
        if (col(i) .gt. nbnd) then
            col(i) = col(i)-nbnd
        else if (col(i) .lt. (-nbnd)) then
            col(i) = col(i)+nbnd
        else
        end if
    end do
end subroutine
