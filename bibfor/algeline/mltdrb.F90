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
subroutine mltdrb(nbloc, ncbloc, decal, seq, nbsn, &
                  nbnd, supnd, adress, global, lgsn, &
                  factol, factou, x, temp, invp, &
                  perm, ad, trav, typsym, nbsm, &
                  s)
!
! aslint: disable=W1504
    use superv_module
    implicit none
#include "jeveux.h"
#include "asterc/llbloc.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mlfmlt.h"
#include "asterfort/mlfmul.h"
!
    integer(kind=8) :: nbsn, nbnd, nbloc, ncbloc(nbnd), decal(nbsn)
    integer(kind=4) :: global(*)
    integer(kind=8) :: seq(nbsn), supnd(nbsn+1), lgsn(nbsn)
    integer(kind=8) :: adress(nbsn+1), invp(nbnd), perm(nbnd), ad(nbnd)
    integer(kind=8) :: typsym, nbsm
    real(kind=8) :: temp(nbnd), x(nbnd, nbsm), trav(nbnd, nbsm), s(nbsm)
    character(len=24) :: factol, factou, factor
    integer(kind=8) :: ib, nc, isnd, long, l, i, ndj, p
!
    integer(kind=8) :: deb1
    integer(kind=8) :: sni, k, j, deb, ifac, ism
    integer(kind=8) :: seuil, tranch, nproc, larg
    integer(kind=8) :: opta, optb, nb
    nb = llbloc()
    call jemarq()
    optb = 1
    nproc = asthread_getmax()
    tranch = (nbsm+nproc-1)/nproc
    seuil = nproc-mod(tranch*nproc-nbsm, nproc)
    do ism = 1, nbsm
        do j = 1, nbnd
            temp(invp(j)) = x(j, ism)
        end do
        do j = 1, nbnd
            x(j, ism) = temp(j)
        end do
    end do
!
!     DESCENTE  L * Y = B
    isnd = 0
    do ib = 1, nbloc
        call jeveuo(jexnum(factol, ib), 'L', ifac)
        do nc = 1, ncbloc(ib)
            isnd = isnd+1
            sni = seq(isnd)
            long = adress(sni+1)-adress(sni)
            l = lgsn(sni)
            do ism = 1, nbsm
                k = 1
                do i = adress(sni), adress(sni+1)-1
                    trav(k, ism) = x(global(i), ism)
                    k = k+1
                end do
            end do
            ad(1) = decal(sni)
            ndj = supnd(sni)-1
            do j = 1, l-1
                ndj = ndj+1
!     CALCUL DU BLOC  DIAGONAL
                temp(ndj) = zr(ifac-1+ad(j))
                do ism = 1, nbsm
                    k = 1
                    do i = j+1, l
                        trav(i, ism) = trav(i, ism)-zr(ifac-1+ad(j)+k)*trav(j, ism)
                        k = k+1
                    end do
                end do
                ad(j+1) = ad(j)+long+1
                ad(j) = ad(j)+l-j+1
            end do
            ndj = ndj+1
!     RANGEMENT DU TERME DIAGONAL
            temp(ndj) = zr(ifac-1+ad(l))
            ad(l) = ad(l)+1
            if (long .gt. l) then
                p = l
                opta = 1
                do ism = 1, nproc
                    if (ism .gt. seuil) then
                        larg = tranch-1
                        deb = seuil*tranch+(ism-seuil-1)*larg+1
                    else
                        deb = (ism-1)*tranch+1
                        larg = tranch
                    end if
! APPEL AU PRODUIT PAR BLOCS
                    call mlfmul(trav(p+1, deb), zr(ifac+ad(1)-1), trav(1, deb), nbnd, long, &
                                p, larg, opta, optb, nb)
                end do
            end if
            do ism = 1, nbsm
                k = 1
                do i = adress(sni), adress(sni+1)-1
                    x(global(i), ism) = trav(k, ism)
                    k = k+1
                end do
            end do
        end do
        call jelibe(jexnum(factol, ib))
    end do
!
    if (typsym .ne. 0) then
        factor = factol
        deb1 = 1
    else
        factor = factou
!     ON DIVISE PAR LE TERME DIAGONAL DANS LA REMONTEE EN NON-SYMETRIQUE
        deb1 = nbnd+1
    end if
!=======================================================================
!     D * Z = Y
    do ism = 1, nbsm
        do j = deb1, nbnd
            x(j, ism) = x(j, ism)/temp(j)
        end do
    end do
!=======================================================================
!     REMONTEE  U * X = Z
    isnd = nbsn+1
    opta = 0
    do ib = nbloc, 1, -1
        call jeveuo(jexnum(factor, ib), 'L', ifac)
        do nc = 1, ncbloc(ib)
            isnd = isnd-1
            sni = seq(isnd)
            l = lgsn(sni)
            long = adress(sni+1)-adress(sni)
            deb = adress(sni)+lgsn(sni)
            do ism = 1, nbsm
                k = 1
                do i = adress(sni), adress(sni+1)-1
                    trav(k, ism) = x(global(i), ism)
                    k = k+1
                end do
            end do
            p = l
            larg = nbsm
            deb = 1
            if (long .gt. p) then
! APPEL AU PRODUIT PAR BLOCS
                call mlfmlt(trav(1, deb), zr(ifac-1+decal(sni)+p), trav(p+1, deb), nbnd, &
                            long, p, larg, opta, optb, &
                            nb)
            end if
!     PARTIE DIAGONALE
            ad(1) = decal(sni)
            do j = 1, l-1
                ad(j+1) = ad(j)+long+1
            end do
            do j = l, 1, -1
!
                do ism = 1, nbsm
                    k = 1
                    s(ism) = 0.d0
                    do i = j+1, l
                        s(ism) = s(ism)+zr(ifac-1+ad(j)+k)*trav(i, &
                                                                ism)
                        k = k+1
                    end do
                    trav(j, ism) = trav(j, ism)-s(ism)
                    if (typsym .eq. 0) then
                        trav(j, ism) = trav(j, ism)/zr(ifac-1+ad(j))
                    end if
                end do
            end do
            do ism = 1, nbsm
                k = 1
                do i = adress(sni), adress(sni+1)-1
                    x(global(i), ism) = trav(k, ism)
                    k = k+1
                end do
            end do
        end do
        call jelibe(jexnum(factor, ib))
    end do
!     ON RANGE DANS SM  LA SOLUTION DANS LA NUMEROTATION INITIALE
    do ism = 1, nbsm
        do j = 1, nbnd
!
            temp(perm(j)) = x(j, ism)
        end do
!
        do j = 1, nbnd
            x(j, ism) = temp(j)
        end do
    end do
    call jedema()
end subroutine
