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
subroutine mltpos(nbsn, parent, fils, frere, pile, &
                  lfront, seq, flag, estim, u, &
                  w, tab, liste)
! person_in_charge: olivier.boiteau at edf.fr
! aslint: disable=
    implicit none
#include "asterfort/blimax.h"
#include "asterfort/tri.h"
    integer(kind=8) :: nbsn, parent(*), fils(*), frere(*), pile(*), lfront(*)
    integer(kind=8) :: seq(*), estim
    integer(kind=8) :: flag(*)
    integer(kind=8) :: u(nbsn), w(nbsn), tab(nbsn), liste(nbsn)
!
    integer(kind=8) :: init, filsi, nd, iq, md, m, i, k, sni, itemp, lp, q1, q2, sn
!-----------------------------------------------------------------------
!     CALCUL DES TABLEAUX U ET W (VOIR NOTES RESP. DE  ASHCRAFT ET YANG)
!
    do sni = 1, nbsn
        u(sni) = (lfront(sni)*(lfront(sni)+1))/2
    end do
    do sni = 1, nbsn
        q1 = u(sni)
        sn = fils(sni)
        if (sn .eq. 0) then
            w(sni) = q1
        else
            m = 1
            sn = fils(sni)
            liste(1) = sn
            tab(1) = u(sn)
            sn = frere(sn)
!          DO WHILE (SN.NE.0)
120         continue
            if (sn .ne. 0) then
                m = m+1
                liste(m) = sn
                tab(m) = tab(m-1)+u(m)
                sn = frere(sn)
                goto 120
! FIN DO WHILE
            end if
            do k = 1, m
                tab(k) = tab(k)+w(liste(k))
            end do
            q2 = tab(blimax(m, tab, 1))
            do i = 1, m
                q1 = q1+u(liste(i))
            end do
            w(sni) = max(q1, q2)
        end if
    end do
!-----------------------------------------------------------------------
!      MODIFICATION DE FILS ET FRERE POUR MINIMISER LA PILE
    do sni = 1, nbsn
        sn = fils(sni)
        if (sn .ne. 0) then
            m = 1
            liste(m) = sn
            tab(m) = w(liste(m))-u(liste(m))
            sn = frere(sn)
!          DO WHILE (SN.NE.0)
160         continue
            if (sn .ne. 0) then
                m = m+1
                liste(m) = sn
                tab(m) = w(liste(m))-u(liste(m))
                sn = frere(sn)
                goto 160
! FIN DO WHILE
            end if
            call tri(tab, liste, 1, m)
            fils(sni) = liste(m)
            sn = fils(sni)
            k = m-1
!          DO WHILE (K.GE.1)
170         continue
            if (k .ge. 1) then
                frere(sn) = liste(k)
                sn = liste(k)
                k = k-1
                goto 170
! FIN DO WHILE
            end if
            frere(liste(1)) = 0
        end if
    end do
!-----------------------------------------------------------------------
!      CALCUL DE LA SEQUENCE D'EXECUTION
!
    do i = 1, nbsn
        flag(i) = 0
    end do
    iq = 0
    do init = 1, nbsn
        if (parent(init) .eq. 0) then
            lp = 0
            filsi = init
!          DO WHILE (FILSI.NE.0)
200         continue
            if (filsi .ne. 0) then
!             ND = FILSI
                lp = lp+1
                pile(lp) = filsi
                filsi = fils(filsi)
                goto 200
! FIN DO WHILE
            end if
!          DO WHILE (LP.GT.0)
210         continue
            if (lp .gt. 0) then
220             continue
                nd = pile(lp)
                md = fils(nd)
!            DO WHILE (MD.NE.0)
230             continue
                if (md .ne. 0) then
                    if (flag(md) .eq. 0) then
                        if (fils(md) .eq. 0) then
                            iq = iq+1
                            seq(iq) = md
                            flag(md) = 1
                        else
                            lp = lp+1
                            pile(lp) = md
                            goto 220
                        end if
                    end if
                    md = frere(md)
                    goto 230
! FIN DO WHILE
                end if
                iq = iq+1
                seq(iq) = nd
                flag(nd) = 1
                lp = lp-1
                goto 210
! FIN DO WHILE
            end if
        end if
    end do
!      ESTIMATION DE LA PILE
    estim = 1
    itemp = 1
    do i = 1, nbsn
        sni = seq(i)
        m = lfront(sni)
        if (fils(sni) .eq. 0) then
            pile(sni) = itemp
            itemp = itemp+(m*(m+1))/2
            estim = max(estim, itemp-1)
        else
            itemp = itemp+(m*(m+1))/2
            estim = max(estim, itemp-1)
            pile(sni) = pile(fils(sni))
            itemp = pile(fils(sni))+(m*(m+1))/2
        end if
    end do
end subroutine
