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
subroutine rc36fs(nbsig1, noc1, sit1, nbsig2, noc2, &
                  sit2, saltij, ns, nscy, matse, &
                  mse, sn, nommat, c, k, &
                  cara, ug)
    implicit none
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/limend.h"
#include "asterfort/rc36sa.h"
#include "asterfort/rcvale.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbsig1, noc1(*), sit1(*), nbsig2, noc2(*), sit2(*), ns, nscy
    real(kind=8) :: saltij(*), matse(*), mse(*), sn(*), c(*), k(*), cara(*), ug
    character(len=8) :: nommat
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_B3600
!
!     CALCUL DU FACTEUR D'USAGE
!
!     ------------------------------------------------------------------
    integer(kind=8) :: is1, is2, is3, i, i1, i2, ind1, ind2, ifm, l, niv, ns2, icmp
    integer(kind=8) :: icomp
    real(kind=8) :: salt, saltm, nadm(1), u1kl, u2kl, sp, snkl, saltkl, mij, sm
    real(kind=8) :: vale(2)
    aster_logical :: trouve, endur
    real(kind=8) :: typeke, spmeca, spther
    integer(kind=8) :: icodre(1)
    character(len=2) :: k2c, k2l
    character(len=8) :: kbid
!     ------------------------------------------------------------------
!
    call infniv(ifm, niv)
!
    if (niv .ge. 2) then
        write (ifm, *) 'MATRICE SALT INITIALE (SEISME)'
        write (ifm, 1012) (sit2(2*(l-1)+1), sit2(2*(l-1)+2), l=1, nbsig2)
        write (ifm, 1010) (noc2(2*(l-1)+1), noc2(2*(l-1)+2), l=1, nbsig2)
        do i = 1, nbsig1
            i1 = 4*nbsig2*(i-1)
            write (ifm, 1000) sit1(2*(i-1)+1), noc1(2*(i-1)+1), (saltij( &
                                                   i1+4*(l-1)+1), saltij(i1+4*(l-1)+3), l=1, nbsig2)
            write (ifm, 1002) sit1(2*(i-1)+2), noc1(2*(i-1)+2), (saltij( &
                                                   i1+4*(l-1)+2), saltij(i1+4*(l-1)+4), l=1, nbsig2)
        end do
    end if
!
    ug = 0.d0
    ns2 = ns/2
    icomp = 0
!
    mij = 0.d0
    do icmp = 1, 3
        mij = mij+mse(icmp)**2
    end do
    mij = sqrt(mij)
!
    sp = k(2)*c(2)*cara(2)*mij/4/cara(1)
!
10  continue
    saltm = 0.d0
    trouve = .false.
    icomp = icomp+1
    if (icomp .gt. ns2) goto 999
!
    do i1 = 1, nbsig1
!
        ind1 = 4*nbsig2*(i1-1)
!
        do i2 = 1, nbsig2
!
            ind2 = 4*(i2-1)
!
            do i = 1, 4
                salt = saltij(ind1+ind2+i)
                if (salt .gt. saltm) then
                    is1 = i1
                    is2 = i2
                    is3 = i
                    saltm = salt
                    trouve = .true.
                end if
            end do
!
        end do
!
    end do
!
    if (trouve) then
!
        call limend(nommat, saltm, 'WOHLER', kbid, endur)
        if (endur) then
            u1kl = 0.d0
        else
            call rcvale(nommat, 'FATIGUE', 1, 'SIGM    ', [saltm], &
                        1, 'WOHLER  ', nadm(1), icodre(1), 2)
            if (nadm(1) .lt. 0) then
                vale(1) = saltm
                vale(2) = nadm(1)
                call utmess('A', 'POSTRCCM_32', nr=2, valr=vale)
            end if
            u1kl = 1.d0/nadm(1)
        end if
!
        snkl = sn(nbsig2*(is1-1)+(is2-1))
        typeke = -1.d0
        spmeca = 0.d0
        spther = 0.d0
        call rc36sa(nommat, matse, matse, snkl, sp, &
                    typeke, spmeca, spther, saltkl, sm)
        call limend(nommat, saltkl, 'WOHLER', kbid, endur)
        if (endur) then
            u2kl = 0.d0
        else
            call rcvale(nommat, 'FATIGUE', 1, 'SIGM    ', [saltkl], &
                        1, 'WOHLER  ', nadm(1), icodre(1), 2)
            if (nadm(1) .lt. 0) then
                vale(1) = saltkl
                vale(2) = nadm(1)
                call utmess('A', 'POSTRCCM_32', nr=2, valr=vale)
            end if
            u2kl = dble(2*nscy-1)/nadm(1)
        end if
!
        if (niv .ge. 2) then
            if (is3 .eq. 1 .or. is3 .eq. 3) then
                k2l = '_A'
            else
                k2l = '_B'
            end if
            if (is3 .eq. 1 .or. is3 .eq. 2) then
                k2c = '_A'
            else
                k2c = '_B'
            end if
            write (ifm, 1040) '=> SALT MAXI = ', saltm, sit1(2*(is1-1)+1) &
                , k2l, sit2(2*(is2-1)+1), k2c
            write (ifm, 1020) '        U1KL = ', u1kl
            write (ifm, 1020) '        SNKL = ', snkl
            write (ifm, 1020) '          SP = ', sp
            write (ifm, 1020) '      SALTKL = ', saltkl
            write (ifm, 1020) '        U2KL = ', u2kl
        end if
!
        ind1 = 4*nbsig2*(is1-1)
        ind2 = 4*(is2-1)
        saltij(ind1+ind2+is3) = 0.d0
        ind1 = 4*nbsig2*(is2-1)
        ind2 = 4*(is1-1)
        if (is3 .eq. 2) then
            saltij(ind1+ind2+3) = 0.d0
        else if (is3 .eq. 3) then
            saltij(ind1+ind2+2) = 0.d0
        else
            saltij(ind1+ind2+is3) = 0.d0
        end if
!
        if (niv .ge. 2) then
            write (ifm, *) 'MATRICE SALT MODIFIEE (SEISME)'
            write (ifm, 1012) (sit2(2*(l-1)+1), sit2(2*(l-1)+2), l=1, &
                               nbsig2)
            write (ifm, 1010) (noc2(2*(l-1)+1), noc2(2*(l-1)+2), l=1, &
                               nbsig2)
            do i = 1, nbsig1
                i1 = 4*nbsig2*(i-1)
                write (ifm, 1000) sit1(2*(i-1)+1), noc1(2*(i-1)+1), &
                    (saltij(i1+4*(l-1)+1), saltij(i1+4*(l-1)+3), l=1, &
                     nbsig2)
                write (ifm, 1002) sit1(2*(i-1)+2), noc1(2*(i-1)+2), &
                    (saltij(i1+4*(l-1)+2), saltij(i1+4*(l-1)+4), l=1, &
                     nbsig2)
            end do
        end if
!
        ug = ug+u1kl+u2kl
        goto 10
!
    end if
!
999 continue
!
1000 format(1p, i7, '_A', i9, '|', 40(e9.2, 1x, e9.2, '|'))
1002 format(1p, i7, '_B', i9, '|', 40(e9.2, 1x, e9.2, '|'))
1010 format(1p, 9x, 'NB_OCCUR ', '|', 40(i9, 1x, i9, '|'))
1012 format(1p, 9x, 'SITUATION', '|', 40(i7, '_A', 1x, i7, '_B|'))
1040 format(1p, a15, e12.5, ', LIGNE:', i4, a2, ', COLONNE:', i4, a2)
1020 format(1p, a15, e12.5)
!
end subroutine
