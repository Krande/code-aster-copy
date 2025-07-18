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
subroutine rc36fu(nbsigr, nocc, situ, saltij, nommat, &
                  ug, factus)
    implicit none
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/limend.h"
#include "asterfort/rc36f0.h"
#include "asterfort/rc36f2.h"
#include "asterfort/rcvale.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbsigr, nocc(*), situ(*)
    real(kind=8) :: saltij(*), ug, factus(*)
    character(len=*) :: nommat
!
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_B3600
!
!     CALCUL DU FACTEUR D'USAGE
!
!     ------------------------------------------------------------------
    integer(kind=8) :: isk, isl, k, l, nk, nl, n0, i1, i1a4, ifm, niv, icompt
    real(kind=8) :: saltm, nadm(1), ukl, vale(2)
    aster_logical :: trouve, endur
    integer(kind=8) :: icodre(1)
    character(len=2) :: k2c, k2l
    character(len=8) :: kbid
!     ------------------------------------------------------------------
!
    call infniv(ifm, niv)
!
    if (niv .ge. 2) then
        write (ifm, *) 'MATRICE SALT INITIALE'
        write (ifm, 1012) (situ(2*(l-1)+1), situ(2*(l-1)+2), l=1, nbsigr)
        write (ifm, 1010) (nocc(2*(l-1)+1), nocc(2*(l-1)+2), l=1, nbsigr)
        do k = 1, nbsigr
            i1 = 4*nbsigr*(k-1)
            write (ifm, 1000) situ(2*(k-1)+1), nocc(2*(k-1)+1), (saltij( &
                                                   i1+4*(l-1)+1), saltij(i1+4*(l-1)+3), l=1, nbsigr)
            write (ifm, 1002) situ(2*(k-1)+2), nocc(2*(k-1)+2), (saltij( &
                                                   i1+4*(l-1)+2), saltij(i1+4*(l-1)+4), l=1, nbsigr)
        end do
    end if
!
    icompt = 0
    ug = 0.d0
!
10  continue
    saltm = 0.d0
    trouve = .false.
!
! --- RECHERCHE DU SALT MAXI
!
    call rc36f0(nbsigr, nocc, saltij, saltm, trouve, &
                isk, isl, i1a4, nk, nl)
!
    if (.not. trouve) goto 999
!
    n0 = min(nk, nl)
    call limend(nommat, saltm, 'WOHLER', kbid, endur)
    if (endur) then
        ukl = 0.d0
    else
        call rcvale(nommat, 'FATIGUE', 1, 'SIGM    ', [saltm], &
                    1, 'WOHLER  ', nadm(1), icodre(1), 2)
        if (nadm(1) .lt. 0) then
            vale(1) = saltm
            vale(2) = nadm(1)
            call utmess('A', 'POSTRCCM_32', nr=2, valr=vale)
        end if
        ukl = dble(n0)/nadm(1)
    end if
!
    if (icompt .le. 49) then
        icompt = icompt+1
        factus(4*(icompt-1)+1) = i1a4
        factus(4*(icompt-1)+2) = situ(2*(isk-1)+1)
        factus(4*(icompt-1)+3) = situ(2*(isl-1)+1)
        factus(4*(icompt-1)+4) = ukl
    end if
!
    if (niv .ge. 2) then
        if (i1a4 .eq. 1 .or. i1a4 .eq. 3) then
            k2l = '_A'
        else
            k2l = '_B'
        end if
        if (i1a4 .eq. 1 .or. i1a4 .eq. 2) then
            k2c = '_A'
        else
            k2c = '_B'
        end if
        write (ifm, 1040) '=> SALT MAXI = ', saltm, situ(2*(isk-1)+1), &
            k2l, situ(2*(isl-1)+1), k2c
        write (ifm, 1030) '          N0 = ', n0
        write (ifm, 1020) '        nadm(1) = ', nadm(1)
        write (ifm, 1020) '         UKL = ', ukl
    end if
!
! --- MISE A ZERO DES LIGNES ET COLONNES DE LA MATRICE SALT SUIVANT
!     LE NOMBRE D'OCCURENCE EGAL A ZERO
!
    call rc36f2(nbsigr, nocc, saltij, i1a4, isk, &
                isl, nk, nl, n0)
!
    if (niv .ge. 2) then
        write (ifm, *) 'MATRICE SALT MODIFIEE'
        write (ifm, 1012) (situ(2*(l-1)+1), situ(2*(l-1)+2), l=1, nbsigr)
        write (ifm, 1010) (nocc(2*(l-1)+1), nocc(2*(l-1)+2), l=1, nbsigr)
        do k = 1, nbsigr
            i1 = 4*nbsigr*(k-1)
            write (ifm, 1000) situ(2*(k-1)+1), nocc(2*(k-1)+1), (saltij( &
                                                   i1+4*(l-1)+1), saltij(i1+4*(l-1)+3), l=1, nbsigr)
            write (ifm, 1002) situ(2*(k-1)+2), nocc(2*(k-1)+2), (saltij( &
                                                   i1+4*(l-1)+2), saltij(i1+4*(l-1)+4), l=1, nbsigr)
        end do
    end if
!
    ug = ug+ukl
    goto 10
!
999 continue
!
1000 format(1p, i7, '_A', i9, '|', 40(e9.2, 1x, e9.2, '|'))
1002 format(1p, i7, '_B', i9, '|', 40(e9.2, 1x, e9.2, '|'))
1010 format(1p, 9x, 'NB_OCCUR ', '|', 40(i9, 1x, i9, '|'))
1012 format(1p, 9x, 'SITUATION', '|', 40(i7, '_A', 1x, i7, '_B|'))
1040 format(1p, a15, e12.5, ', LIGNE:', i4, a2, ', COLONNE:', i4, a2)
! 1040 FORMAT(1P,A15,E12.5,' LIGNE ',I4,' COLONNE ',I4,' INDICE ',I4)
1030 format(1p, a15, i12)
1020 format(1p, a15, e12.5)
!
end subroutine
