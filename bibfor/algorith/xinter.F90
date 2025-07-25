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

subroutine xinter(ndim, ndime, elrefp, geom, lsn, ia, ib, &
                  im, pintt, pmitt, lsna, lsnb, lsnm, inref, inter)
    implicit none
!
#include "jeveux.h"
#include "asterc/r8gaem.h"
#include "asterfort/assert.h"
#include "asterfort/reeref.h"
#include "asterfort/reerel.h"
#include "asterfort/xelrex.h"
#include "asterfort/xveri0.h"
#include "asterfort/xnewto.h"
    character(len=8) :: elrefp
    integer(kind=8) :: ndim, ndime, ia, ib, im
    real(kind=8) :: lsn(*), geom(*), inter(3), inref(3), pintt(*), pmitt(*)
    real(kind=8) :: lsna, lsnb, lsnm
!
!                      TROUVER LE PT D'INTERSECTION ENTRE L'ARETE
!                      ET LA FISSURE
!
!     ENTREE
!
!     SORTIE
!
!     ----------------------------------------------------------------
!
    character(len=6) :: name
    real(kind=8) :: ksi(ndime), ptxx(3*ndime), x(81), ff(27)
    real(kind=8) :: epsmax, a, b, c, pta(ndim), ptb(ndim), newpt(ndim)
    real(kind=8) :: ptm(ndim), dekker(4*ndime)
    integer(kind=8) :: itemax, ibid, n(3), j, nno, iret, exit(2)
!
!---------------------------------------------------------------------
!     DEBUT
!---------------------------------------------------------------------
!
    itemax = 100
    epsmax = 1.d-9
    name = 'XINTER'
    n(1) = ia
    n(2) = ib
    n(3) = 0
!   COORDONNEES DANS L ELEMENT DE REFERENCE PARENT
    call xelrex(elrefp, nno, x)
    if (ia .lt. 1000) then
        do j = 1, ndime
            pta(j) = x(ndime*(ia-1)+j)
        end do
    else
        do j = 1, ndim
            newpt(j) = pintt(ndim*(ia-1001)+j)
        end do
        call reeref(elrefp, nno, geom, newpt, ndim, &
                    pta, ff)
    end if
!
    if (ib .lt. 1000) then
        do j = 1, ndime
            ptb(j) = x(ndime*(ib-1)+j)
        end do
    else
        do j = 1, ndim
            newpt(j) = pintt(ndim*(ib-1001)+j)
        end do
        call reeref(elrefp, nno, geom, newpt, ndim, &
                    ptb, ff)
    end if
!
    if (im .eq. 0) then
        do j = 1, ndime
            ptm(j) = (pta(j)+ptb(j))/2.d0
        end do
    elseif (im .lt. 2000) then
        do j = 1, ndime
            ptm(j) = x(ndime*(im-1)+j)
        end do
    else
        do j = 1, ndim
            newpt(j) = pmitt(ndim*(im-2001)+j)
        end do
        call reeref(elrefp, nno, geom, newpt, ndim, &
                    ptm, ff)
    end if
!  ON STOCKE LES COORDONEES DE REFERENCE DE A ET B DANS <ptxx>
    do j = 1, ndime
        ptxx(j) = pta(j)
        ptxx(j+ndime) = ptb(j)
        ptxx(j+2*ndime) = ptm(j)
    end do
!!!!!ATTENTION INITIALISATION DU NEWTON: INTERPOLATION LINEAIRE DE LSN
    ksi(:) = 0.d0
!   INITIALISATION DU NEWTON
    ASSERT(abs(lsna-lsnb) .gt. 1.d0/r8gaem())
    a = (lsna+lsnb-2*lsnm)/2.d0
    b = (lsnb-lsna)/2.d0
    c = lsnm
    ASSERT(b**2 .ge. (4*a*c))
    if (abs(a) .lt. 1.d-8) then
        ksi(1) = lsna/(lsna-lsnb)
    else
        ksi(1) = (-b-sqrt(b**2-4*a*c))/(2.d0*a)
        if (abs(ksi(1)) .gt. 1) ksi(1) = (-b+sqrt(b**2-4*a*c))/(2.d0*a)
        ASSERT(abs(ksi(1)) .le. 1)
        ksi(1) = (ksi(1)+1)/2.d0
    end if
    dekker(:) = 0.d0
    exit(1:2) = 0
    call xnewto(elrefp, name, n, &
                ndime, ptxx, ndim, geom, lsn, &
                ibid, ibid, itemax, &
                epsmax, ksi, exit, dekker)
!  FIN DE RECHERCHE SUR SEGMENT AB
    do j = 1, ndime
        inref(j) = 2.d0*(1.d0-ksi(1))*(5.d-1-ksi(1))*ptxx(j)+4.d0*ksi(1)*(1.d0-ksi(1))* &
                   ptxx(j+2*ndime)+2.d0*ksi(1)*(ksi(1)-5.d-1)*ptxx(j+ndime)
    end do
!
    call xveri0(ndime, elrefp, inref, iret)
    ASSERT(iret .eq. 0)
!
    call reerel(elrefp, nno, ndim, geom, inref, &
                inter)
!
!---------------------------------------------------------------------
!     FIN
!---------------------------------------------------------------------
end subroutine
