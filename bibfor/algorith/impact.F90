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
subroutine impact(nmtab, nbpt, fn, vn, wk3, &
                  offset, t, elapse, nbchoc, fnmaxa, &
                  fnmmoy, fnmety, npari, lpari, valek)
    implicit none
#include "asterfort/tbajli.h"
    integer(kind=8) :: nbpt, nbchoc, npari
    real(kind=8) :: fn(*), t(*), vn(*), offset, elapse, wk3(*), fnmaxa, fnmety
    real(kind=8) :: fnmmoy
    character(len=16) :: lpari(*)
    character(len=24) :: valek(*)
    character(len=*) :: nmtab
!        COMPTAGE DES CHOCS AMV
!
! IN  : NBPT   : NB DE POINTS DU SIGNAL
! IN  : FN     : TABLEAU DU SIGNAL
! IN  : T      : TABLEAU DU TEMPS
! IN  : OFFSET : VALEUR DU SEUIL DE DETECTION D UN CHOC
! IN  : ELAPSE : TEMPS MINIMUM POUR VRAI FIN DE CHOC
! OUT : NBCHOC : NB DE CHOC GLOBAUX ( CRITERE ELAPSE )
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ipar(2), irebo, ichoc, idebut, ifin, nbpas, i, j, idech, nbrebo
    integer(kind=8) :: k
    real(kind=8) :: impuls, para(5), zero, fnmax, tchoc, vinit, dt, tfnmax
    real(kind=8) :: fnmmo2
    complex(kind=8) :: c16b
! ----------------------------------------------------------------------
!
    c16b = (0.d0, 0.d0)
    zero = 0.0d0
    nbchoc = 0
    nbrebo = 0
    fnmax = zero
    fnmaxa = zero
    fnmmoy = zero
    fnmety = zero
    tfnmax = zero
    impuls = zero
    tchoc = zero
    vinit = zero
    irebo = 0
    ichoc = 0
    idebut = 1
    ifin = 1
    dt = t(4)-t(3)
    nbpas = nint(elapse/dt)
!
    k = 0
    do i = 1, nbpt
!
        if (abs(fn(i)) .le. offset) then
!
!           SI SOUS LE SEUIL DE FORCE
!
            if (irebo .eq. 1) then
!
!              ET QUE ETAIT EN REBOND ALORS COMPTER REBOND
!
                nbrebo = nbrebo+1
            end if
!
            idech = 0
!
            do j = i, min(i+nbpas, nbpt)
!
!              EST CE QUE C'EST LA FIN D'UN CHOC GLOBAL
!
                if (abs(fn(j)) .gt. offset) idech = 1
            end do
!
            if (idech .eq. 0 .and. ichoc .eq. 1) then
!
!              OUI C'EST LA FIN D'UN CHOC GLOBAL
!
                ifin = i
                tchoc = t(ifin)-t(idebut)
                fnmmoy = fnmmoy+fnmax
!                    FNMMOY EST PROVISOIREMENT LE CUMUL DES FNMAX, ON
!                    DIVISE A LA FIN PAR NBCHOC POUR AVOIR LA MOYENNE
                fnmety = fnmety+fnmax*fnmax
                nbchoc = nbchoc+1
                ichoc = 0
                k = k+1
                wk3(k) = fnmax
                para(1) = tfnmax
                para(2) = fnmax
                para(3) = impuls
                para(4) = tchoc
                para(5) = vinit
                ipar(1) = nbchoc
                ipar(2) = nbrebo
                call tbajli(nmtab, npari, lpari, ipar, para, &
                            [c16b], valek, 0)
            end if
!
            irebo = 0
!
        else
!
            if (ichoc .eq. 0) then
!              DEBUT D'UN CHOC GLOBAL
                idebut = i
                vinit = vn(idebut)
                fnmax = zero
                impuls = zero
                nbrebo = 0
            end if
            if (i .eq. 1) then
                impuls = impuls+fn(i)*t(i)/2.d0
            else if (i .lt. nbpt) then
                j = i-1
                impuls = impuls+fn(i)*(t(i+1)-t(j))/2.d0
            else
                impuls = impuls+fn(i)*t(i)/2.d0
            end if
            if (fn(i) .ge. fnmax) then
                fnmax = fn(i)
                tfnmax = t(i)
            end if
            if (fnmax .ge. fnmaxa) fnmaxa = fnmax
            irebo = 1
            ichoc = 1
!
        end if
!
    end do
!
    if (nbchoc .ne. 0) then
!      ON PASSE PAR UNE VARIABLE INTERMEDIAIRE FNMMO2
!      POUR EVITER LES PROBLEMES DE PRECISION
        fnmmoy = fnmmoy/nbchoc
        fnmmo2 = fnmmoy*fnmmoy
        fnmety = sqrt(abs(fnmety/nbchoc-fnmmo2))
    else
        k = k+1
        wk3(k) = fnmax
        fnmmoy = zero
        fnmety = zero
        para(1) = tfnmax
        para(2) = fnmax
        para(3) = impuls
        para(4) = tchoc
        para(5) = vinit
        nbrebo = 0
        ipar(1) = nbchoc
        ipar(2) = nbrebo
        call tbajli(nmtab, npari, lpari, ipar, para, &
                    [c16b], valek, 0)
    end if
!
end subroutine
