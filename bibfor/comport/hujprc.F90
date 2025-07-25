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

subroutine hujprc(kk, k, tin, vin, mater, &
                  yf, p, q, toud)
    implicit none
!  LOI DE HUJEUX: PROJECTION DANS LE PLAN DEVIATEUR K
!  POUR UN MECANISME CYCLIQUE
!  IN  KK       :  NUMERO D'ORDRE DU MECANISME
!      K        :  MECANISME K = 1 A 3
!      TIN( )   :  TENSEUR DES CONTRAINTES DE CAUCHY
!      VIN      :  VARIABLES INTERNES ASSOCIEES
!      MATER    :  COEFFICIENTS MATERIAU
!      YF       :  VECTEUR SOLUTION DONNANT LA VALEUR DE R ET EPSVP
!
!  OUT
!      P     :  PRESSION ISOTROPE 2D DANS LE PLAN K
!      Q     :  NORME DEVIATEUR CYCLIQUE K
!      TOUD  :  TENSEUR DEVIATOIRE CYCLIQUE DES CONTRAINTES
!  --------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/infniv.h"
#include "asterfort/tecael.h"
    integer(kind=8) :: ndt, ndi, i, j, k, kk, nmod
    parameter(nmod=18)
    integer(kind=8) :: ifm, niv, iadzi, iazk24
    real(kind=8) :: yf(nmod), d12, dd, deux, vin(*)
    real(kind=8) :: r, x(2), th(2), pa, ptrac
    real(kind=8) :: tin(6), tou(3), toud(3), p, pp, q
    real(kind=8) :: epsvp, beta, b, phi, pcref, pcr
    real(kind=8) :: m, un, degr, mater(22, 2), tole
    aster_logical :: debug
    character(len=8) :: nomail
!
    parameter(degr=0.0174532925199d0)
    parameter(tole=1.d-7)
!
    common/tdim/ndt, ndi
    common/meshuj/debug
!
    data d12, deux, un/0.5d0, 2.d0, 1.d0/
!
    call infniv(ifm, niv)
!
!
! ==================================================================
! --- VARIABLES INTERNES -------------------------------------------
! ==================================================================
    epsvp = yf(7)
    r = yf(kk+7)
    x(1) = vin(4*k+5)
    x(2) = vin(4*k+6)
    th(1) = vin(4*k+7)
    th(2) = vin(4*k+8)
!
! ==================================================================
! --- CARACTERISTIQUES MATERIAU ------------------------------------
! ==================================================================
!
    beta = mater(2, 2)
    b = mater(4, 2)
    phi = mater(5, 2)
    pcref = mater(7, 2)
    pa = mater(8, 2)
    pcr = pcref*exp(-beta*epsvp)
    ptrac = mater(21, 2)
    m = sin(degr*phi)
!
!
! ==================================================================
! ----------------- CONSTRUCTION DU DEVIATEUR DES CONTRAINTES ------
! ==================================================================
    j = 1
    do i = 1, ndi
        if (i .ne. k) then
            tou(j) = tin(i)
            j = j+1
        end if
    end do
!
    tou(3) = tin(ndt+1-k)
    dd = d12*(tou(1)-tou(2))
!
!
! ==================================================================
! ----------------- CONSTRUCTION DU DEVIATEUR CYCLIQUE -------------
! ==================================================================
    pp = d12*(tou(1)+tou(2))-ptrac
!
    tou(1) = dd
    tou(2) = -dd
!
    if ((pp/pa) .le. tole) then
!
        if (debug) then
            call tecael(iadzi, iazk24)
            nomail = zk24(iazk24-1+3) (1:8)
            write (ifm, '(10(A))')&
     &    'HUJPRC :: LOG(P/PA) NON DEFINI DANS LA MAILLE ', nomail
        end if
!
        q = 0.d0
        toud(1) = 0.d0
        toud(2) = 0.d0
        toud(3) = 0.d0
        goto 999
    end if
!
    toud(1) = tou(1)-(x(1)-r*th(1))*pp*(un-b*log(pp/pcr))*m
    toud(2) = -toud(1)
    toud(3) = tou(3)-(x(2)-r*th(2))*pp*(un-b*log(pp/pcr))*m
    q = toud(1)**deux+(toud(3)**deux)/deux
    q = sqrt(q)
    p = pp+ptrac
!
999 continue
end subroutine
