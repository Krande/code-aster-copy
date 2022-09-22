! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine nmedpi(spg, sdg, qg, d, npg,&
                  typmod, mate, up, ud, geom,&
                  nno, def)
!
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
    integer :: nno, npg, mate
    character(len=8) :: typmod(*)
    real(kind=8) :: spg(2), sdg(2), qg(2, 2), d(4, 2)
    real(kind=8) :: geom(2, nno)
    real(kind=8) :: up(8), ud(8), def(4, nno, 2)
!
!--------------------------------------------------------
!  CALCUL DE SPG, SDG et QG AU POINT DE GAUSS COURANT
!  POUR LE PILOTAGE.
!--------------------------------------------------------
!
    integer :: i, j, k, kl, n
    real(kind=8) :: valres(2)
    real(kind=8) :: long, e, nu, deuxmu, troisk, val
    real(kind=8) :: mtemp(2, 4), mats(2, 8), dsidep(6, 6)
    real(kind=8) :: xa, xb, ya, yb
    character(len=8) :: fami, poum
    character(len=16) :: nomres(2)
    integer :: icodre(2), kpg, spt
    aster_logical :: axi
!
    axi = typmod(1) .eq. 'AXIS'
!
    call r8inir(2, 0.d0, spg, 1)
    call r8inir(2, 0.d0, sdg, 1)
    call r8inir(4, 0.d0, qg, 1)
    call r8inir(36, 0.d0, dsidep, 1)
!
! SOIT A ET B LES MILIEUX DES COTES [14] ET [23]
! t TANGENT AU COTE [AB]
!
    xa = ( geom(1,1) + geom(1,4) ) / 2
    ya = ( geom(2,1) + geom(2,4) ) / 2
!
    xb = ( geom(1,2) + geom(1,3) ) / 2
    yb = ( geom(2,2) + geom(2,3) ) / 2
!
!     LONGUEUR DE L'ELEMENT : NORME DU COTE [AB]
    long = sqrt( (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) )
    if (axi) then
        long = long * (xa + xb)/2.d0
    endif
!
! CALCUL DE SG ET QG :
! --------------------
!
! RECUPERATION DES CARACTERISTIQUES
    nomres(1)='E'
    nomres(2)='NU'
    fami='FPG1'
    kpg=1
    spt=1
    poum='+'
!
    call rcvalb(fami, kpg, spt, poum, mate,&
                ' ', 'ELAS', 0, ' ', [0.d0],&
                2, nomres(1), valres(1), icodre(1), 2)
    e = valres(1)
    nu = valres(2)
!
    deuxmu = e/(1.d0+nu)
    troisk = e/(1.d0-2.d0*nu)
!
!
! CALCUL DE DSIDEP
    call r8inir(36, 0.d0, dsidep, 1)
    do k = 1, 3
        do j = 1, 3
            dsidep(k,j) = troisk/3.d0 - deuxmu/(3.d0)
        end do
    end do
    do k = 1, 4
        dsidep(k,k) = dsidep(k,k) + deuxmu
    end do
!
! CALCUL DE MATS = Dt E B = Dt*DSIDEP*DEF
    call r8inir(8, 0.d0, mtemp, 1)
    do i = 1, 2
        do j = 1, 4
            val=0.d0
            do k = 1, 4
                val = val + d(k,i)*dsidep(k,j)
            end do
            mtemp(i,j)=val
        end do
    end do
!
    call r8inir(16, 0.d0, mats, 1)
    do k = 1, 2
        do n = 1, nno
            do i = 1, 2
                kl=2*(n-1)+i
                do j = 1, 4
                    mats(k,kl) = mats(k,kl) + mtemp(k,j)*def(j,n,i)
                end do
            end do
        end do
    end do
!
!
! CALCUL DE SPG ET SDG :
!
    do i = 1, 2
        do kl = 1, 8
            spg(i) = spg(i) - mats(i,kl)*up(kl)/long
            sdg(i) = sdg(i) - mats(i,kl)*ud(kl)/long
        end do
    end do
!
!
! CALCUL DE QG = Dt E D (Dt*DSIDEP*D)
!
    do i = 1, 2
        do j = 1, 2
            val=0.d0
            do k = 1, 4
                val = val - mtemp(i,k)*d(k,j)/long
            end do
            qg(i,j)=val
        end do
    end do
!
!
end subroutine
