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
subroutine stapu2(nbobst, nbpt, nbpair, temps, fcho,&
                  vgli, dloc, coef, ang, wk1,&
                  wk2, wk3, wk4, wk5, wk6,&
                  idebut, nbloc, nbval, inoe, isupp,&
                  nbinst, tmpvu, pusurn, vustub, vusob,&
                  pus, pmoye, pourpu, poupre)
!
!     CALCUL DE LA PUISSANCE D'USURE AU SENS D'ARCHARD
!
! IN  : NBOBST : NB DE NOEUDS DE CHOC
! IN  : NBPT   : NB DE PAS DE TEMPS TEMPORELS ARCHIVES
! IN  : TEMPS  : INSTANTS DE CALCUL
! IN  : FCHO   : VECTEUR DES FORCES DE CHOC
! IN  : VGLI   : VECTEUR DES VITESSES DE GLISSEMENT
! IN  : NBLOC  : NB DE BLOCS POUR LE MOYENNAGE
! IN  : INOE   : NUMERO DE NOEUD TRAITE
! OUT : PUSURN : PUISSANCE D'USURE MOYENNEE
!-----------------------------------------------------------------------
! aslint: disable=W1504
    implicit none
#include "asterc/r8rddg.h"
#include "asterfort/pusur2.h"
#include "blas/dcopy.h"
    integer :: occure(100)
    integer :: nbinst, nbloc, nbobst, nbpair, nbpt, nbval
    real(kind=8) :: temps(*), fcho(*), vgli(*), wk1(*), wk2(*), wk3(*), wk4(*)
    real(kind=8) :: wk5(*), wk6(*), dloc(*), ang(*), coef(*), tmpvu(*)
    real(kind=8) :: vustub(nbpair, nbinst), vusob(nbpair, nbinst), pmoye
    real(kind=8) :: poupre(*), pourpu(*), pus(*)
    real(kind=8) :: zero, pusee
!     REAL*8   TMPV(10)
!-----------------------------------------------------------------------
    integer :: ibl, idebut, in, inoe, isupp, j, jdg, k
    integer :: noccur, ntot
    real(kind=8) :: pusurn, rad
!-----------------------------------------------------------------------
    rad = r8rddg()
    zero = 0.d0
    pusurn = 0.d0
    ntot = 0
    pmoye = zero
    call dcopy(nbpt, dloc(3* (inoe-1)+2), 3*nbobst, wk4, 1)
    call dcopy(nbpt, dloc(3* (inoe-1)+3), 3*nbobst, wk5, 1)
    do in = 1, nbpt
        if ((wk4(in).ne.zero) .or. (wk5(in).ne.zero)) then
            wk6(in) = rad*atan2(wk5(in),wk4(in))
        else
            wk6(in) = zero
        endif
    end do
!
    call dcopy(nbpt, fcho(3* (inoe-1)+1), 3*nbobst, wk1, 1)
    call dcopy(nbpt, vgli(3* (inoe-1)+2), 3*nbobst, wk2, 1)
    call dcopy(nbpt, vgli(3* (inoe-1)+3), 3*nbobst, wk3, 1)
    do jdg = 1, nbpair
        pus(jdg) = zero
        occure(jdg) = 0
        do ibl = 1, nbloc
            call pusur2(jdg, nbval, ang, wk1((ibl-1)*nbval+idebut), wk2((ibl-1)*nbval+idebut),&
                        wk3((ibl-1)*nbval+idebut), wk6((ibl-1)*nbval+idebut),&
                        temps((ibl-1)*nbval+idebut), pusee, noccur)
            pus(jdg) = pus(jdg) + pusee
            occure(jdg) = occure(jdg) + noccur
        end do
        ntot = ntot + occure(jdg)
    end do
!
    do jdg = 1, nbpair
!
        pus(jdg) = pus(jdg)/nbloc
        pmoye = pmoye + pus(jdg)*occure(jdg)
    end do
    pmoye = pmoye/ntot
    do jdg = 1, nbpair
        poupre(jdg) = occure(jdg)*1.d0/ntot
        pourpu(jdg) = (pus(jdg)*occure(jdg)*100.d0/ (ntot*pmoye))
        if (isupp .eq. 1) then
            do k = 1, nbinst
                vustub(jdg,k) = coef(jdg)*tmpvu(k)*pus(jdg)*poupre( jdg)
!         TMPV(K)=TMPVU(K)/365.25/24/3600
            end do
        else if (isupp.eq.2) then
            do k = 1, nbinst
                vusob(jdg,k) = coef(jdg)*tmpvu(k)*pus(jdg)*poupre(jdg)
            end do
        endif
    end do
!
    do j = 1, nbpair
        pusurn = pusurn + coef(j)*occure(j)*pus(j)
    end do
    pusurn = pusurn/ntot
end subroutine
