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
subroutine te0530(option, nomte)
!
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
    character(len=16) :: option, nomte
!
! ......................................................................
!  CALCUL DES VARIABLES DE COMMANDE UTILISEES DANS LES CALCULS
!  MECANIQUES
! ......................................................................
!
    real(kind=8) :: r1, rvid
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfde, jgano, ipg, iret
    integer(kind=8) :: jpvarc, ivrc
    integer(kind=8) :: nbvarc
    integer(kind=8) :: jtab(7), idx, nbsp, isp
    parameter(nbvarc=8)
    character(len=8) :: nomvrc(nbvarc)
!
! ---------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)

    call tecach('OOO', 'PVARC_R', 'E', iret, nval=7, itab=jtab)
    jpvarc = jtab(1)
    nbsp = jtab(7)
    ASSERT(npg .eq. jtab(3))
    rvid = r8vide()
!
!     VARC_R   = R    TEMP HYDR SECH IRRA CORR PTOT NEUT1 NEUT2
    nomvrc(1) = 'TEMP'
    nomvrc(2) = 'HYDR'
    nomvrc(3) = 'SECH'
    nomvrc(4) = 'IRRA'
    nomvrc(5) = 'CORR'
    nomvrc(6) = 'PTOT'
    nomvrc(7) = 'NEUT1'
    nomvrc(8) = 'NEUT2'
!
    do ipg = 1, npg
        do isp = 1, nbsp
!
            do ivrc = 1, nbvarc
                call rcvarc(' ', nomvrc(ivrc), '+', 'RIGI', ipg, isp, r1, iret)
                idx = jpvarc+nbvarc*((ipg-1)*nbsp+isp-1)+ivrc-1
                zr(idx) = merge(r1, rvid, iret .eq. 0)

            end do
!
        end do
    end do
!
end subroutine
