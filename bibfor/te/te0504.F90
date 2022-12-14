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
subroutine te0504(option, nomte)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/foderi.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/vff2dn.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'RIGI_THER_FLUTNL'
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    real(kind=8) :: poids, r, nx, ny, alpha, alpha0
    integer :: nno, kp, npg, ipoids, ivf, idfde, igeom
    integer :: imattt, i, j, ij, l, li, lj, iflux
    aster_logical :: laxi
!
!-----------------------------------------------------------------------
    integer :: itemp, itemps, jgano, ndim, nnos
    real(kind=8) :: tpgi
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg,&
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    laxi = .false.
    if (lteatt('AXIS','OUI')) laxi = .true.
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PFLUXNL', 'L', iflux)
    call jevech('PTEMPER', 'L', itemp)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PMATTTR', 'E', imattt)
!
    if (zk8(iflux) (1:7) .eq. '&FOZERO') goto 50
!
    do kp = 1, npg
        call vff2dn(ndim, nno, kp, ipoids, idfde,&
                    zr(igeom), nx, ny, poids)
        r = 0.d0
        tpgi = 0.d0
        do i = 1, nno
            l = (kp-1)*nno + i
            r = r + zr(igeom+2*i-2)*zr(ivf+l-1)
            tpgi = tpgi + zr(itemp+i-1)*zr(ivf+l-1)
        end do
        if (laxi) poids = poids*r
        call foderi(zk8(iflux), tpgi, alpha, alpha0)
        ij = imattt - 1
        do i = 1, nno
            li = ivf + (kp-1)*nno + i - 1
!
            do j = 1, i
                lj = ivf + (kp-1)*nno + j - 1
                ij = ij + 1
                zr(ij) = zr(ij) - poids*alpha0*zr(li)*zr(lj)
            end do
        end do
    end do
 50 continue
end subroutine
