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
subroutine te0178(option, nomte)
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
    implicit none
!                          D'AMORTISSEMENT ACOUSTIQUE SUR DES ARETES
!                          D'ELEMENTS 2D
!                          OPTION : 'AMOR_ACOU'
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/rcvalb.h"
#include "asterfort/vff2dn.h"
!
    complex(kind=8) :: rhosz
    character(len=8) :: fami, poum
    character(len=16) :: option, nomte
    integer :: icodre(1)
    real(kind=8) :: poids, r, nx, ny, rho(1)
    integer :: nno, kp, npg, ipoids, ivf, idfde, igeom
    integer :: imattt, i, j, ij, l, li, lj
    integer :: imate, iimpe, kpg, spt
    aster_logical :: laxi
!
!
!-----------------------------------------------------------------------
    integer :: jgano, mater, ndim, nnos
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg,&
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    laxi = .false.
    if (lteatt('AXIS','OUI')) laxi = .true.
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PIMPEDC', 'L', iimpe)
    call jevech('PMATERC', 'L', imate)
    call jevech('PMATTTC', 'E', imattt)
!
    fami='FPG1'
    kpg=1
    spt=1
    poum='+'
    mater = zi(imate)
    call rcvalb(fami, kpg, spt, poum, mater,&
                ' ', 'FLUIDE', 0, ' ', [0.d0],&
                1, 'RHO', rho, icodre, 1)
!
    if (zc(iimpe) .ne. (0.d0,0.d0)) then
        rhosz = rho(1)/zc(iimpe)
    else
        goto 50
    endif
!
    do kp = 1, npg
        call vff2dn(ndim, nno, kp, ipoids, idfde,&
                    zr(igeom), nx, ny, poids)
        if (laxi) then
            r = 0.d0
            do i = 1, nno
                l = (kp-1)*nno + i
                r = r + zr(igeom+2*i-2)*zr(ivf+l-1)
            end do
            poids = poids*r
        endif
        ij = imattt - 1
        do i = 1, nno
            li = ivf + (kp-1)*nno + i - 1
            do j = 1, i
                lj = ivf + (kp-1)*nno + j - 1
                ij = ij + 1
                zc(ij) = zc(ij) + poids*rhosz*zr(li)*zr(lj)
            end do
        end do
    end do
 50 continue
end subroutine
