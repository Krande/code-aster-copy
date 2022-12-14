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
subroutine te0506(option, nomte)
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
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_FLUTNL'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
    character(len=8) :: coef
    real(kind=8) :: poids, r, tpg
    real(kind=8) :: alpha, alphap, nx, ny
    integer :: nno, nnos, jgano, ndim, kp, npg, i, k, itemps, itemp, itempi
    integer :: iflux
    integer :: ipoids, ivf, idfde, igeom
    integer :: iveres
    aster_logical :: laxi
!
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg,&
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    laxi = .false.
    if (lteatt('AXIS','OUI')) laxi = .true.
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PTEMPER', 'L', itemp)
    call jevech('PTEMPEI', 'L', itempi)
    call jevech('PFLUXNL', 'L', iflux)
    call jevech('PRESIDU', 'E', iveres)
!
    coef = zk8(iflux)
    if (coef(1:7) .eq. '&FOZERO') goto 40
!
!
    do kp = 1, npg
        k = (kp-1)*nno
        call vff2dn(ndim, nno, kp, ipoids, idfde,&
                    zr(igeom), nx, ny, poids)
        r = 0.d0
        tpg = 0.d0
        do i = 1, nno
            r = r + zr(igeom+2* (i-1))*zr(ivf+k+i-1)
            tpg = tpg + zr(itempi+i-1)*zr(ivf+k+i-1)
        end do
        call foderi(coef, tpg, alpha, alphap)
        if (laxi) poids = poids*r
!
!
        do i = 1, nno
            zr(iveres+i-1) = zr(iveres+i-1) + poids*zr(ivf+k+i-1)* ( alpha-alphap*tpg)
        end do
    end do
 40 continue
! FIN ------------------------------------------------------------------
end subroutine
