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

subroutine te0086(option, nomte)
!
! --------------------------------------------------------------------------------------------------
!
!    ELEMENT MECABL2
!       OPTION : 'RIGI_GEOM'
!
!       si pas de champe DEPL en entrée alors MATGEOM = 0
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
    character(len=16) :: option, nomte
!
#include "jeveux.h"
#include "asterfort/biline.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/matvec.h"
#include "asterfort/rcvalb.h"
#include "asterfort/get_value_mode_local.h"
#include "asterfort/tecach.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)             :: icodre(2)
    real(kind=8)        :: valres(2)
    character(len=16)   :: nomres(2)
!
    integer(kind=8) :: nno, kp, ii, jj, imatuu, ipoids, ivf, igeom, imate, iforce
    integer(kind=8) :: ideplp, idfdk, imat, iyty, iret
    integer(kind=8) :: jgano, kk, ndim, nelyty, nnos, nbval, nordre, npg
!
    real(kind=8) :: aire, coef1, coef2, demi, etraction, ecompress, ecable, nx, jacobi
    real(kind=8) :: ytywpq(9), w(9)
    real(kind=8) :: r8bid
!
    real(kind=8)        :: valr(2)
    character(len=8)    :: valp(2)
!
! --------------------------------------------------------------------------------------------------
    demi = 0.5d0
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
    call jevete('&INEL.CABPOU.YTY', 'L', iyty)
!   3 efforts par noeud
    nordre = 3*nno
!
! --------------------------------------------------------------------------------------------------
!   Parametres en sortie
    call jevech('PMATUUR', 'E', imatuu)
    nbval = nordre*(nordre+1)/2
    zr(imatuu:imatuu-1+nbval) = 0.0d0
! --------------------------------------------------------------------------------------------------
!   Parametres en entree : si pas de DEPL MATGEOM=0
    call tecach('ONO', 'PDEPLPR', 'L', iret, iad=ideplp)
    if (iret .ne. 0) goto 999
!
! --------------------------------------------------------------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTRR', 'L', iforce)
    call jevech('PMATERC', 'L', imate)
!
    nomres(1) = 'E'
    nomres(2) = 'EC_SUR_E'
    r8bid = 0.0d0
    call rcvalb('RIGI', 1, 1, '+', zi(imate), ' ', 'ELAS', 0, '  ', [r8bid], &
                1, nomres, valres, icodre, 1)
    call rcvalb('RIGI', 1, 1, '+', zi(imate), ' ', 'CABLE', 0, '  ', [r8bid], &
                1, nomres(2), valres(2), icodre(2), 1)
    etraction = valres(1)
    ecompress = etraction*valres(2)
    ecable = etraction
!
    valp(1) = 'SECT'
    call get_value_mode_local('PCACABL', valp, valr, iret, nbpara_=1)
    aire = valr(1)
!
    imat = imatuu-1
    do ii = 1, 3*nno
        w(ii) = zr(ideplp-1+ii)
    end do
    do kp = 1, npg
        kk = (kp-1)*nordre*nordre
        jacobi = sqrt(biline(nordre, zr(igeom), zr(iyty+kk), zr(igeom)))
!
        nx = zr(iforce-1+kp)
!       Le cable a un module plus faible en compression qu'en traction
!       Le module de compression peut même être nul.
        if (nx .lt. 0.0d0) then
            ecable = ecompress
        end if
!
        coef1 = ecable*aire*zr(ipoids-1+kp)/jacobi**3
        coef2 = nx*zr(ipoids-1+kp)/jacobi
        call matvec(nordre, zr(iyty+kk), 2, zr(igeom), w, ytywpq)
        nelyty = iyty-1-nordre+kk
        do ii = 1, nordre
            nelyty = nelyty+nordre
            do jj = 1, ii
                imat = imat+1
                zr(imat) = zr(imat)+coef1*ytywpq(ii)*ytywpq(jj)+coef2*zr(nelyty+jj)
            end do
        end do
    end do
!
999 continue
end subroutine
