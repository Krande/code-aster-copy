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

subroutine te0220(option, nomte)
!
    use calcul_module, only: ca_jvcnom_, ca_nbcvrc_
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:
!                         CALCUL DE L'ENERGIE THERMIQUE A L'EQUILIBRE
!                         OPTION : 'ETHE_ELEM'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!
    integer(kind=8) :: icodre(2), kpg, spt
    character(len=8) :: nompar(ca_nbcvrc_+1), fami, poum, novrc
    character(len=16) :: nomres(2)
    character(len=32) :: phenom
    real(kind=8) :: valres(2), valpar(ca_nbcvrc_+1)
    real(kind=8) :: dfdx(9), dfdy(9), poids, flux, fluy, epot
    real(kind=8) :: angmas(3), fluglo(2), fluloc(2), p(2, 2)
    integer(kind=8) :: ndim, nno, nnos, npg, kp, j, itempe, itemp, iener
    integer(kind=8) :: ipoids, ivf, idfde, jgano, igeom, imate, iret, nbpar, ipar
    aster_logical :: aniso
!     ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PTEMPER', 'L', itempe)
    call jevech('PENERDR', 'E', iener)
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call tecach('ONO', 'PINSTR', 'L', iret, iad=itemp)
    if (itemp .eq. 0) then
        nbpar = 0
        nompar(1) = ' '
        valpar(1) = 0.d0
    else
        nbpar = 1
        nompar(1) = 'INST'
        valpar(1) = zr(itemp)
    end if
!
    do ipar = 1, ca_nbcvrc_
        novrc = zk8(ca_jvcnom_-1+ipar)
        nbpar = nbpar+1
        nompar(nbpar) = novrc
        call rcvarc(' ', nompar(nbpar), poum, fami, kpg, spt, valpar(nbpar), iret)
        ASSERT(iret .eq. 0)
    end do
!
    call rccoma(zi(imate), 'THER', 1, phenom, iret)
    if (phenom .eq. 'THER') then
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'THER', nbpar, nompar, [valpar], &
                    1, 'LAMBDA', valres, icodre, 1)
        aniso = .false.
    else if (phenom .eq. 'THER_ORTH') then
        nomres(1) = 'LAMBDA_L'
        nomres(2) = 'LAMBDA_T'
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'THER_ORTH', nbpar, nompar, [valpar], &
                    2, nomres, valres, icodre, 1)
        aniso = .true.
        call getElemOrientation(ndim, nno, igeom, angmas)
        p(1, 1) = cos(angmas(1))
        p(2, 1) = sin(angmas(1))
        p(1, 2) = -sin(angmas(1))
        p(2, 2) = cos(angmas(1))
    else
        call utmess('F', 'ELEMENTS2_68')
    end if
!
    epot = 0.d0
    do kp = 1, npg
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids, dfdx, dfdy)
        flux = 0.d0
        fluy = 0.d0
        do j = 1, nno
            flux = flux+zr(itempe+j-1)*dfdx(j)
            fluy = fluy+zr(itempe+j-1)*dfdy(j)
        end do
        if (.not. aniso) then
            fluglo(1) = valres(1)*flux
            fluglo(2) = valres(1)*fluy
        else
            fluglo(1) = flux
            fluglo(2) = fluy
!
            fluloc(1) = p(1, 1)*fluglo(1)+p(2, 1)*fluglo(2)
            fluloc(2) = p(1, 2)*fluglo(1)+p(2, 2)*fluglo(2)
!
            fluloc(1) = valres(1)*fluloc(1)
            fluloc(2) = valres(2)*fluloc(2)
!
            fluglo(1) = p(1, 1)*fluloc(1)+p(1, 2)*fluloc(2)
            fluglo(2) = p(2, 1)*fluloc(1)+p(2, 2)*fluloc(2)
        end if
!
        epot = epot-(flux*fluglo(1)+fluy*fluglo(2))*poids
    end do
    zr(iener) = epot/2.d0
!
end subroutine
